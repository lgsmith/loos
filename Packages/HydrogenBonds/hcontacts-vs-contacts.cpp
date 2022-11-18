/*
  hcontacts.cpp

  Constructs a matrix representing time series for multiple for inter and/or intra-molecular hbonds
*/

/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2022, Louis G. Smith
  Department of Biochemistry and Biophysics
  School of Medicine, University of Pennsylvania

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <loos.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <fstream>

#include "hcore.hpp"

using namespace std;
using namespace loos;
using namespace HBonds;
namespace po = boost::program_options;
namespace opts = loos::OptionsFramework;

// @cond TOOLS_INTERNAL

// clang-format off
string fullHelpMessage = 
  "XXX";
// clang-format on

class ToolOptions : public opts::OptionsPackage
{
public:
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("contact", po::value<hreal>(&contact_threshold)->default_value(-1), "If > 0, use for contact. Otherwise use value of 'bhi'")
      ("blow", po::value<hreal>(&length_low)->default_value(1.5), "Low cutoff for bond length")
      ("bhi", po::value<hreal>(&length_high)->default_value(3.0), "High cutoff for bond length")
      ("angle", po::value<hreal>(&max_angle)->default_value(30.0), "Max bond angle deviation from linear")
      ("periodic", po::value<bool>(&use_periodicity)->default_value(false), "Use periodic boundary")
      ("polymer-donors", po::value<string>(&polymer_donor_selection)-> default_value(""), 
      "If provided, use instead of 'donor' to select the atoms in 'polymer'")
      ("polymer-acceptors", po::value<string>(&polymer_acceptor_selection)-> default_value(""), 
      "If provided, use instead of 'acceptor' to select the atoms in 'polymer'");
  }
  void addHidden(po::options_description &o)
  {
    o.add_options()
    ("donor", po::value<string>(&donor_selection), "donor selection")
    ("acceptor", po::value<string>(&acceptor_selection), "acceptor selection")
    ("polymer", po::value<string>(&polymer_selection), "polymer selection")
    ("solvent", po::value<string>(&solvent_selection), "solvent selection")
    ("out_prefix", po::value<string>(&out_prefix), "Prefix to use to write output");
  }
  // clang-format on

  void addPositional(po::positional_options_description &opts)
  {
    opts.add("donor", 1);
    opts.add("acceptor", 1);
    opts.add("polymer", 1);
    opts.add("solvent", 1);
    opts.add("out_prefix", 1);
  }

  string help() const
  {
    return ("donor-selection acceptor-selection polymer-selection solvent-selection output-prefix");
  }

  string print() const
  {
    ostringstream oss;
    oss << boost::format("contact=%f,blow=%f,bhi=%f,angle=%f,"
                         "periodic=%d,acceptor=\"%s\",donor=\"%s\",polymer=\"%s\",solvent=\"%s\",out=\"%s\"") %
               contact_threshold % length_low % length_high % max_angle % use_periodicity % acceptor_selection % donor_selection % polymer_selection % solvent_selection % out_prefix;

    return (oss.str());
  }
  // tool options data members
  bool inter_bonds, intra_bonds;

  hreal contact_threshold;

  hreal length_low, length_high, max_angle;
  bool use_periodicity;
  string donor_selection, acceptor_selection;
  string polymer_selection, solvent_selection;
  string polymer_donor_selection, polymer_acceptor_selection;
  string model_name, traj_name;
  string out_prefix;
};

// @endcond

int main(int argc, char *argv[])
{
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage);
  opts::MultiTrajOptions *mtopts = new opts::MultiTrajOptions;
  ToolOptions *topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(mtopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = mtopts->model;
  pTraj traj = mtopts->trajectory;

  if (topts->use_periodicity && !traj->hasPeriodicBox())
  {
    cerr << "Error- trajectory has no periodic box information\n";
    exit(-1);
  }

  // Get coords if required...
  if (!model.hasCoords())
  {
    traj->readFrame();
    traj->updateGroupCoords(model);
    if (mtopts->skip > 0)
      traj->readFrame(mtopts->skip - 1);
    else
      traj->rewind();
  }

  // set up for the calculations
  AtomicGroup polymer = selectAtoms(model, topts->polymer_selection);
  AtomicGroup solvent = selectAtoms(model, topts->solvent_selection);
  vGroup vPoly = polymer.splitByResidue();
  vGroup vSolv = solvent.splitByMolecule();
  // will fill now, but refill over the course of the loop.
  // Will be organized by residues, then by coordinates within each AG
  vector<vector<GCoord>> poly_crds;
  for (auto res_ix = 0; res_ix < vPoly.size(); ++res_ix)
  {
    poly_crds.emplace_back(vector<GCoord>(vPoly[res_ix].size(), GCoord(0, 0, 0)));
  }

  // output file streams
  ofstream out_contacts(topts->out_prefix + "-contacts.out");
  ofstream out_hbs(topts->out_prefix + "-hbcounts.out");

  // contact threshold squared
  const hreal cthresh = topts->contact_threshold;

  // Build list of HBpairs we will search for.
  vGroup poly_donors, poly_accept;
  if (topts->polymer_donor_selection.empty())
  {
    poly_donors = splitSelection(vPoly, topts->donor_selection);
  }
  else
  {
    poly_donors = splitSelection(vPoly, topts->polymer_donor_selection);
  }
  if (topts->polymer_acceptor_selection.empty())
  {
    poly_accept = splitSelection(vPoly, topts->acceptor_selection);
  }
  else
  {
    poly_accept = splitSelection(vPoly, topts->polymer_acceptor_selection);
  }
  vGroup solv_donors = splitSelectionKeepEmpties(vSolv, topts->donor_selection);
  vGroup solv_accept = splitSelectionKeepEmpties(vSolv, topts->acceptor_selection);

  // This will be a vector (water indexes) of vectors (residue indexes) of vectors of bonds.
  vvvBond sp_pothbs_bymol;

  for (auto j = 0; j < vSolv.size(); ++j)
  {
    vvBond wat_bonds;
    for (auto i = 0; i < vPoly.size(); ++i)
    {
      // find all potential hbs between this solvent molecule's acceptors and the residue's donors.
      vBond pDsA_potentials = findPotentialBonds(poly_donors[i], solv_accept[j], model);
      // find all potential hbs between this molecule's donors and the residue's acceptors;
      // construct directly in vector of vector of bonds.
      wat_bonds.emplace_back(findPotentialBonds(solv_donors[j], poly_accept[i], model));
      // Concatenate polymer-donor-solvent-acceptor bonds onto stored vector of polyer-acceptor-solvent-donor bonds.
      wat_bonds[i].insert(wat_bonds[i].end(), pDsA_potentials.begin(), pDsA_potentials.end());
    }
  }

  // Generate the metadata for the output...
  out_contacts << hdr << "\n";
  out_hbs << hdr << "\n";
  for (auto res : vPoly)
  {
    out_contacts << res[0]->resid();
    out_hbs << res[0]->resid();
  }

  // loos Matricies initialized to zero
  loos::Math::Matrix<uint> HydrogenBondCounts(mtopts->frameList().size(), sp_pothbs_bymol[0].size() + 1);
  loos::Math::Matrix<uint> Contacts(mtopts->frameList().size(), sp_pothbs_bymol[0].size() + 1);

  uint frame_count = 0;
  uint total_contacts;
  for (auto frame_index : mtopts->frameList())
  {
    traj->readFrame(frame_index);
    traj->updateGroupCoords(model);

    HydrogenBondCounts(frame_count, 0) = frame_index;
    Contacts(frame_count, 0) = frame_index;
    // update the polymer coord vector, for easier scanning. Inner loop is over poly residues.
    for (auto j = 0; j < vPoly.size(); ++j)
      for (auto i = 0; i < vPoly[j].size(); ++i)
        poly_crds[j][i] = vPoly[j][i]->coords();

    for (auto solv_index = 0; solv_index < solv_donors.size(); ++solv_index)
    {
      for (auto res_index = 0; res_index < poly_accept.size(); ++res_index)
      {
        total_contacts = vSolv[solv_index].totalContacts(cthresh, poly_crds[res_index]);
        Contacts(frame_count, res_index) = total_contacts;
        // if the residue makes contact, then check for hydrogen bonds.
        if (total_contacts > 0)
        {
          for (auto pothb : sp_pothbs_bymol[solv_index][res_index])
            HydrogenBondCounts(frame_index, res_index) += pothb.first.hydrogenBond(pothb.second);
        }
      }
    }
    ++frame_count;
  }

  writeAsciiMatrix(out_contacts, Contacts, "polymer-solvent contact counts per residue");
  writeAsciiMatrix(out_hbs, HydrogenBondCounts, "polymer-solvent hydrogen bond counts per residue");
}
