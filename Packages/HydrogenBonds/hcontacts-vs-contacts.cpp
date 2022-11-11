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

#include "hcore.hpp"

using namespace std;
using namespace loos;
using namespace HBonds;
namespace po = boost::program_options;
namespace opts = loos::OptionsFramework;

bool inter_bonds, intra_bonds;
hreal putative_threshold;
hreal contact_threshold;

hreal length_low, length_high, max_angle;
bool use_periodicity;
string donor_selection, acceptor_selection;
string model_name, traj_name;

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
      ("search", po::value<hreal>(&putative_threshold)->default_value(10.0), "Threshold for initial bond search")
      ("contact", po::value<hreal>(&contact_threshold)->default_value(-1), "If > 0, use for contact. Otherwise use value of 'bhi'")
      ("blow", po::value<hreal>(&length_low)->default_value(1.5), "Low cutoff for bond length")
      ("bhi", po::value<hreal>(&length_high)->default_value(3.0), "High cutoff for bond length")
      ("angle", po::value<hreal>(&max_angle)->default_value(30.0), "Max bond angle deviation from linear")
      ("periodic", po::value<bool>(&use_periodicity)->default_value(false), "Use periodic boundary")
      ("inter", po::value<bool>(&inter_bonds)->default_value(true), "Inter-molecular bonds")
      ("intra", po::value<bool>(&intra_bonds)->default_value(false), "Intra-molecular bonds");
  }
  // clang-format on
  void addHidden(po::options_description &o)
  {
    o.add_options()("donor", po::value<string>(&donor_selection), "donor selection")("acceptor", po::value<string>(&acceptor_selection), "acceptor selection");
  }

  void addPositional(po::positional_options_description &opts)
  {
    opts.add("donor", 1);
    opts.add("acceptor", 1);
  }

  bool check(po::variables_map &map)
  {
    if (!(inter_bonds || intra_bonds))
    {
      cerr << "Error- must specify at least some kind of bond (inter/intra) to calculate.\n";
      return (true);
    }
    return (false);
  }

  string help() const
  {
    return ("donor-selection acceptor-selection");
  }

  string print() const
  {
    ostringstream oss;
    oss << boost::format("search=%f,inter=%d,intra=%d,blow=%f,bhi=%f,angle=%f,periodic=%d,acceptor=\"%s\",donor=\"%s\"") % putative_threshold % inter_bonds % intra_bonds % length_low % length_high % max_angle % use_periodicity % acceptor_selection % donor_selection;

    return (oss.str());
  }
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

  if (use_periodicity && !traj->hasPeriodicBox())
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

  vGroup mols = model.splitByMolecule();

  // First, build list of pairs we will search for...
  vGroup raw_donors = splitSelection(mols, donor_selection);
  vGroup raw_acceptors = splitSelection(mols, acceptor_selection);

  vBond bond_list;

  for (uint j = 0; j < raw_donors.size(); ++j)
  {

    if (intra_bonds)
    {
      vBond bonds = findPotentialBonds(raw_donors[j], raw_acceptors[j], model);
      copy(bonds.begin(), bonds.end(), back_inserter(bond_list));
    }

    if (inter_bonds)
    {
      for (uint i = 0; i < raw_acceptors.size(); ++i)
      {
        if (j == i)
          continue;
        vBond bonds = findPotentialBonds(raw_donors[j], raw_acceptors[i], model);
        copy(bonds.begin(), bonds.end(), back_inserter(bond_list));
      }
    }
  }

  // Generate the metadata for the output...
  ostringstream oss;
  oss << hdr << endl;
  for (uint i = 0; i < bond_list.size(); ++i)
  {
    formatBond(oss, i + 1, bond_list[i]);
    if (i < bond_list.size() - 1)
      oss << endl;
  }

  loos::Math::Matrix<uint> M(traj->nframes() - mtopts->skip, bond_list.size() + 1);

  uint j = 0;
  for (auto frame_index : mtopts->frameList())
  {
    traj->readFrame(frame_index);
    traj->updateGroupCoords(model);

    M(j, 0) = j + mtopts->skip;
    for (uint i = 0; i < bond_list.size(); ++i)
      M(j, i + 1) = bond_list[i].first.hydrogenBond(bond_list[i].second);

    ++j;
  }

  writeAsciiMatrix(cout, M, oss.str());
}
