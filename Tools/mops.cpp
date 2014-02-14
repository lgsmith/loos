/*
  Simple molecular-order parameters for comparing CG to AA MD
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo and Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

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

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


typedef vector<string>      vString;
typedef vector<double>      vecDbl;
typedef vector<AtomicGroup> vGroup;
typedef vector<vGroup>      vvGroup;


const double minp = 1e-3;
const double maxp = 100;
ulong nplanar = 0;
ulong ntotal = 0;


// @cond TOOLS_INTERNAL

struct ToolOptions : public opts::OptionsPackage 
{
 
  void addGeneric(po::options_description& o) 
  {
    o.add_options()
      ("skip", po::value<uint>(&skip)->default_value(0), "Skip these frames at the start of each trajectory");
  }

  void addHidden(po::options_description& o) 
  {
    o.add_options()
      ("selection", po::value<string>(&selection), "Atoms to use")
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value< vector<string> >(&traj_names), "Trajectory filenames");
  }

  void addPositional(po::positional_options_description& o) 
  {
    o.add("selection", 1);
    o.add("model", 1);
    o.add("traj", -1);
  }
  

  bool check(po::variables_map& vm) 
  {
    return( selection.empty() || model_name.empty() || traj_names.empty() );
  }

  string help() const 
  {
    return("selection model trajectory [trajectory ...]");
  }
  

  string print() const 
  {
    ostringstream oss;
    oss << boost::format("skip=%d, selection='%s', model='%s', traj='%s'")
      % skip
      % selection
      % model_name
      % vectorAsStringWithCommas(traj_names);
    

    return(oss.str());
  }

  uint skip;
  string selection;
  string model_name;
  vector<string> traj_names;
};


// @endcond



double rowStats(const RealMatrix& M, const uint row, double *std) {
  double avg = 0.0;
  for (uint i=0; i<M.cols(); ++i)
    avg += M(row, i);
  avg /= M.cols();
  
  double var = 0.0;
  for (uint i=0; i<M.cols(); ++i) {
    double d = M(row, i) - avg;
    var += d*d;
  }
  var /= (M.cols() - 1);

  *std = sqrt(var/M.cols());
  return(avg);
}


ulong calculateSize(AtomicGroup& model, const vString& names) {
  ulong n = 0;
  for (vString::const_iterator i = names.begin(); i != names.end(); ++i) {
    pTraj traj = createTrajectory(*i, model);
    n += traj->nframes();
  }
  return(n);
}


// Note: hist is really an average binned on distance
void principalComponentsOrder(dTimeSeries& order_parameters,
                              const vGroup& residues) {

  for (vGroup::const_iterator i = residues.begin(); i != residues.end(); ++i) {
    AtomicGroup residue = i->copy();
    residue.centerAtOrigin();
    residue.mergeImage();
    vector<GCoord> axes = residue.principalAxes();
    bool planar = false;

    if (axes[3].z() < minp) {
      if (nplanar == 0) {
        PDB pdb = PDB::fromAtomicGroup(residue);
        cerr << "Warning- PCA magnitudes out of bounds " << axes[3] << endl;
        cerr << pdb;
      }
      planar = true;
      ++nplanar;
    }

    double order1 = 0.5 - 1.5 * axes[1].z() * axes[1].z();
    double order2 = 0.5 - 1.5 * axes[2].z() * axes[2].z();

    order_parameters.push_back(order1);
    ++ntotal;
    if (!planar) {
      order_parameters.push_back(order2);
      ++ntotal;
    }
  }
}



vGroup extractSelections(const AtomicGroup& model, const string& selection) {
  AtomicGroup subset = selectAtoms(model, selection);
  vGroup residues = subset.splitByUniqueSegid();

  if (residues.empty()) {
    cerr << boost::format("ERROR- could not split group using selection '%s'\n") % selection;
    exit(EXIT_FAILURE);
  }
  
  // Autodetect whether we should use segid or residue to split...
  if (residues[0].size() == subset.size()) {
    cerr << "WARNING- apparent GROMACS source data...switching to splitByResidue() mode\n";
    residues = subset.splitByResidue();
  }
  return(residues);
}



int main(int argc, char *argv[]) {
  

  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);


  uint skip = topts->skip;
  string selection = topts->selection;
  AtomicGroup model = createSystem(topts->model_name);
  vGroup subset = extractSelections(model, selection);
  vString traj_names = topts->traj_names;

  cout << "# " << hdr << endl;



  // Calculate size of matrix to use...
  ulong n = calculateSize(model, traj_names);
  n -= traj_names.size() * skip;


  dTimeSeries order;

  for (vString::const_iterator i = traj_names.begin(); i != traj_names.end(); ++i) {
    dTimeSeries suborder;

    pTraj traj = createTrajectory(*i, model);
    if (skip != 0)
      traj->readFrame(skip-1);
    
    while (traj->readFrame()) {
      traj->updateGroupCoords(model);
      principalComponentsOrder(suborder, subset);
    }
    order.push_back(suborder.average());
  }

  

  cout << "Avg = " << order.average() << endl;
  cout << "Std = " << order.stdev() << endl;
  cout << boost::format("OB Data = %d out of %d (%.2f%%)\n") % nplanar % ntotal % (nplanar * 100.0 / ntotal);

  exit(EXIT_SUCCESS);
}
