/*
  
  Cosine content for varying windows of a trajectory,
  based on:
    Hess, B.  "Convergence of sampling in protein simulations."
      Phys Rev E (2002) 65(3):031910

*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo
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


#include "bcomlib.hpp"

using namespace std;
using namespace loos;
using namespace Convergence;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// @cond TOOLS_INTERNAL

typedef vector<AtomicGroup>                               vGroup;



// Convenience structure for aggregating results
struct Datum {
  Datum(const double avg, const double var, const uint nblks) : avg_cosine(avg),
                                                                var_cosine(var),
                                                                nblocks(nblks) { }


  double avg_cosine;
  double var_cosine;
  uint nblocks;
};



string fullHelpMessage() {

  string s = 
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Calculate the cosine constent of a simulation\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\n"
    "Calculate the cosine content in the PCA of a simulation.\n"
    "This is done in a block averaged manner.  Smaller contiguous\n"
    "blocks of the trajectory are used for the calcualation.  This\n"
    "tool assesses convergence of a trajectory by calulating how\n"
    "\"cosine-like\" the principal components appear.  The technique\n"
    "was described by Hess who noted that random diffusion showed a \n"
    "close resemblance to a cosine in:\n"
    "\n"
    "Hess, B. \"Convergence of sampling in protein \n"
    "     simulations.\" Phys Rev E (2002) 65(3):031910\n"
    "\n"
    "As the blocks increase in size more of the trajectory is used\n"
    "in this calulation.  The output is formatted in 4 columns:\n"
    "\n"
    "\t   n     - current block size (nanoseconds)\n"
    "\tCos cont - avg cosine content in a block\n"
    "\tVariance - variance across all (N_blocks)\n"
    "\tN_blocks - number of blocks of a given length\n"
    "\n"
    "\n"
    //
    "EXAMPLES\n"
    "\n"
    "coscon -s 'name==\"CA\"' model.pdb traj.dcd\n"
    "\tCalculate the block averaged cosine content in the\n"
    "\ttrajectory using only CA's in the calculation.\n"
    "\n"
    "coscon --pc=1 -s 'name==\"CA\"' model.pdb traj.dcd\n"
    "\tSame as above, but only compute the cosine content\n"
    "\tof the first principal component.\n"
    "\n"
    "coscon --block 100:100:500 -s 'name==\"CA\"' model.pdb traj.dcd\n"
    "\tHere we specify the range over which the calculation\n"
    "\tis performed.  The smallest block size used is 100 frames.\n"
    "\tThe block size will then increase by 100 frames upto a max\n"
    "\tblock size of 500 frames (ie...100,200,300,400,500).\n"
    "\n"
    "\n"
    "SEE ALSO\n"
    "\n"
    "Packages/Convergence/qcoscon - \n"
    "\tPerform a quick cos content analysis on a simulation.\n"
    "\tSimilar to coscon, but only performs the analysis on\n"
    "\tthe full length simulation.\n"
    "\n"
    "Packages/Convergence/rsv-coscon - \n"
    "\tCalculate the cos content of the RSVs from a simulation\n"
    "\tPCA.\n"
    "\t\n"
    "\n"
    "Tools/svd - \n"
    "\tCompute the principal components via the SVD.\n"
    "\tThis results in several matrix files including\n"
    "\tthe RSVs used as input to the current tool. \n"
    "\tThe file [prefix]_V.asc contains the RSV matrix.\n"
    "\n"
    "\n";
  
  
  return(s);
  
}

// @endcond

// Configuration

const uint nsteps = 50;



// Global options
bool local_average;
vector<uint> blocksizes;
string model_name, traj_name, selection;
uint principal_component;


// @cond TOOLS_INTERAL
class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("pc", po::value<uint>(&principal_component)->default_value(0), "Which principal component to use")
      ("blocks", po::value<string>(&blocks_spec), "Block sizes (MATLAB style range)")
      ("local", po::value<bool>(&local_average)->default_value(true), "Use local avg in block PCA rather than global");

  }

  bool postConditions(po::variables_map& vm) {
    if (!blocks_spec.empty())
      blocksizes = parseRangeList<uint>(blocks_spec);
    for (vector<uint>::iterator i = blocksizes.begin(); i != blocksizes.end(); ++i)
      if (*i == 0) {
        cerr << "Error- a block-size must be > 0\n";
        return(false);
      }

    return(true);
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("blocks='%s', local=%d, pc=%d")
      % blocks_spec
      % local_average
      % principal_component;
    return(oss.str());
  }

  string blocks_spec;
};
// @endcond







vGroup subgroup(const vGroup& A, const uint a, const uint b) {
  vGroup B;

  for (uint i=a; i<b; ++i)
    B.push_back(A[i]);

  return(B);
}



// Breaks the ensemble up into blocks and computes the RSV for each
// block and the statistics for the cosine content...

template<class ExtractPolicy>
Datum blocker(const uint pc, vGroup& ensemble, const uint blocksize, ExtractPolicy& policy) {


  TimeSeries<double> cosines;

  for (uint i=0; i<ensemble.size() - blocksize; i += blocksize) {
    vGroup subset = subgroup(ensemble, i, i+blocksize);
    RealMatrix V = rsv(subset, policy);

    double val = cosineContent(V, pc);
    cosines.push_back(val);
  }

  return( Datum(cosines.average(), cosines.variance(), cosines.size()) );
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;
  
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << hdr << endl;
  cout << "# " << vectorAsStringWithCommas<string>(options.print()) << endl;

  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;


  if (blocksizes.empty()) {
    uint n = traj->nframes();
    uint step = n / nsteps;
    if (step < 1)
      step = 1;
    cout << "# Auto block-sizes - " << step << ":" << step << ":" << n-1 << endl;

    for (uint i = step; i < n; i += step)
      blocksizes.push_back(i);
  }


  AtomicGroup subset = selectAtoms(model, sopts->selection);


  vector<AtomicGroup> ensemble;
  readTrajectory(ensemble, subset, traj);
 
  // First, read in and align trajectory
  boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(ensemble);
  AtomicGroup avg = averageStructure(ensemble);
  NoAlignPolicy policy(avg, local_average);


  // Now iterate over all requested block sizes

  // Provide user-feedback since this can be a slow computation
  PercentProgress watcher;
  ProgressCounter<PercentTrigger, EstimatingCounter> slayer(PercentTrigger(0.1), EstimatingCounter(blocksizes.size()));
  slayer.attach(&watcher);
  slayer.start();

  for (vector<uint>::iterator i = blocksizes.begin(); i != blocksizes.end(); ++i) {
    Datum result = blocker(principal_component, ensemble, *i, policy);
    cout << *i << "\t" << result.avg_cosine << "\t" << result.var_cosine << "\t" << result.nblocks << endl;
    slayer.update();
  }

  slayer.finish();

}
