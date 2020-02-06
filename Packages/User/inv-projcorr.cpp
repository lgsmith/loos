/*
  internuclear-vector-corr-projections
  Writes the z projection of all selected internuclear vectors out as a
  timeseries, Optionally finding each autocorrelation as well.
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2020, Louis
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

#include "Packages/Clustering/Clustering.hpp"
#include "loos.hpp"
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace loos;
using namespace eigen;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

const string fullHelpMsg = "XXX";

// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("no-time-series,N", po::bool_switch(&ts)->default_value(false),
       "If thrown, suppresses writing timeseries to stdout.")
      ("auto-correlation,c", po::value<string>(&autocorrs)->default_value(""), 
      "If provided, write autocorrelations for each inv-projection to filename.");
  }
  // clang-format on
  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("autocorrs=%s,ts=%b") % autocorrs % ts;
    return (oss.str());
  }
  string autocorrs;
  bool ts;
};
// @endcond

int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);

  opts::BasicOptions *bopts = new opts::BasicOptions;
  opts::BasicSelection *sopts = new opts::BasicSelection;
  opts::MultiTrajOptions *mtopts = new opts::MultiTrajOptions;
  ToolOptions *topts = new ToolOptions;

  // combine options
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);

  // Parse the command-line.  If an error occurred, help will already
  // be displayed and it will return a FALSE value.
  if (!options.parse(argc, argv))
    exit(-1);

  // Pull the model from the options object (it will include coordinates)
  AtomicGroup model = mtopts->model;

  // Pull out the trajectory...
  pTraj traj = mtopts->trajectory;

  // Select the desired atoms to operate over...
  AtomicGroup nuclei = selectAtoms(model, sopts->selection);
  VectorXd zcoords(nuclei.size());
  MatrixXd zdists(nuclei.size(), nuclei.size());
  vector<matrixX
  // Now iterate over all frames in the skipped & strided trajectory
  while (traj->readFrame()) {

    // Update the coordinates ONLY for the subset
    traj->updateGroupCoords(nuclei);
    // pick out zcoords
    for (auto i=0; i < nuclei.size(); i++){
      zcoords(i) = nuclei[i]->coords();
    }
    zdists = Clustering::pairwiseDists(zcoords);

  }

}
