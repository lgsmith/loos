/*
  internuclear-vector-corr-projections
  Writes the z projection of all selected internuclear vectors out as a timeseries,
  Optionally finding each autocorrelation as well.
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

#include "loos.hpp"
#include <eigen3/Eigen/Dense>


using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


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
    return(oss.str());
  }
  string autocorrs;
  bool ts;

};
// @endcond


// ***EDIT***
void calculate(const AtomicGroup& structure) {
  // Do something here with an atom...
}




int main(int argc, char *argv[]) {
  
  // Store the invocation information for logging later
  string header = invocationHeader(argc, argv);
  
  // Build up the command-line options for this tool by instantiating
  // the appropriate OptionsPackage objects...

  // Basic options should be used by all tools.  It provides help,
  // verbosity, and the ability to read options from a config file
  opts::BasicOptions* bopts = new opts::BasicOptions;

  // This tool can operate on a subset of atoms.  The BasicSelection
  // object provides the "--selection" option.
  opts::BasicSelection* sopts = new opts::BasicSelection;

  // The BasicTrajectory object handles specifying a trajectory as
  // well as a "--skip" option that lets the tool skip the first
  // number of frames (i.e. equilibration).  It creates a pTraj object
  // that is already primed for reading...
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;

  // ***EDIT***
  // Tool-specific options can be included here...
  ToolOptions* topts = new ToolOptions;

  // ***EDIT***
  // All of the OptionsPackages are combined via the AggregateOptions
  // object.  First instantiate it, then add the desired
  // OptionsPackage objects.  The order is important.  We recommend
  // you progress from general (Basic and Selection) to more specific
  // (model) and finally the tool options.
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);

  // Parse the command-line.  If an error occurred, help will already
  // be displayed and it will return a FALSE value.
  if (!options.parse(argc, argv))
    exit(-1);

  // Pull the model from the options object (it will include coordinates)
  AtomicGroup model = tropts->model;
  
  // Pull out the trajectory...
  pTraj traj = tropts->trajectory;

  // Select the desired atoms to operate over...
  AtomicGroup subset = selectAtoms(model, sopts->selection);

  // Now iterate over all frames in the trajectory (excluding the skip
  // region)
  while (traj->readFrame()) {

    // Update the coordinates ONLY for the subset of atoms we're
    // interested in...
    traj->updateGroupCoords(subset);

    // ***EDIT***
    // Now calculate something with the AtomicGroup
    calculate(subset);
  }

  // ***EDIT***
  // Output results...

}
