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

#include "loos.hpp"
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/CXX11/TensorSymmetry>

using namespace std;
using namespace loos;
using namespace Eigen;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

const string fullHelpMsg = "XXX";

// binary log 2. from bit-twiddling hacks naive impl.
inline uint binlog2(uint x){
  uint twoLog = 0;
  while (x >>= 1) ++twoLog;
  return twoLog;
}
// binary 2^x
inline uint binexp2(uint x){
  return 1 << x;
}
// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("time-series,t", po::value<string>(&ts)->default_value(""),
       "If thrown, write timeseries of instantaneous correlation function to stdout.")
      ("gamma,g", po::value<double>(&w)->default_value(42.58),
      "Set the larmor frequency to arg.")
      ("field-strength,B", po::value<double>(&B)->default_value(14.1),
       "Set the value of the experimental field strength.")
      ("sampling-freq,f", po::value<double>(&f)->default_value(1),
       "time-spacing of samples from trajectory, in GHz (frames/ns). A frequency of zero is an error.")
      ;
  }
  // clang-format on
  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("ts=%s,w=%d,B=%d,f=%d") % ts % w % B % f;
    return (oss.str());
  }
  string ts;
  double w,B,f;
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
  options.add(bopts).add(sopts).add(mtopts).add(topts);

  // Parse the command-line.  If an error occurred, help will already
  // be displayed and it will return a FALSE value.
  if (!options.parse(argc, argv))
    exit(-1);

  // Pull the model from the options object (it will include coordinates)
  AtomicGroup model = mtopts->model;

  // Pull out the trajectory...
  pTraj traj = mtopts->trajectory;
  cout << "got through parsing\n";

  // Select the desired atoms to operate over...
  AtomicGroup nuclei = selectAtoms(model, sopts->selection);
  Eigen::array<int, 3> bcBack({nuclei.size(), 3, 1});
  Eigen::array<int, 3> bcMiddle({nuclei.size(), 1, 3});
  Eigen::array<int, 3> bcFront({1, nuclei.size(), 3});
  Eigen::array<int, 3> bcFlipMiddle({1,3,nuclei.size()});
  Eigen::array<int, 3> bcBig({nuclei.size(), nuclei.size(), 3});
  Eigen::array<int, 2> flip({1,0});
  Eigen::array<int, 3> bigFlip({0,2,1});
  Tensor<double, 2> coords(nuclei.size(), 3);
  Tensor<double, 3> invs(nuclei.size(), nuclei.size(), 3);
  // SGroup<Symmetry<0,1>> sym;

  // pad batch dimension with zeros so that autocorrelation isn't circularly 
  if (mtopts->mtraj.nframes() == 0){
    cout << "ERROR: can't work with zero frames.\n";
    exit(-1);
  }
  // Now iterate over all frames in the skipped & strided trajectory
  while (traj->readFrame()) {

    // Update the coordinates ONLY for the subset
    traj->updateGroupCoords(nuclei);
    // Get coords into a tensor (presuming unrolling below)
    for (auto i = 0; i < nuclei.size(); i++)
      for (auto j = 0; j < 3; j++)
        coords(i, j) = nuclei[i]->coords()[j];
    // all_zdists.chip(mtopts->mtraj.currentFrame(), 2) = zcoords.broadcast(bc) - zcoords.broadcast(bc).shuffle(flip);
    // compute interatomic vectors

    cout << "Frame: " << mtopts->mtraj.currentFrame() << "\n";
    cout << "coords: \n";
    cout << coords << endl;
    cout << "broadcast back: \n";
    cout << coords.broadcast(bcBack) << endl;
    cout << "broadcast middle: \n";
    cout << coords.broadcast(bcMiddle)<< endl;
    cout << "broadcast front: \n";
    cout << coords.broadcast(bcFront)<< endl;
    cout << "broadcast flipMiddle: \n";
    cout << coords.broadcast(bcFlipMiddle)<< endl;
    cout << "broadcast big: \n";
    cout << coords.broadcast(bcBig) << endl;
    cout << "broadcast big flip: \n";
    cout << coords.broadcast(bcBig).shuffle(bigFlip) << endl;
    cout << "Internuclear Vectors:\n";
    invs = coords.broadcast(bcBig).shuffle(bigFlip) - coords.broadcast(bcBig);
    cout << invs << endl;
    

  }

  // cout << all_zdists << endl;
  // for (auto i = 0; i < mtopts->mtraj.nframes(); i++)
  //   cout << "frame " << i << "\n" << all_zdists.chip(i, 2) << "\n";

}
