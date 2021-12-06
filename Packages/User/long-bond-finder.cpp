/*
  long-bond-finder.cpp

  (c) 2021 Louis Smith, Bowman Lab
           Department of Biochemistry & Biophysics
           Washington University School of Medicine

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011 Tod D. Romo
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
#include <iostream>
#include <vector>
#include <utility>
#include <unordered_set>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

// The following conditional prevents this class from appearing in the
// Doxygen documentation for LOOS:
//
// @cond TOOLS_INTERNAL
// clang-format off
class ToolOptions : public opts::OptionsPackage {
public:

  // Change these options to reflect what your tool needs
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("infer-connectivity", po::value<float>(&bondlength)->default_value(-1), 
       "Infer connectivity using provided distance for models lacking this. "
       "ALERT: uses provided value as hard distance cutoff on first frame of traj to infer connectivity. "
       "Only does this for values greater than zero.")
      ("max-bond,M", po::value<float>(&max_bond)->default_value(2.5),
       "Maximum permissible distance for plausible bond.")
      ("timeseries,t", po::value<string>(&timeseries)->default_value(""), 
       "Write bond-pairs in violation of cutoff per-frame to file name provided.");
  }

  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("bondlength=%f,max_bond=%f,timeseries=%s") % bondlength % max_bond % timeseries;
    return(oss.str());
  }

  float bondlength, max_bond;
  string timeseries;

};
// clang-format on
// @endcond
// ----------------------------------------------------------------


// clang-format off
const string msg = 
"XXX";
// clang-format on

std::vector<std::pair<int, int>> getBondList(AtomicGroup &g) {
    // should hash pairs in an unordered way, that is hash(pair(a,b)) == hash(pair(b,a))
    struct unordered_pair_hash
    {
      std::size_t operator () (std::pair<int, int> const &pair) const {
        std::size_t h0 = std::hash<int>()(pair.first);
        std::size_t h1 = std::hash<int>()(pair.second);
 
        return h0 ^ h1;
      }
    };

    std::unordered_set<std::pair<int, int>, unordered_pair_hash> bond_set;
    AtomicGroup::const_iterator ci;

    for (ci = g.begin(); ci != g.end(); ++ci){
      std::vector<int> bonds = (*ci)->getBonds();
      int id = (*ci)->id();
      for(auto b : bonds){
        bond_set.emplace(std::make_pair(id, b));
      }
    }
    std::vector<std::pair<int, int>> bond_list(bond_set.begin(), bond_set.end());
    return bond_list;
  }

int main(int argc, char *argv[]) {
  
  string header = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions(msg);
  opts::BasicSelection* sopts = new opts::BasicSelection("all");
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);

  if (!options.parse(argc, argv))
    exit(-1);

  // set up system for looping. Load coords from frame 0 into scope.
  AtomicGroup model = tropts->model;
  if (model.hasBonds()) {
  } else if (topts->bondlength > 0)
    model.findBonds(topts->bondlength);
  else
    throw(LOOSError(
        "Model does not appear to have chemical connectivity, and "
        "infer-connectivity has not been set to a positive value.\n"));
  AtomicGroup scope = selectAtoms(model, sopts->selection);
  pTraj traj = tropts->trajectory;
  traj->updateGroupCoords(model);
  vector<pair<int, int>> bond_list = getBondList(scope);
  double max_bond2 = topts->max_bond * topts->max_bond;

  // Operating in scanning mode; don't report anything except the presence of an unacceptable bond
  if (topts->timeseries.empty()){
    cout << "# " << header << "\n";
    for (auto frame_index : tropts->frameList()){
        traj->readFrame(frame_index);
        traj->updateGroupCoords(scope);
        for(auto b : bond_list){
          if (scope.getAtom(b.first)->coords().distance2(scope.getAtom(b.second)->coords()) > max_bond2){
            return EXIT_FAILURE;
          }
        }
      }
  } else {
    ofstream tsf(topts->timeseries);
    bool found_viol = false;
    tsf << "# " << header << "\n";
    for (auto frame_index : tropts->frameList()){
      traj->readFrame(frame_index);
      traj->updateGroupCoords(scope);
      for(auto b : bond_list){
        if (scope.getAtom(b.first)->coords().distance2(scope.getAtom(b.second)->coords()) > max_bond2){
          found_viol = true;
          tsf << frame_index << " " << b.first << " " << b.second << "\n";
        }
      }
    }
    if (found_viol)
      return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
