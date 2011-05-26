/*
  water-sides.cpp
  (c) 2008, 2009 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School


  Desc:
    Classifies a water as being on one side of the membrane or the
    other or inside the membrane (1 = upper, 0 = inside, -1 = lower).

*/


#include <loos.hpp>
#include <limits>
#include <boost/format.hpp>
#include <DensityTools.hpp>


using namespace std;
using namespace loos;


typedef std::pair<double,double> Range;
typedef Math::Matrix<int, Math::ColMajor> Matrix;

Range membrane(0.0, 0.0);
string model_name, traj_name, selection_string;

enum Location { UPPER = 1, MEMBRANE = 0, LOWER = -1 };

/// @cond TOOLS_INTERNAL

class WaterSidesOptions : public opts::OptionsPackage {
public:

  WaterSidesOptions() : lower_bounds(0.0), upper_bounds(0.0) { }

  void addHidden(po::options_description& options) {
    options.add_options()
      ("lower", po::value<double>(&lower_bounds), "Lower leaflet bounds")
      ("upper", po::value<double>(&upper_bounds), "Upper leaflet bounds");
  }

  void addPositional(po::positional_options_description& options) {
    options.add("lower", 1);
    options.add("upper", 1);
  }

  bool check(po::variables_map& map) {
    return(! (map.count("lower") && map.count("upper")) );
  }

  string help() const { return("membrane-lower-bounds membrane-upper-bounds"); }
  string print() const {
    ostringstream oss;
    oss << boost::format("lower=%f, upper=%f") % lower_bounds % upper_bounds;
    return(oss.str());
  }

  double lower_bounds, upper_bounds;
};

/// @endcond


Range parseRange(const string& s) {
  double a, b;
  int i = sscanf(s.c_str(), "%lf:%lf", &a, &b);
  if (i != 2) {
    cerr << "Parse error with " << s << endl;
    exit(-1);
  }
  return(Range(a,b));
}




int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions *basic_opts = new opts::BasicOptions;
  opts::BasicSelection *basic_selection = new opts::BasicSelection("name == 'OH2'");
  opts::TrajectoryWithFrameIndices *basic_traj = new opts::TrajectoryWithFrameIndices;
  WaterSidesOptions *my_opts = new WaterSidesOptions;

  opts::AggregateOptions options;
  options.add(basic_opts).add(basic_selection).add(basic_traj).add(my_opts);
  if (!options.parse(argc, argv)) {
    options.showHelp();
    exit(0);
  }

  AtomicGroup model = basic_traj->model;
  pTraj traj = basic_traj->trajectory;
  AtomicGroup subset = selectAtoms(model, basic_selection->selection);
  vector<uint> frames = basic_traj->frameList();

  uint n = subset.size();
  uint m = frames.size();
  
  Matrix M(m,n+1);
  for (uint j=0; j<m; ++j) {
    traj->readFrame(frames[j]);
    traj->updateGroupCoords(subset);
    M(j, 0) = frames[j];
    for (uint i=0; i<n; ++i) {
      GCoord c = subset[i]->coords();
      Location l;
      if (c[2] > membrane.second)
        l = UPPER;
      else if (c[2] >= membrane.first)
        l = MEMBRANE;
      else
        l = LOWER;
      M(j, i+1) = l;
    }
  }

  writeAsciiMatrix(cout, M, hdr);
}
