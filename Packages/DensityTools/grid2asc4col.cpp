/*
  grid2asc4col.cpp


  Converts a grid (with a number of types) into a four column text
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include <DensityGrid.hpp>


namespace opts = loos::OptionsFramework;
namespace po = boost::program_options;

using namespace std;
using namespace loos;
using namespace loos::DensityTools;


// @cond TOOLS_INTERNAL

enum GridType { CHAR, INT, FLOAT, DOUBLE };

GridType gtype;
double scaling;
double threshold;

string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "\tConvert a LOOS grid into an ASCII four column format: x   y   z   density\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool converts a LOOS density grid into an x,y,z,density four vector\n"
    "that can be readily scripted in numpy, perl, awk, etc.  By default, the grid is\n"
    "assumed to contain double-precision floating point data (i.e. what is normally written\n"
    "out by the various LOOS tools).  Different data types can be converted by specifying\n"
    "what the grid contains on the command-line.\n"
    "\tThis tool can also threshold the grid for you, using the --threshold x option.\n"
    "\nEXAMPLES\n"
    "\tgrid2asc4col <foo.grid >foo.asc\n"
    "This converts a typical LOOS grid into a four column tab-delimited\n"
    "list of realspace points and densities\n\n"
    "\tgrid2asc4col --type int <foo_id.grid >foo.asc\n"
    "This converts an int-grid (from blobid, for example) into a tab-delimited\n"
    "list of realspace points and densities\n";

  return(msg);
}



class ToolOptions : public opts::OptionsPackage {
public:

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("type", po::value<string>(&type)->default_value("double"), "Set the grid type (char, int, float, double)")
      ("scale", po::value<double>(&scaling)->default_value(1.0), "Scale the grid data")
      ("threshold", po::value<double>(&threshold)->default_value(-1.0), "Threshold output to the provided level (none by default)");
  }

  bool postConditions(po::variables_map& map) {
    if (type == "double")
      gtype = DOUBLE;
    else if (type == "float")
      gtype = FLOAT;
    else if (type == "int")
      gtype = INT;
    else if (type == "char")
      gtype = CHAR;
    else {
      cerr << "Error- unknown grid type " << type << endl;
      return(false);
    }
    return(true);
  }

  string help() const {
    return(" <foo.grid >foo.asc");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("type='%s',scale='%f'") % type % scaling;
    return(oss.str());
  }

private:
  string type;

};



template<typename T>
DensityGrid<double> scaleGrid(DensityGrid<T>& g, const double scale) {
  DensityGridpoint dims = g.gridDims();
  long k = dims[0] * dims[1] * dims[2];
  DensityGrid<double> out(g.minCoord(), g.maxCoord(), g.gridDims());

  for (long i = 0; i<k; i++)
    out(i) = g(i) * scale;

  out.metadata(g.metadata());
  return(out);
}


// @endcond

int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);
  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
  ToolOptions* topts = new ToolOptions();
  
  opts::AggregateOptions options;
  options.add(bopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  DensityGrid<double> edm;
  if (gtype == CHAR) {
    DensityGrid<char> grid;
    cin >> grid;
    edm = scaleGrid(grid, scaling);

  } else if (gtype == INT) {
    DensityGrid<int> grid;
    cin >> grid;
    edm = scaleGrid(grid, scaling);

  } else if (gtype == FLOAT) {
    DensityGrid<float> grid;
    cin >> grid;
    edm = scaleGrid(grid, scaling);

  } else if (gtype == DOUBLE) {
    DensityGrid<double> grid;
    cin >> grid;
    edm = scaleGrid(grid, scaling);

  } else {
    cerr << "ERROR- bad grid type internally...\n";
    exit(-2);
  }

  edm.addMetadata(header);
  GCoord min = edm.minCoord();
  GCoord max = edm.maxCoord();
  DensityGridpoint dim = edm.gridDims();
  cerr << "Read in a grid of size " << dim << endl;
  cerr << "Grid range is from " << min << " to " << max << endl;
  if (threshold < 0) {
    for (int k = 0; k < dim.z(); ++k)
      for (int j = 0; j < dim.y(); ++j)
        for (int i = 0; i < dim.x(); ++i){
          DensityGridpoint gp(k, j, i);
          GCoord wp = edm.gridToWorld(gp);
          cout << boost::format("%f\t%f\t%f\t%f\n") % wp[0] % wp[1] % wp[2] % edm(k, j, i);
        }
  } else{
    for (int k = 0; k < dim.z(); ++k)
      for (int j = 0; j < dim.y(); ++j)
        for (int i = 0; i < dim.x(); ++i)
        {
          if (edm(k,j,i) < threshold)
            continue;
          else{
              DensityGridpoint gp(k, j, i);
              GCoord wp = edm.gridToWorld(gp);
              cout << boost::format("%f\t%f\t%f\t%f\n") % wp[0] % wp[1] % wp[2] % edm(k, j, i);
          }
        }
  }
}
