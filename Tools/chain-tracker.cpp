#include <iostream>
#include <loos.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

const string fullHelpMessage =
    // clang-format off
"XXX"
;
// clang-format on

class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() {}
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("bond-atom-selection,B", po::value<string>(&bond_atom_selection)->
      default_value("name == 'CA' || name == 'P'"),
      "Selection of atoms to compute the OCF across")
      ("group-centroids", po::bool_switch(&group_centroids)->default_value(false),
       "If thrown, split bond-atom-selection by molecule and compute BVs between centroids.")
      ("residue-centroids", po::bool_switch(&residue_centroids)->default_value(false),
       "Split bond-atom-selection by residue, then track centroids for bond-vectors.")
      ("center-of-mass,c", po::bool_switch(&com)->default_value(false),
       "Instead of using centroids, use centers of mass for groups/residues.")
      ("infer-connectivity", po::value<float>(&bondlength)->default_value(-1), 
       "Infer connectivity using provided distance for models lacking this. ALERT: "
       "uses hard distance cutoff on first frame of traj to infer connectivity. "
       "Only does this for values greater than zero.")
    ;
  }
  // clang-format on
  string print() const {
    ostringstream oss;
    oss << boost::format("bond_atom_selection=%s,group_centroids=%b,bondlength=%d,residue_centroids%b,com=%b") %
               bond_atom_selection % group_centroids % bondlength %
               residue_centroids % com;
    return (oss.str());
  }

  bool postConditions(po::variables_map &map) {
    if (group_centroids && residue_centroids) {
      cerr << "ERROR: --group-centroids and --residue-centroids flags are "
              "mutually exclusive.\n";
      return (false);
    } else if (com && !(group_centroids || residue_centroids)) {
      cerr << "ERROR: --center-of-mass must be used with --group-centroids or"
              "--residue-centroids.\n";
      return (false);
    } else
      return (true);
  }
  string bond_atom_selection;
  bool group_centroids;
  bool residue_centroids;
  bool com;
  float bondlength;
};

inline void print_centroid_dists(vector<AtomicGroup> &chain,
                                 ostream& outs) {
  for (uint i = 0; i < chain.size() - 1; i++)
    outs << (chain[i].centroid() - chain[i + 1].centroid()).length() << " ";
  outs << "\n";
}

inline void print_com_dists(vector<AtomicGroup> &chain,
                             ostream& outs) {
  for (uint i = 0; i < chain.size() - 1; i++) {
   outs << (chain[i].centerOfMass() - chain[i + 1].centerOfMass()).length() 
        << " ";
  }
  outs << "\n";
}

int main(int argc, char *argv[]) {

  // parse the command line options
  string hdr = invocationHeader(argc, argv);
  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage);
  opts::BasicSelection *sopts = new opts::BasicSelection("all");
  opts::MultiTrajOptions *mtopts = new opts::MultiTrajOptions;
  ToolOptions *topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mtopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << hdr << "\n";
  // establish system, and subsystems
  AtomicGroup model = mtopts->model;
  if (model.hasBonds()) {
  } else if (topts->bondlength > 0)
    model.findBonds(topts->bondlength);
  else
    throw(LOOSError(
        "Model does not appear to have chemical connectivity, and "
        "infer-connectivity has not been set to a positive value.\n"));
  AtomicGroup scope = selectAtoms(model, sopts->selection);
  pTraj traj = mtopts->trajectory;
  // Set up chain to track
  vector<AtomicGroup> chain;
  // Define pointer to the function that will update the bond sites -- will be
  // either com_bond_vectors or centroid_bond_vectors, depending on user input
  void (*length_printer)(vector<AtomicGroup> &, ostream &);
  // determine what points to use for the centers of each link in the chain
  // based on user input
  if (topts->com) {
    length_printer = &print_com_dists;
    if (topts->group_centroids) {
      chain = scope.splitByMolecule(topts->bond_atom_selection);
    } else if (topts->residue_centroids) {
      chain = selectAtoms(scope, topts->bond_atom_selection).splitByResidue();
    }
  } else {
    length_printer = &print_centroid_dists;
    if (topts->group_centroids) {
      chain = scope.splitByMolecule(topts->bond_atom_selection);
    } else if (topts->residue_centroids) {
      chain = selectAtoms(scope, topts->bond_atom_selection).splitByResidue();
    } else {
      chain = selectAtoms(scope, topts->bond_atom_selection).splitByMolecule();
    }
  }
  // loop over trajectory
  for (auto frame_index : mtopts->frameList()) {
    traj->readFrame(frame_index);
    traj->updateGroupCoords(scope);
    cout << frame_index << " ";
    length_printer(chain, cout);
  }
  exit(EXIT_SUCCESS);
}
