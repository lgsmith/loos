#include <iostream>
#include <loos.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

const string fullHelpMessage =
    // clang-format off
"SYNOPSIS \n"
" \n"
"This is a tool for writing the lengths of vectors between atoms or groups of \n"
"adjacent residues in primary sequence to standard out, as a time-series. Here \n"
"I'll call these vectors 'link-vectors' or 'links' in the chain. By default, \n"
"this tool looks for either CAs or Ps, and draws vectors between the neighboring\n"
" such vectors along the length of the chain. The options permit some more \n"
"complex selections for defining links. \n"
" \n"
"DESCRIPTION \n"
" \n"
"This tool defines the  groups between which vectors are being drawn in a \n"
"fashion similar to the 'ocf' tool. The expected output will be Nframes rows x N\n"
" link-vectors +1 columns, where the +1 emerges from the frame index occupying \n"
"the first column. The i+1 column will be the length of the ith link. The groups\n"
" are defined by splitting up the selection into separate atomic groups. If \n"
"these groups are multi-atomic, then they will be reduced to either their \n"
"centroid or center of mass. \n"
" \n"
"The '--selection' determines which part of the broader structure the analysis \n"
"is performed within. The '--link-atom-selection' determines which atom(s) to \n"
"use within that scope. Cutting the scope up using 'group-centroids' assumes \n"
"that the link-atom-selection will specify non-contiguous parts of the same \n"
"molecule, or discontinuous molecules, that may not correspond to a division \n"
"along residue lines. For example, the distances between each phosphorus and \n"
"C1\' on the same residue could be tracked in a nucleic acid in this way. \n"
"Cutting up the scope using 'residue-centroids' breaks the scope into a \n"
"collection of residues, then applies the link selection within each residue and\n"
" uses the centroids. If the '--center-of-mass' flag is thrown, centers of mass \n"
"will be used in place of centroids to reduce any multiatomic link selection to \n"
"a point so that the distance to the neighboring link can be computed. \n"
" \n"
"Note that infer-connectivity uses the coordinates of the first frame to \n"
"determine which atoms in the model as a whole are bonded. It does this with a \n"
"hard distance cutoff. If none is provided, or if what is provided is not a \n"
"positive number, whatever connectivity is present in the model file will be \n"
"used. This tool requires connectivity because of how the aforementioned \n"
"splitting schemes work. \n"
" \n"
"EXAMPLES \n"
" \n"
"The following most simple example will draw link-vectors between adjacent \n"
"residues in the scope that contain either a 'CA' or a 'P'. So if model.psf has \n"
"a protein only within it, the following will write all the CA distances to \n"
"file. \n"
" \n"
"link-distances model.psf traj.dcd > distances.dat \n"
" \n"
"If there are other groups in the system, but the first 120 residues are \n"
"protein, one could reduce the scope of the model to just that subset like so: \n"
" \n"
"link-distances -s 'resid < 121' model.psf traj.dcd > distances.dat \n"
" \n"
"In both cases, because of the default link selection string, these examples \n"
"could also work on a nucleic acid and report on neighboring inter-phosphorus \n"
"distances instead. \n"
" \n"
"POTENTIAL COMPLICATIONS \n"
" \n"
"The order in which the splitting operations are done is probably the greatest \n"
"possible source of complication. Doing simple things, like trying to track a \n"
"particular backbone atom (like CA) that is defined for each residue within the \n"
"scope, are the primary reason this tool was written. One could potentially do \n"
"something like use a regex to select every nitrogen (they'd be disjoint in most\n"
" biophysically relevant systems) then track the distances between each \n"
"consecutive pair of nitrogens across the entire structure. If one does that, \n"
"doing some work with model-select or PyLOOS to ensure that the subselection is \n"
"getting you what you want is likely a good idea. \n"
" \n"
"Also, in regard to inferring connectivity, this has been mentioned elsewhere \n"
"but choosing to use a model without connectivity and inferring it in this crude \n"
"way can produce surprising results if you're not pretty sure of what's \n"
"happening with your system. For example, if the _first_ frame of your \n"
"trajectory has extremely distorted bonds for some reason, you could introduce \n"
"artificial discontinuities by trying to infer connectivity using the flag \n"
"provided by this tool. Depending on how generous with the margin you are, and \n"
"what kind of distortions or hydrogen bonds the first structure contains, you \n"
"could also accidentally introduce connectivity where none should exist. \n"
;
// clang-format on

class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() {}
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("link-atom-selection,L", po::value<string>(&bond_atom_selection)->
      default_value("name == 'CA' || name == 'P'"),
      "Selection of atoms to compute the link length across.")
      ("group-centroids", po::bool_switch(&group_centroids)->default_value(false),
       "If thrown, split link-atom-selection by molecule and compute link-lengths between centroids.")
      ("residue-centroids", po::bool_switch(&residue_centroids)->default_value(false),
       "Split link-atom-selection by residue and compute link lengths between residue centroids.")
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
  cout << "# frame length0 length1 length2 ...\n";
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
