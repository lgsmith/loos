
#include <loos.hpp>
#include "boost/tuple/tuple.hpp"

using namespace std;
using namespace loos;
// using namespace ::boost::tuples;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

string fullHelpMessage(void)
{
  string msg =
      "\nNothing here yet\n";
  return (msg);
}

typedef tuple<pAtom, pAtom, double> AtomPair;

bool cmpAtomPairs(const AtomPair &a1, const AtomPair &a2)
{
  return (get<2>(a1) < get<2>(a2));
}

int main(int argc, char *argv[])
{

  cout << "# " << invocationHeader(argc, argv) << endl;
  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage());
  opts::MultiTrajOptions *mtopts = new opts::MultiTrajOptions;
  opts::BasicSelection *sopts = new opts::BasicSelection("all");
  opts::WeightsOptions *wopts = new opts::WeightsOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mtopts).add(wopts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup system = mtopts->model;
  pTraj traj = mtopts->trajectory;

  AtomicGroup molecule = selectAtoms(system, sopts->selection);
  ;
  AtomicGroup hydrogens = molecule.select(HydrogenSelector());
  cout << "# Number of hydrogens = " << hydrogens.size() << endl;

  // TODO: Add a check for connectivity

  vector<AtomPair> hydrogen_pairs;

  // Build up a list of the pairs to check
  for (uint i = 0; i < hydrogens.size() - 1; i++)
  {
    int bondedTo = hydrogens[i]->getBonds()[0]; // hydrogens only have 1 bond
    for (uint j = i + 1; j < hydrogens.size(); j++)
    {
      // skip hydrogens bound to the same heavy atom
      // TODO: other pairs we should skip
      if (bondedTo != hydrogens[j]->getBonds()[0])
      {
        hydrogen_pairs.push_back(
            AtomPair(hydrogens[i], hydrogens[j], 0.0));
      }
    }
  }
  cout << "# Number of pairs = " << hydrogen_pairs.size() << endl;
  // track frame count for ranking after traj loop.
  uint frame = 0;
  if (wopts->has_weights)
  {
    wopts->weights.add_traj(traj);

    while (traj->readFrame())
    {
      traj->updateGroupCoords(system);

      for (vector<AtomPair>::iterator p = hydrogen_pairs.begin();
           p != hydrogen_pairs.end();
           ++p)
      {
        pAtom a1 = get<0>(*p);
        pAtom a2 = get<1>(*p);

        double d2 = a1->coords().distance2(a2->coords(), system.periodicBox());

        // TODO: if we store this as a time series instead of just
        //       accumulating, we could then do PCA and pick out
        //       structural transitions
        // Here incorporate frame weight as numerator of 1/r^6
        get<2>(*p) += wopts->weights() / (d2 * d2 * d2);
      }

      frame++;
    }
  }
  else
  { // If there are no frame weights, don't worry about reweighting in calx.
    while (traj->readFrame())
    {
      traj->updateGroupCoords(system);

      for (vector<AtomPair>::iterator p = hydrogen_pairs.begin();
           p != hydrogen_pairs.end();
           ++p)
      {
        pAtom a1 = get<0>(*p);
        pAtom a2 = get<1>(*p);

        double d2 = a1->coords().distance2(a2->coords(), system.periodicBox());

        // TODO: if we store this as a time series instead of just
        //       accumulating, we could then do PCA and pick out
        //       structural transitions
        get<2>(*p) += 1.0 / (d2 * d2 * d2);
      }

      frame++;
    }
  }
  // sort into ascending order by noe amplitude
  sort(hydrogen_pairs.begin(), hydrogen_pairs.end(), &cmpAtomPairs);

  const double threshold = 1. / (8 * 8 * 8 * 8 * 8 * 8);

  for (vector<AtomPair>::reverse_iterator p = hydrogen_pairs.rbegin();
       p != hydrogen_pairs.rend();
       ++p)
  {
    pAtom a1 = get<0>(*p);
    pAtom a2 = get<1>(*p);

    double val = get<2>(*p) / frame;
    double val6 = pow(1 / val, 1 / 6.);

    if (val > threshold)
    {
      cout << val << "\t" << val6 << "\t"
           << a1->resname() << a1->resid() << ":" << a1->name() << "\t"
           << " ( " << a1->index() << " ) "
           << "\t"
           << a2->resname() << a2->resid() << ":" << a2->name()
           << " ( " << a2->index() << " ) "
           << "\t"
           << endl;
    }
  }
}
