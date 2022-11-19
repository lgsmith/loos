/*
  Core code for hbonds tools
*/

/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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

#include <boost/format.hpp>

#include "hcore.hpp"

using namespace loos;
using namespace HBonds;

double SimpleAtom::inner = 0.0;
double SimpleAtom::outer = 3.5;
double SimpleAtom::deviation = 20.0;

// Reports distance^2 between hydrogen and heavy atom
// D-H ... X
//   |-----|
double SimpleAtom::distance2(const SimpleAtom &s) const {
  double d;

  if (usePeriodicity)
    d = atom->coords().distance2(s.atom->coords(), sbox.box());
  else
    d = atom->coords().distance2(s.atom->coords());

  return (d);
}

// Returns angle between atoms in degrees...
//
//  D-H ... X
//   \---/
double SimpleAtom::angle(const SimpleAtom &s) const {
  loos::GCoord left, middle, right;

  if (isHydrogen) {
    if (s.isHydrogen)
      throw(std::runtime_error("Cannot take the angle between two hydrogens"));

    left = attached_to->coords();
    middle = atom->coords();
    right = s.atom->coords();

  } else {
    if (!s.isHydrogen)
      throw(std::runtime_error(
          "Cannot take the angle between two non-hydrogens"));

    left = atom->coords();
    middle = s.atom->coords();
    right = s.attached_to->coords();
  }

  if (usePeriodicity) {
    loos::GCoord box = sbox.box();

    left.reimage(box);
    middle.reimage(box);
    right.reimage(box);
  }

  return (loos::Math::angle(left, middle, right));
}

// Converts a selection into a vector of SimpleAtoms.  All new
// SimpleAtom's get assigned the use_periodicity flag

std::vector<SimpleAtom>
SimpleAtom::processSelection(const std::string &selection,
                             const loos::AtomicGroup &system,
                             const bool use_periodicity) {
  std::vector<SimpleAtom> processed;

  // We don't want to force model to be sorted (esp since it's const)
  // So, create a lightweight copy but use the same pAtoms...
  loos::AtomicGroup searchable(system);
  searchable.sort();

  loos::AtomicGroup model = selectAtoms(system, selection);

  for (loos::AtomicGroup::const_iterator i = model.begin(); i != model.end();
       ++i) {
    SimpleAtom new_atom(*i, system.sharedPeriodicBox(), use_periodicity);
    std::string n = (*i)->name();
    if (n[0] == 'H') {
      new_atom.isHydrogen = true;
      std::vector<int> bond_list = (*i)->getBonds();
      if (bond_list.size() == 0)
        throw(ErrorWithAtom(
            *i, "Detected a hydrogen bond that has no connectivity"));
      if (bond_list.size() != 1)
        throw(ErrorWithAtom(
            *i, "Detected a hydrogen bond that has more than one atom bound"));

      loos::pAtom pa = searchable.findById(bond_list[0]);
      if (pa == 0)
        throw(ErrorWithAtom(*i, "Cannot find atom hydrogen is bound too"));

      new_atom.attached_to = pa;
    }

    processed.push_back(new_atom);
  }

  return (processed);
}

// Check for a hdyrogen bond between two SimpleAtom's

bool SimpleAtom::hydrogenBond(const SAtom &o) const {
  double dist = distance2(o);
  double angl = angle(o);

  if (dist >= inner && dist <= outer &&
      fmod(fabs(angl - 180.0), 360.0) <= deviation)
    return (true);

  return (false);
}

// Returns an AtomicGroup containing the atoms that are hydrogen
// bonded to self.  If findFirstOnly is true, than the first hydrogen
// bond found causes the function to return...  (i.e. it may be a
// small optimization in performance)

loos::AtomicGroup
SimpleAtom::findHydrogenBonds(const std::vector<SimpleAtom> &group,
                              const bool findFirstOnly) {
  loos::AtomicGroup results;

  for (SAGroup::const_iterator i = group.begin(); i != group.end(); ++i) {
    bool b = hydrogenBond(*i);
    if (b) {
      results.append(i->atom);
      if (findFirstOnly)
        break;
    }
  }

  return (results);
}

// Returns a vector of flags indicating which SimpleAtoms form a
// hydrogen bond to self.

std::vector<uint>
SimpleAtom::findHydrogenBondsVector(const std::vector<SimpleAtom> &group) {
  std::vector<uint> results;

  for (SAGroup::const_iterator i = group.begin(); i != group.end(); ++i)
    results.push_back(hydrogenBond(*i));

  return (results);
}

// Returns a matrix of flags indicating which SimpleAtoms form a
// hydrogen bond to self of a trajectory

BondMatrix
SimpleAtom::findHydrogenBondsMatrix(const std::vector<SimpleAtom> &group,
                                    loos::pTraj &traj, loos::AtomicGroup &model,
                                    const uint maxt) const {

  if (maxt > traj->nframes()) {
    std::cerr << boost::format(
                     "Error- row clip (%d) exceeds trajectory size (%d)\n") %
                     maxt % traj->nframes();
    exit(-10);
  }

  BondMatrix bonds(maxt, group.size());

  for (uint t = 0; t < maxt; ++t) {
    traj->readFrame(t);
    traj->updateGroupCoords(model);

    for (uint i = 0; i < group.size(); ++i)
      bonds(t, i) = hydrogenBond(group[i]);
  }

  return (bonds);
}


vBond findPotentialBonds(const AtomicGroup &donors,
                         const AtomicGroup &acceptors,
                         const AtomicGroup &system,
                         const bool use_periodicity=false) {
  vBond bonds;

  for (AtomicGroup::const_iterator j = donors.begin(); j != donors.end(); ++j) {
    GCoord u = (*j)->coords();
    for (AtomicGroup::const_iterator i = acceptors.begin();
         i != acceptors.end(); ++i) {
      // Manually build simple atoms
      SimpleAtom new_donor(*j, system.sharedPeriodicBox(), use_periodicity);
      string name = (*j)->name();
      if (name[0] != 'H') {
        cerr << boost::format("Error- atom %s was given as a donor, but donors "
                              "can only be hydrogens.\n") %
                    name;
        exit(-10);
      }

      std::vector<int> bond_list = (*j)->getBonds();
      if (bond_list.size() != 1) {
        cerr << "Error- The following hydrogen atom has more than one bond to "
                "it...woops...\n";
        cerr << *j;
        exit(-10);
      }

      pAtom pa = system.findById(bond_list[0]);
      if (pa == 0) {
        cerr << boost::format("Error- cannot find atomid %d in system.\n") %
                    bond_list[0];
        exit(-10);
      }
      new_donor.attach(pa);

      SimpleAtom new_acceptor(*i, system.sharedPeriodicBox(), use_periodicity);

      Bond new_bond(new_donor, new_acceptor);
      bonds.push_back(new_bond);
    }
  }

  return (bonds);
}


vBond findPotentialBonds(const AtomicGroup &donors,
                         const AtomicGroup &acceptors,
                         const AtomicGroup &system,
                         const hreal putative_threshold,
                         const bool use_periodicity = false) {
  vBond bonds;

  for (AtomicGroup::const_iterator j = donors.begin(); j != donors.end(); ++j) {
    GCoord u = (*j)->coords();
    for (AtomicGroup::const_iterator i = acceptors.begin();
         i != acceptors.end(); ++i) {
      if (u.distance((*i)->coords()) <= putative_threshold) {

        // Manually build simple atoms
        SimpleAtom new_donor(*j, system.sharedPeriodicBox(), use_periodicity);
        string name = (*j)->name();
        if (name[0] != 'H') {
          cerr << boost::format("Error- atom %s was given as a donor, but "
                                "donors can only be hydrogens.\n") %
                      name;
          exit(-10);
        }

        vector<int> bond_list = (*j)->getBonds();
        if (bond_list.size() != 1) {
          cerr << "Error- The following hydrogen atom has more than one bond "
                  "to it...woops...\n";
          cerr << *j;
          exit(-10);
        }

        pAtom pa = system.findById(bond_list[0]);
        if (pa == 0) {
          cerr << boost::format("Error- cannot find atomid %d in system.\n") %
                      bond_list[0];
          exit(-10);
        }
        new_donor.attach(pa);

        SimpleAtom new_acceptor(*i, system.sharedPeriodicBox(),
                                use_periodicity);

        Bond new_bond(new_donor, new_acceptor);
        bonds.push_back(new_bond);
      }
    }
  }

  return (bonds);
}

std::ostringstream &formatBond(std::ostringstream &oss, const uint i, const Bond &bond) {
  pAtom a = bond.first.rawAtom();
  pAtom b = bond.second.rawAtom();

  oss << boost::format("# %d : %d-%s-%s-%d-%s => %d-%s-%s-%d-%s") % i %
             a->id() % a->name() % a->resname() % a->resid() % a->segid() %
             b->id() % b->name() % b->resname() % b->resid() % b->segid();

  return (oss);
}

bool SimpleAtom::divineHydrogen(const std::string &name) {
  if (name[0] == 'H')
    return (true);

  return (false);
}

vGroup splitSelection(const vGroup &molecules, const std::string &selection) {
  vGroup results;

  for (vGroup::const_iterator i = molecules.begin(); i != molecules.end();
       ++i) {
    AtomicGroup subset;
    try {
      subset = selectAtoms(*i, selection);
    } catch (...) { // Ignore exceptions
      ;
    }
    if (!subset.empty())
      results.push_back(subset);
  }

  if (results.empty()) {
    cerr << "Error- The selection '" << selection
         << "' resulted in nothing being selected.\n";
    exit(-1);
  }

  return (results);
}

vGroup splitSelectionKeepEmpties(const vGroup &molecules,
                                 const std::string &selection) {
  vGroup results;

  for (vGroup::const_iterator i = molecules.begin(); i != molecules.end();
       ++i) {
    AtomicGroup subset;
    try {
      subset = selectAtoms(*i, selection);
    } catch (...) { // Ignore exceptions
      ;
    }
    // Always push back, even empty AGs, as placeholders.
    results.push_back(subset);
  }

  return (results);
}

