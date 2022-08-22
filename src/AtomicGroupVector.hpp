/*
  AtomicGroupVector.hpp

  functions that operate on vectors of AGs.
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

#include <AtomicGroup.hpp>
#include <vector>

namespace loos {

// Applies 'hardContact' to each of the elements in contacted vs
// each of the elements in contactor.
std::vector<uint>
hardContactHistogram(const std::vector<AtomicGroup> &contacted_vec,
                  const std::vector<AtomicGroup> &contactor_vec, greal radius,
                  const GCoord &box) {
  std::vector<uint> contact_hist(contacted_vec.size());
  uint contacts;
  for (auto contacted : contacted_vec) {
    contacts = 0;
    for (auto contactor : contactor_vec) {
      contacts += contacted.hardContact(contactor, radius, box);
    }
    contact_hist.emplace_back(contacts);
  }
  return (contact_hist);
}

// Applies 'logisticContact' to each of the elements in contacted vs
// each of the elements in contactor.
std::vector<greal>
logisticContactHistogram(const std::vector<AtomicGroup> &contacted_vec,
                        const std::vector<AtomicGroup> &contactor_vec,
                        greal radius, int sigma, const GCoord &box) {
  std::vector<greal> contact_hist(contacted_vec.size());
  uint contacts;
  for (auto contacted : contacted_vec) {
    contacts = 0;
    for (auto contactor : contactor_vec) {
      contacts += contacted.logisticContact(contactor, radius, sigma, box);
    }
    contact_hist.emplace_back(contacts);
  }
  return (contact_hist);
}

// Applies hardContact without box info
std::vector<uint>
hardContactHistogram(const std::vector<AtomicGroup> &contacted_vec,
                  const std::vector<AtomicGroup> &contactor_vec, greal radius) {
  std::vector<uint> contact_hist(contacted_vec.size());
  uint contacts;
  for (auto contacted : contacted_vec) {
    contacts = 0;
    for (auto contactor : contactor_vec) {
      contacts += contacted.hardContact(contactor, radius);
    }
    contact_hist.push_back(contacts);
  }
  return (contact_hist);
}

// Applies logisticContact without box info
std::vector<greal>
logisticContactHistogram(const std::vector<AtomicGroup> &contacted_vec,
                        const std::vector<AtomicGroup> &contactor_vec,
                        greal radius, int sigma) {
  std::vector<greal> contact_hist(contacted_vec.size());
  uint contacts;
  for (auto contacted : contacted_vec) {
    contacts = 0;
    for (auto contactor : contactor_vec) {
      contacts += contacted.logisticContact(contactor, radius, sigma);
    }
    contact_hist.emplace_back(contacts);
  }
  return (contact_hist);
}

// Record logistic contacts as a matrix: 
// (all possible contact groups)X(all contactors)
std::vector<std::vector<greal>>
hardContactMatrix(const std::vector<AtomicGroup> &contacted_vec,
                     const std::vector<AtomicGroup> &contactor_vec,
                     greal radius, int sigma) {
  std::vector<std::vector < greal> > contact_matrix(contacted_vec.size(), std::vector<greal>(contactor_vec.size())); 
  for (auto contacted : contacted_vec) {
    std::vector<greal> contacts(size(contacted_vec));
    for (auto contactor : contactor_vec) {
      contacts.emplace_back(contacted.hardContact(contactor, radius));
    }
    contact_matrix.emplace_back(contacts);
  }
  return (contact_matrix);
}

// Record hard contacts (with box info) as a matrix: 
// (all possible contact groups)X(all contactors)
std::vector<std::vector<greal>>
hardContactMatrix(const std::vector<AtomicGroup> &contacted_vec,
                     const std::vector<AtomicGroup> &contactor_vec,
                     greal radius, const GCoord & box) {
  std::vector<std::vector < greal> > contact_matrix(contacted_vec.size(), std::vector<greal>(contactor_vec.size())); 
  for (auto contacted : contacted_vec) {
    std::vector<greal> contacts(size(contacted_vec));
    for (auto contactor : contactor_vec) {
      contacts.emplace_back(contacted.hardContact(contactor, radius, box));
    }
    contact_matrix.emplace_back(contacts);
  }
  return (contact_matrix);
}

// write hard contacts to a matrix: 
// (all possible contact groups)X(all contactors)
std::vector<std::vector<greal>>
logisticContactMatrix(const std::vector<AtomicGroup> &contacted_vec,
                     const std::vector<AtomicGroup> &contactor_vec,
                     greal radius, int sigma, const GCoord & box) {
  std::vector<std::vector < greal> > contact_matrix(contacted_vec.size(), std::vector<greal>(contactor_vec.size())); 
  for (auto contacted : contacted_vec) {
    std::vector<greal> contacts(size(contacted_vec));
    for (auto contactor : contactor_vec) {
      contacts.emplace_back(contacted.logisticContact(contactor, radius, sigma, box));
    }
    contact_matrix.emplace_back(contacts);
  }
  return (contact_matrix);
}

// write logistic contacts to a matrix: 
// (all possible contact groups)X(all contactors)
std::vector<std::vector<greal>>
logisticContactMatrix(const std::vector<AtomicGroup> &contacted_vec,
                     const std::vector<AtomicGroup> &contactor_vec,
                     greal radius, int sigma) {
  std::vector<std::vector < greal> > contact_matrix(contacted_vec.size(), std::vector<greal>(contactor_vec.size())); 
  for (auto contacted : contacted_vec) {
    std::vector<greal> contacts(size(contacted_vec));
    for (auto contactor : contactor_vec) {
      contacts.emplace_back(contacted.logisticContact(contactor, radius, sigma));
    }
    contact_matrix.emplace_back(contacts);
  }
  return (contact_matrix);
}

} // namespace loos