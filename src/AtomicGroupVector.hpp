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
#if !defined(AG_VEC)
#define AG_VEC
#include <AtomicGroup.hpp>
#include <vector>

namespace loos {

// Applies 'hardContact' to each of the elements in contacted vs
// each of the elements in contactor.
std::vector<int>
hardContactHistogram(const std::vector<AtomicGroup> &contacted_vec,
                     const std::vector<AtomicGroup> &contactor_vec,
                     greal radius, const GCoord &box) {
  std::vector<int> contact_hist(contacted_vec.size());
  int contacts;
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
  int contacts;
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
std::vector<int>
hardContactHistogram(const std::vector<AtomicGroup> &contacted_vec,
                     const std::vector<AtomicGroup> &contactor_vec,
                     greal radius) {
  std::vector<int> contact_hist(contacted_vec.size());
  int contacts;
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
  greal contacts;
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
std::vector<std::vector<int>>
hardContactMatrix(const std::vector<AtomicGroup> &contacted_vec,
                  const std::vector<AtomicGroup> &contactor_vec, greal radius,
                  int sigma) {
  std::vector<std::vector<int>> contact_matrix(
      contacted_vec.size(), std::vector<int>(contactor_vec.size()));
  for (auto contacted : contacted_vec) {
    std::vector<int> contacts(contacted_vec.size());
    for (auto contactor : contactor_vec) {
      contacts.emplace_back(contacted.hardContact(contactor, radius));
    }
    contact_matrix.emplace_back(contacts);
  }
  return (contact_matrix);
}

// Record hard contacts (with box info) as a matrix:
// (all possible contact groups)X(all contactors)
std::vector<std::vector<int>>
hardContactMatrix(const std::vector<AtomicGroup> &contacted_vec,
                  const std::vector<AtomicGroup> &contactor_vec, greal radius,
                  const GCoord &box) {
  std::vector<std::vector<int>> contact_matrix(
      contacted_vec.size(), std::vector<int>(contactor_vec.size()));
  for (auto contacted : contacted_vec) {
    std::vector<int> contacts(contacted_vec.size());
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
                      greal radius, int sigma, const GCoord &box) {
  std::vector<std::vector<greal>> contact_matrix(
      contacted_vec.size(), std::vector<greal>(contactor_vec.size()));
  for (auto contacted : contacted_vec) {
    std::vector<greal> contacts(contacted_vec.size());
    for (auto contactor : contactor_vec) {
      contacts.emplace_back(
          contacted.logisticContact(contactor, radius, sigma, box));
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
  std::vector<std::vector<greal>> contact_matrix(
      contacted_vec.size(), std::vector<greal>(contactor_vec.size()));
  for (auto contacted : contacted_vec) {
    std::vector<greal> contacts(contacted_vec.size());
    for (auto contactor : contactor_vec) {
      contacts.emplace_back(
          contacted.logisticContact(contactor, radius, sigma));
    }
    contact_matrix.emplace_back(contacts);
  }
  return (contact_matrix);
}


  // All this stuff is for type-uplifting in numpy
  struct IntArr1DCan {
    std::vector<int> data;
    void to_numpy(int **out_buff, int *len) { *out_buff = data.data(); }
  };
  struct GrealArr1DCan {
    std::vector<greal> data;
    void to_numpy(double**out_buff, int *len) { *out_buff = data.data(); }
  };
  IntArr1DCan hard_contact_histogram(const std::vector<AtomicGroup> &contacted_vec,
                              const std::vector<AtomicGroup> &contactor_vec,
                              greal radius) {
    IntArr1DCan can;
    can.data = hardContactHistogram(contacted_vec, contactor_vec, radius);
    return (can);
  }

}; // namespace loos
#endif