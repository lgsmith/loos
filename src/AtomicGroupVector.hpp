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

#include <vector>
#include <AtomicGroup.hpp>

namespace loos {
  
  // Applies 'hardContact' to each of the elements in contacted vs
  // each of the elements in contactor.
  std::vector<uint> contact_histogram(const std::vector<AtomicGroup> & contacted_vec, const std::vector<AtomicGroup>& contactor_vec, greal radius, const GCoord& box) {
    std::vector<uint> contact_hist(contacted_vec.size());
    uint contacts;
    for (auto contacted : contacted_vec) {
      contacts = 0;
      for (auto contactor : contactor_vec){
        contacts += contacted.hardContact(contactor, radius, box);
      }
      contact_hist.push_back(contacts);
    }
    return(contact_hist);
  }

  // Applies 'logisticContact' to each of the elements in contacted vs
  // each of the elements in contactor.
  std::vector<uint> logit_contact_histogram(const std::vector<AtomicGroup> & contacted_vec, const std::vector<AtomicGroup>& contactor_vec, greal radius, const GCoord& box) {
    std::vector<uint> contact_hist(contacted_vec.size());
    uint contacts;
    for (auto contacted : contacted_vec) {
      contacts = 0;
      for (auto contactor : contactor_vec){
        contacts += contacted.logisticContact(contactor, radius, box);
      }
      contact_hist.push_back(contacts);
    }
    return(contact_hist);
  }
}