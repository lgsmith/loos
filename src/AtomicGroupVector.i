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

%include "numpy.i"

%rename(cpp_hardContactHistogram) loos::hardContactHistogram;
%rename(cpp_logisticContactHistogram) loos::logisticContactHistogram;
%rename(cpp_hardContactMatrix) loos::hardContactMatrix;
%rename(cpp_logisticContactMatrix) loos::logisticContactMatrix;

%header %{
#include <AtomicGroupVector.hpp>
#if defined(SWIGPYTHON)
// namespace loos {
//   struct IntArr1DCan {
//     std::vector<int> data;
//     void to_numpy(int **out_buff, int *len) { *out_buff = data.data(); }
//   };
//   struct GrealArr1DCan {
//     std::vector<greal> data;
//     void to_numpy(double**out_buff, int *len) { *out_buff = data.data(); }
//   };
//   IntArr1DCan hard_contact_histogram(const std::vector<AtomicGroup> &contacted_vec,
//                               const std::vector<AtomicGroup> &contactor_vec,
//                               greal radius) {
//     IntArr1DCan can;
//     can.data = loos::hardContactHistogram(contacted_vec, contactor_vec, radius);
//     return (can);
//   }
// }; // namespace loos

#endif
%}

%apply (double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_buff, int* len)};
%apply (int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(int** out_buff, int* len)};
%apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** outseq, int* m, int* n)};
%apply (int** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(int** outseq, int* m, int* n)};

%include "AtomicGroupVector.hpp"

namespace loos {
#if defined(SWIGPYTHON)
// struct IntArr1DCan {
//   std::vector<int> data;
//   void to_numpy(int **out_buff, int * len);
// };
// IntArr1DCan hard_contact_histogram(const std::vector<AtomicGroup> &contacted_vec,
//                               const std::vector<AtomicGroup> &contactor_vec,
//                               greal radius);

%pythoncode %{
  def hardContactHistogram(contacted_vec, contactor_vec, radius):
    return hard_contact_histogram(contacted_vec, contactor_vec, radius).to_numpy()
%}
#endif
}; // namespace loos