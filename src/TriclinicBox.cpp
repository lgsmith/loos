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

#include <cmath>
#include <stdexcept>
#include <sstream>
// Coord.hpp must come before PeriodicBox.hpp so GCoord is a complete type
// when the inline methods in PeriodicBox.hpp are instantiated.
#include "Coord.hpp"
#include "PeriodicBox.hpp"

namespace loos {

// ---------------------------------------------------------------------------
// Construction helpers
// ---------------------------------------------------------------------------

// Store lattice vectors as the columns of H:
//   H = [ a | b | c ]
// so H[row][col]: H[0][0]=ax, H[1][0]=ay, H[2][0]=az, etc.
void TriclinicBox::setFromVectors(const GCoord& a, const GCoord& b, const GCoord& c) {
    H[0][0] = a.x();  H[0][1] = b.x();  H[0][2] = c.x();
    H[1][0] = a.y();  H[1][1] = b.y();  H[1][2] = c.y();
    H[2][0] = a.z();  H[2][1] = b.z();  H[2][2] = c.z();
    computeInverse();
}

TriclinicBox::TriclinicBox(const GCoord& a, const GCoord& b, const GCoord& c) {
    setFromVectors(a, b, c);
}

// Orthorhombic convenience constructor: a=(lx,0,0), b=(0,ly,0), c=(0,0,lz)
TriclinicBox::TriclinicBox(const GCoord& lengths) {
    setFromVectors(GCoord(lengths.x(), 0.0, 0.0),
                   GCoord(0.0, lengths.y(), 0.0),
                   GCoord(0.0, 0.0, lengths.z()));
}

// Construct from a raw 9-element row-major float array as stored by XTC/TRR.
// GROMACS convention: the matrix rows are the lattice vectors, so
//   a = (box[0], box[1], box[2])
//   b = (box[3], box[4], box[5])
//   c = (box[6], box[7], box[8])
// scale converts units (e.g. nm->Angstrom = 10.0).
TriclinicBox::TriclinicBox(const float box[9], double scale) {
    GCoord a(box[0]*scale, box[1]*scale, box[2]*scale);
    GCoord b(box[3]*scale, box[4]*scale, box[5]*scale);
    GCoord c(box[6]*scale, box[7]*scale, box[8]*scale);
    setFromVectors(a, b, c);
}

// ---------------------------------------------------------------------------
// Matrix inverse (analytic 3x3 via cofactors)
// ---------------------------------------------------------------------------
void TriclinicBox::computeInverse() {
    double a00=H[0][0], a01=H[0][1], a02=H[0][2];
    double a10=H[1][0], a11=H[1][1], a12=H[1][2];
    double a20=H[2][0], a21=H[2][1], a22=H[2][2];

    double det = a00*(a11*a22 - a12*a21)
               - a01*(a10*a22 - a12*a20)
               + a02*(a10*a21 - a11*a20);

    if (std::fabs(det) < 1e-10) {
        std::ostringstream oss;
        oss << "TriclinicBox: singular cell matrix (det=" << det << ")";
        throw std::runtime_error(oss.str());
    }

    double inv_det = 1.0 / det;

    Hinv[0][0] =  (a11*a22 - a12*a21) * inv_det;
    Hinv[0][1] = -(a01*a22 - a02*a21) * inv_det;
    Hinv[0][2] =  (a01*a12 - a02*a11) * inv_det;

    Hinv[1][0] = -(a10*a22 - a12*a20) * inv_det;
    Hinv[1][1] =  (a00*a22 - a02*a20) * inv_det;
    Hinv[1][2] = -(a00*a12 - a02*a10) * inv_det;

    Hinv[2][0] =  (a10*a21 - a11*a20) * inv_det;
    Hinv[2][1] = -(a00*a21 - a01*a20) * inv_det;
    Hinv[2][2] =  (a00*a11 - a01*a10) * inv_det;
}

// ---------------------------------------------------------------------------
// Accessors
// ---------------------------------------------------------------------------
GCoord TriclinicBox::a() const { return GCoord(H[0][0], H[1][0], H[2][0]); }
GCoord TriclinicBox::b() const { return GCoord(H[0][1], H[1][1], H[2][1]); }
GCoord TriclinicBox::c() const { return GCoord(H[0][2], H[1][2], H[2][2]); }

GCoord TriclinicBox::lengths() const {
    return GCoord(a().length(), b().length(), c().length());
}

// A box is orthorhombic when all off-diagonal elements are near zero.
// Tolerance is relative to the smallest diagonal.
bool TriclinicBox::isOrthorhombic(double tol) const {
    double scale = std::min({std::fabs(H[0][0]),
                             std::fabs(H[1][1]),
                             std::fabs(H[2][2])});
    double thresh = tol * scale;
    return (std::fabs(H[0][1]) < thresh &&
            std::fabs(H[0][2]) < thresh &&
            std::fabs(H[1][0]) < thresh &&
            std::fabs(H[1][2]) < thresh &&
            std::fabs(H[2][0]) < thresh &&
            std::fabs(H[2][1]) < thresh);
}

// ---------------------------------------------------------------------------
// Core reimaging: fractional coordinate wrap to [-0.5, 0.5)
// ---------------------------------------------------------------------------
void TriclinicBox::reimage(GCoord& r) const {
    // s = H^{-1} * r  (fractional coordinates)
    double sx = Hinv[0][0]*r.x() + Hinv[0][1]*r.y() + Hinv[0][2]*r.z();
    double sy = Hinv[1][0]*r.x() + Hinv[1][1]*r.y() + Hinv[1][2]*r.z();
    double sz = Hinv[2][0]*r.x() + Hinv[2][1]*r.y() + Hinv[2][2]*r.z();

    // Wrap each fractional coordinate to [-0.5, 0.5)
    sx -= std::floor(sx + 0.5);
    sy -= std::floor(sy + 0.5);
    sz -= std::floor(sz + 0.5);

    // r = H * s  (back to Cartesian)
    r = GCoord(H[0][0]*sx + H[0][1]*sy + H[0][2]*sz,
               H[1][0]*sx + H[1][1]*sy + H[1][2]*sz,
               H[2][0]*sx + H[2][1]*sy + H[2][2]*sz);
}

} // namespace loos

// ---------------------------------------------------------------------------
// Coord<T>::reimage(TriclinicBox) — out-of-line definition.
// Declared in Coord.hpp with only a forward declaration of TriclinicBox;
// the full definitions live here where TriclinicBox is complete.
// greal == double, so only the double specialization is needed for GCoord;
// the float specialization covers any float-typed Coord users.
// ---------------------------------------------------------------------------
template<>
void loos::Coord<double>::reimage(const loos::TriclinicBox& cell) {
    loos::GCoord r(v[0], v[1], v[2]);
    cell.reimage(r);
    v[0] = r.x();  v[1] = r.y();  v[2] = r.z();
}

template<>
void loos::Coord<float>::reimage(const loos::TriclinicBox& cell) {
    loos::GCoord r(v[0], v[1], v[2]);
    cell.reimage(r);
    v[0] = static_cast<float>(r.x());
    v[1] = static_cast<float>(r.y());
    v[2] = static_cast<float>(r.z());
}
