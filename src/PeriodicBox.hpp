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




#if !defined(LOOS_PERIODICBOX_HPP)
#define LOOS_PERIODICBOX_HPP



#include <boost/shared_ptr.hpp>

#include <loos_defs.hpp>
#include <Coord.hpp>



namespace loos {

  // ---------------------------------------------------------------------------
  //! General triclinic unit cell described by three lattice vectors a, b, c.
  /**
   * Stores the cell as a 3x3 matrix H whose columns are the lattice vectors,
   * along with its precomputed inverse.  Provides reimaging via the fractional
   * coordinate approach, which works for any valid (non-singular) parallelotope
   * unit cell including orthorhombic, triclinic, rhombic dodecahedron, and
   * truncated octahedron cells.
   *
   * The primary cell is defined as fractional coordinates in [-0.5, 0.5),
   * matching the GROMACS convention.
   */
  class TriclinicBox {
  public:
    //! Construct from three explicit lattice vectors
    TriclinicBox(const GCoord& a, const GCoord& b, const GCoord& c);

    //! Orthorhombic convenience constructor from box edge lengths
    explicit TriclinicBox(const GCoord& lengths);

    //! Construct from a 9-element row-major float array (XTC/TRR format).
    //! @param box   9-element array: rows are a, b, c vectors
    //! @param scale unit conversion factor (e.g. 10.0 for nm->Angstrom)
    TriclinicBox(const float box[9], double scale = 1.0);

    //! The a lattice vector
    GCoord a() const;
    //! The b lattice vector
    GCoord b() const;
    //! The c lattice vector
    GCoord c() const;

    //! Lengths of a, b, c (not necessarily the bounding box for triclinic cells)
    GCoord lengths() const;

    //! True when all off-diagonal elements of H are smaller than tol * min(diagonal)
    bool isOrthorhombic(double tol = 1e-6) const;

    //! Reimage r into the primary cell [-0.5, 0.5) in fractional coordinates.
    //! r is modified in place.
    void reimage(GCoord& r) const;

  private:
    void setFromVectors(const GCoord& a, const GCoord& b, const GCoord& c);
    void computeInverse();

    double H[3][3];     //!< columns are a, b, c lattice vectors
    double Hinv[3][3];  //!< precomputed inverse of H
  };


  // ---------------------------------------------------------------------------
  //! Class for managing periodic box information.
  /** This is the fundamental object that gets shared amongst related
   *  groups.  It contains the GCoord representing the box size and a
   *  flag that indicates whether or not the box has actually been set.
   *  The client will not interact with this class/object directly, but
   *  will use the SharedPeriodicBox instead.
   */

  class PeriodicBox {
  public:
    PeriodicBox() : thebox(99999,99999,99999), box_set(false), has_triclinic(false) { }
    explicit PeriodicBox(const GCoord& c) : thebox(c), box_set(true), has_triclinic(false) { }

    GCoord box(void) const { return(thebox); }
    void box(const GCoord& c) {
      thebox = c;
      has_triclinic = false;

      // Because of the way boxes are handled elsewhere, setting an
      // unset box in AtomicGroup can leave PeriodicBox thinking it
      // has been set, when it's really just the default value.  So,
      // check on set for the magic value and if we find it, then
      // assume we're being set by an unset box.
      box_set = (c.x() != 99999 || c.y() != 99999 || c.z() != 99999);
    }

    bool isPeriodic(void) const { return(box_set); }
    void setPeriodic(const bool b) { box_set = b; }

    //! Store a general triclinic cell.  Also updates the diagonal GCoord
    //! (returned by box()) to the lattice vector lengths for backward compat.
    void setTriclinicBox(const TriclinicBox& tb) {
      triclinic_box = tb;
      has_triclinic = true;
      thebox = tb.lengths();
      box_set = true;
    }

    bool hasTriclinicBox(void) const { return(has_triclinic); }

    const TriclinicBox& triclinicBox(void) const { return(triclinic_box); }

  private:
    GCoord thebox;
    bool box_set;
    bool has_triclinic;
    TriclinicBox triclinic_box{GCoord(1,1,1)};  // default placeholder
  };


  //! This class manages a shared Periodicbox
  /** This is what most clients should use.  It maintains a shared
   * resource (the shared PeriodicBox) and forwards member function
   * calls to it.  The copy() member function creates a new
   * (i.e. dissociated PeriodicBox) and returns the associated
   * SharedPeriodicBox
   */

  class SharedPeriodicBox {
  public:
    SharedPeriodicBox() : pbox(new PeriodicBox) { }
    GCoord box(void) const { return(pbox->box()); }
    void box(const GCoord& c) { pbox->box(c); }
    bool isPeriodic(void) const { return(pbox->isPeriodic()); }

    void setTriclinicBox(const TriclinicBox& tb) { pbox->setTriclinicBox(tb); }
    bool hasTriclinicBox(void) const { return(pbox->hasTriclinicBox()); }
    const TriclinicBox& triclinicBox(void) const { return(pbox->triclinicBox()); }


    SharedPeriodicBox copy(void) const {
      SharedPeriodicBox thecopy;

      if (isPeriodic()) {
        if (hasTriclinicBox())
          thecopy.setTriclinicBox(triclinicBox());
        else
          thecopy.box(box());
      }

      return(thecopy);
    }

  private:
    boost::shared_ptr<PeriodicBox> pbox;
  };

}

#endif
