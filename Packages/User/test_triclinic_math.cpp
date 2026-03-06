/*
  test_triclinic_math.cpp

  Unit tests for TriclinicBox that require no external data files.
  Tests pass silently; any failure prints a message and exits non-zero.

  This file is part of LOOS.
*/

#include <loos.hpp>
#include <cmath>
#include <iostream>
#include <string>

using namespace loos;

// ---------------------------------------------------------------------------
// Tiny test harness
// ---------------------------------------------------------------------------
static int failures = 0;

static void check(bool cond, const std::string& label) {
    if (!cond) {
        std::cerr << "FAIL: " << label << "\n";
        ++failures;
    }
}

static void checkClose(double a, double b, const std::string& label, double tol = 1e-8) {
    if (std::fabs(a - b) > tol) {
        std::cerr << "FAIL: " << label
                  << " (got " << a << ", expected " << b << ")\n";
        ++failures;
    }
}

static void checkCoordClose(const GCoord& got, const GCoord& expected,
                            const std::string& label, double tol = 1e-8) {
    checkClose(got.x(), expected.x(), label + " x", tol);
    checkClose(got.y(), expected.y(), label + " y", tol);
    checkClose(got.z(), expected.z(), label + " z", tol);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Returns true if r is in the primary cell [-L/2, L/2) along each axis
// defined by the lattice H.
static bool inPrimaryCell(const GCoord& r, const TriclinicBox& cell, double tol = 1e-8) {
    // Convert to fractional coordinates
    GCoord tmp = r;
    GCoord origin(0, 0, 0);
    // Use distance2 to check: a point in the primary cell should not move
    // when reimaged.
    GCoord reimaged = r;
    cell.reimage(reimaged);
    return (reimaged - r).length() < tol;
}

// ---------------------------------------------------------------------------
// Test 1: Orthorhombic equivalence
//   TriclinicBox(GCoord(lx, ly, lz)).reimage() must give the same result
//   as the original per-axis algorithm.
// ---------------------------------------------------------------------------
static void testOrthorhombicEquivalence() {
    const double lx = 10.0, ly = 20.0, lz = 30.0;
    GCoord box(lx, ly, lz);
    TriclinicBox cell(box);

    check(cell.isOrthorhombic(), "isOrthorhombic detects diagonal box");

    // Known case: point (13, -24, 37) with box (10, 20, 30) -> (3, -4, 7).
    // Chosen to avoid exact ±L/2 boundaries where floor() rounding can differ
    // between the two algorithms.
    GCoord pt(13.0, -24.0, 37.0);
    GCoord expected(3.0, -4.0, 7.0);

    // Old algorithm result
    GCoord old = pt;
    old.reimage(box);

    // New algorithm result
    GCoord neo = pt;
    cell.reimage(neo);

    checkCoordClose(neo, old,      "orthorhombic: triclinic matches old algorithm");
    checkCoordClose(neo, expected, "orthorhombic: known answer (13,-24,37) -> (3,-4,7)");

    // Several more randomish points
    // Avoid exact ±L/2 boundaries (e.g. 5, 10, 15 for this box) where
    // floor() rounding can differ between the matrix-inverse path and the
    // simple-division path.
    struct { double x, y, z; } pts[] = {
        {  0.0,  0.0,  0.0},   // origin — always at center
        {  4.9,  9.9, 14.9},   // just inside cell
        { -6.0, 11.0, -16.0},
        { 25.0, -41.0, 91.0},
    };
    for (auto& p : pts) {
        GCoord q(p.x, p.y, p.z);
        GCoord qold = q;  qold.reimage(box);
        GCoord qnew = q;  qnew.reimage(cell);
        checkCoordClose(qnew, qold, "orthorhombic equivalence for varied points");
    }
}

// ---------------------------------------------------------------------------
// Test 2: Rhombic dodecahedron (xy-square orientation, GROMACS convention)
//
//   Lattice vectors for box size d:
//     a = (d, 0, 0)
//     b = (0, d, 0)
//     c = (d/2, d/2, d/sqrt(2))
//
//   Verify that:
//   (a) points that differ by a lattice vector reimage to the same place
//   (b) a known point outside the cell maps to its known image
//   (c) the reimaged point stays put when reimaged a second time (idempotence)
// ---------------------------------------------------------------------------
static void testRhombicDodecahedron() {
    const double d = 10.0;
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    GCoord a(d, 0, 0);
    GCoord b(0, d, 0);
    GCoord c(d/2, d/2, d * inv_sqrt2);
    TriclinicBox cell(a, b, c);

    check(!cell.isOrthorhombic(), "rhombic dodecahedron is not orthorhombic");

    // (a) A point and the same point shifted by lattice vector a should
    //     reimage to the same location.
    GCoord p0(1.0, 2.0, 1.0);
    GCoord p1 = p0 + a;   // shifted by one a
    GCoord p2 = p0 - b;   // shifted by minus b
    GCoord p3 = p0 + c;   // shifted by c

    GCoord r0 = p0;  cell.reimage(r0);
    GCoord r1 = p1;  cell.reimage(r1);
    GCoord r2 = p2;  cell.reimage(r2);
    GCoord r3 = p3;  cell.reimage(r3);

    checkCoordClose(r1, r0, "rhombic dodec: p + a -> same image", 1e-7);
    checkCoordClose(r2, r0, "rhombic dodec: p - b -> same image", 1e-7);
    checkCoordClose(r3, r0, "rhombic dodec: p + c -> same image", 1e-7);

    // (b) Point (12, 3, 0) should reimage to (2, 3, 0) because it is
    //     exactly one a-vector outside the cell.
    GCoord pt(12.0, 3.0, 0.0);
    GCoord expected(2.0, 3.0, 0.0);
    GCoord result = pt;
    cell.reimage(result);
    checkCoordClose(result, expected, "rhombic dodec: (12,3,0) -> (2,3,0)", 1e-7);

    // (c) Idempotence: reimaging an already-imaged point is a no-op.
    GCoord r0_again = r0;
    cell.reimage(r0_again);
    checkCoordClose(r0_again, r0, "rhombic dodec: reimage is idempotent", 1e-7);
}

// ---------------------------------------------------------------------------
// Test 3: AtomicGroup::reimage() uses TriclinicBox
//   Build a one-atom group outside the cell, set a TriclinicBox, reimage it,
//   and confirm the atom moved to the correct image.
// ---------------------------------------------------------------------------
static void testAtomicGroupReimage() {
    const double lx = 12.0, ly = 12.0, lz = 12.0;
    GCoord a(lx, 0, 0);
    GCoord b(0, ly, 0);
    GCoord c(lx/2, ly/2, lz/std::sqrt(2.0));
    TriclinicBox cell(a, b, c);

    pAtom pa(new Atom);
    pa->index(0);
    // Place the atom at a+origin, i.e. one box length outside in x
    pa->coords(GCoord(lx + 1.0, 2.0, 1.0));

    AtomicGroup ag;
    ag.append(pa);
    ag.periodicBox(cell);

    ag.reimage();

    GCoord expected(1.0, 2.0, 1.0);
    checkCoordClose(ag[0]->coords(), expected,
                    "AtomicGroup::reimage() with TriclinicBox", 1e-7);
}

// ---------------------------------------------------------------------------
// Test 4: distance2 with TriclinicBox
//   In a box of size 10, a particle at (9,0,0) is only 1 unit from the
//   particle at (0,0,0) when periodicity is considered.
// ---------------------------------------------------------------------------
static void testDistance2() {
    const double L = 10.0;
    TriclinicBox cell(GCoord(L, L, L));

    GCoord p1(0.0, 0.0, 0.0);
    GCoord p2(9.0, 0.0, 0.0);  // 1 away via PBC, 9 away naively

    double d2_pbc  = p1.distance2(p2, cell);
    double d2_naive = p1.distance2(p2);

    checkClose(d2_pbc,   1.0, "distance2 with TriclinicBox: PBC image distance");
    checkClose(d2_naive, 81.0, "distance2 with TriclinicBox: naive distance");

    // Same check for non-orthorhombic: rhombic dodecahedron, d=10
    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    TriclinicBox rdodec(GCoord(L, 0, 0),
                        GCoord(0, L, 0),
                        GCoord(L/2, L/2, L*inv_sqrt2));

    // Point (9, 0, 0) is 1 away from origin via the a-image
    double d2_rd = p1.distance2(p2, rdodec);
    checkClose(d2_rd, 1.0, "distance2 with rhombic dodecahedron", 1e-7);
}

// ---------------------------------------------------------------------------
// Test 5: TriclinicBox construction from raw float[9] array (XTC convention)
//   GROMACS stores rows as lattice vectors, scale nm->Angstrom (x10).
// ---------------------------------------------------------------------------
static void testRawArrayConstruction() {
    // Orthorhombic 5x6x7 nm box stored as GROMACS would (nm units)
    float raw[9] = {5.0f, 0.0f, 0.0f,
                    0.0f, 6.0f, 0.0f,
                    0.0f, 0.0f, 7.0f};
    TriclinicBox cell(raw, 10.0);  // scale nm -> Angstrom

    check(cell.isOrthorhombic(), "raw array orthorhombic box detected");
    checkClose(cell.a().length(), 50.0, "raw array: a length = 50 Angstrom");
    checkClose(cell.b().length(), 60.0, "raw array: b length = 60 Angstrom");
    checkClose(cell.c().length(), 70.0, "raw array: c length = 70 Angstrom");

    // Roundtrip: a point 1 nm outside in x should image correctly
    GCoord pt(55.0, 0.0, 0.0);
    cell.reimage(pt);
    checkClose(pt.x(), 5.0, "raw array: reimage across x boundary");
}

// ---------------------------------------------------------------------------
int main() {
    testOrthorhombicEquivalence();
    testRhombicDodecahedron();
    testAtomicGroupReimage();
    testDistance2();
    testRawArrayConstruction();

    if (failures == 0) {
        std::cout << "All TriclinicBox math tests passed.\n";
        return 0;
    } else {
        std::cerr << failures << " test(s) failed.\n";
        return 1;
    }
}
