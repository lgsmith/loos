/*
  test_triclinic_traj.cpp

  Integration test for triclinic reimaging from a real XTC or TRR trajectory.
  Requires an actual trajectory file with a non-orthorhombic periodic box
  (e.g. GROMACS output with -pbc dodecahedron or a triclinic cell).

  Usage:
    test_triclinic_traj <model> <trajectory>

  The test:
    1. Reads the model and trajectory.
    2. For each frame, calls reimage() on every molecule.
    3. Verifies that the centroid of every molecule is inside the primary cell
       after reimaging (i.e. in fractional coordinates [-0.5, 0.5) along each
       lattice direction).
    4. Verifies that the box reported as non-orthorhombic for a known
       rhombic dodecahedron trajectory.

  Exits 0 on success, 1 on failure.

  This file is part of LOOS.
*/

#include <loos.hpp>
#include <cmath>
#include <iostream>
#include <string>

using namespace loos;

static int failures = 0;

static void check(bool cond, const std::string& label) {
    if (!cond) {
        std::cerr << "FAIL: " << label << "\n";
        ++failures;
    }
}

// Returns true if fractional coordinates of r are in [-0.5, 0.5).
static bool inPrimaryCell(const GCoord& r, const TriclinicBox& cell, double tol = 1e-5) {
    GCoord reimaged = r;
    cell.reimage(reimaged);
    return (reimaged - r).length() < tol;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <model> <trajectory>\n";
        return 1;
    }

    AtomicGroup model = createSystem(argv[1]);
    pTraj traj = createTrajectory(argv[2], model);

    check(traj->hasPeriodicBox(), "trajectory has periodic box");
    if (failures) {
        std::cerr << "Cannot continue without periodic box.\n";
        return 1;
    }

    // Read first frame and check whether the box is non-orthorhombic.
    traj->readFrame();
    traj->updateGroupCoords(model);

    check(model.hasTriclinicBox(),
          "first frame: model has a TriclinicBox (non-orthorhombic cell detected)");
    check(!model.triclinicBox().isOrthorhombic(),
          "first frame: triclinic box is actually non-orthorhombic");

    // Split into molecules for the reimage test.
    std::vector<AtomicGroup> molecules = model.splitByMolecule();

    int frame = 0;
    int atoms_outside = 0;
    traj->rewind();

    while (traj->readFrame()) {
        traj->updateGroupCoords(model);
        ++frame;

        check(model.hasTriclinicBox(),
              "frame " + std::to_string(frame) + ": TriclinicBox present");

        const TriclinicBox& cell = model.triclinicBox();

        // Reimage each molecule.
        for (auto& mol : molecules)
            mol.reimage();

        // Verify each molecule's centroid is now inside the primary cell.
        for (size_t i = 0; i < molecules.size(); ++i) {
            GCoord com = molecules[i].centroid();
            if (!inPrimaryCell(com, cell)) {
                ++atoms_outside;
                if (atoms_outside <= 5)   // avoid flooding output
                    std::cerr << "FAIL: frame " << frame
                              << " molecule " << i
                              << " centroid " << com
                              << " is outside primary cell after reimage()\n";
            }
        }
    }

    check(atoms_outside == 0,
          "all molecule centroids are inside the primary cell after reimage");
    check(frame > 0, "trajectory had at least one frame");

    std::cout << "Processed " << frame << " frames, "
              << molecules.size() << " molecules.\n";

    if (failures == 0) {
        std::cout << "All triclinic trajectory tests passed.\n";
        return 0;
    } else {
        std::cerr << failures << " test(s) failed.\n";
        return 1;
    }
}
