#!/usr/bin/env python3

import sys
import loos
from numpy.random import default_rng

cmdline = " ".join(sys.argv)

system_file = sys.argv[1]
rot_size = float(sys.argv[2])
dcd_file = sys.argv[3]
num_steps = int(sys.argv[4])

dcd = loos.DCDWriter(dcd_file)
dcd.addTitle(cmdline)


rg = default_rng()

system = loos.createSystem(system_file)
max_idx = len(system) - 1

axis = loos.GCoord()
for i in range(num_steps):
    # generate a random axis and displacement
    axis.random()
    rot = rg.uniform(-rot_size, rot_size)
    rotation_center = system[int(rg.integers(0, max_idx))].coords()

    system.rotate(axis, rot, rotation_center)

    system.centerAtOrigin()

    dcd.writeFrame(system)