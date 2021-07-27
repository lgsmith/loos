import loos
from loos import pyloos
import numpy as np

sys = loos.createSystem('sucrose-test/SUCROSE_topology.prmtop')
trj = pyloos.Trajectory('sucrose-test/SUCROSE_trajectory.crd', sys)
glucose_h1 = loos.selectAtoms(sys, 'name == "H1" && resname == "0GA"')
fructose_h11 = loos.selectAtoms(sys, 'name == "H11" && resname == "2CU"')
fructose_h12 = loos.selectAtoms(sys, 'name == "H12" && resname == "2CU"')
dists = []
for frame in trj:
    d1 = glucose_h1.dist(fructose_h11)
    d2 = glucose_h1.dist(fructose_h12)
    print(d1, d2)