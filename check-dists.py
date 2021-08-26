import loos
from loos import pyloos
import numpy as np

sys = loos.createSystem('sucrose-test/SUCROSE_topology.prmtop')
trj = pyloos.Trajectory('sucrose-test/SUCROSE_trajectory.crd', sys)
glucose_h1 = loos.selectAtoms(sys, 'name == "H1" && resname == "0GA"')
fructose_h11 = loos.selectAtoms(sys, 'name == "H11" && resname == "2CU"')
fructose_h12 = loos.selectAtoms(sys, 'name == "H12" && resname == "2CU"')
dists = []
for i, frame in enumerate(trj):
    d1 = np.linalg.norm(glucose_h1.getCoords() - fructose_h11.getCoords())
    d2 = np.linalg.norm(glucose_h1.getCoords() - fructose_h12.getCoords())
    dists.append(np.array([i, d1, d2]))

header = 'frame H1-H11 H1-H12'
np.savetxt('chkdists.out', dists, header=header)