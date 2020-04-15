import numpy as np
from scipy.signal import correlate
import loos
from loos import pyloos

modelfn = '/home/louis/gagug/gagug-tester.pdb'
trajfn = '/home/louis/gagug/gagug-tester.dcd'
selectionstr = '(name == "H8" || name == "H1\'") && resid < 3'

model = loos.createSystem(modelfn)
traj = pyloos.Trajectory(trajfn, model, subset=selectionstr)

zdist_frames = []
zs_frames = []
for frame in traj:
    zdists = []
    zs = []
    for i in range(0, len(frame)):
        zs.append(frame[i].coords()[2])
        for j in range(i+1, len(frame)):
            zdists.append(
                frame[i].coords()[2] - frame[j].coords()[2]
            )
    zdist_frames.append(np.array(zdists))
    zs_frames.append(zs)

zdist_frames = np.array(zdist_frames)
zdist_corrs = np.array([correlate(zdist_frames[:,col], zdist_frames[:,col]) for col in np.arange(zdist_frames.shape[1])])

corrz = correlate(zs_frames, zs_frames)
print(zdist_corrs)
print('just coords')
print(corrz.T[0:6])
print(corrz.T[0]-corrz.T[6])