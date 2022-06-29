import loos
from loos import pyloos
from loos.pyloos import options
import argparse
from pathlib import Path
​
from sys import exit
​
fullhelp = """
XXX
"""


def list_difference_by_name(ag1, ag2):
    unmatched = []
    if len(ag1) == len(ag2):
        for i, a in enumerate(ag1):
            if a.name() != ag2[i].name():
                unmatched.append(i)
    else:
        longer_len = max(len(ag1), len(ag2))
        shorter_len = min(len(ag1), len(ag2))
        for i in range(shorter_len):
            if ag1[i].name() != ag2[i].name():
                unmatched.append(i)
        for j in range(shorter_len, longer_len):
            unmatched.append(j)
    return unmatched




​
​
def boolprint(*args, **kwargs):
    if args[0]:
        print(*args[1:], **kwargs)
​
​
p = argparse.ArgumentParser()
​
if __name__ == '__main__':
    lo = options.LoosOptions("Find a frame based on RMSD to provided structure",
                             fullhelp)
    lo.modelSelectionOptions()
    lo.trajOptions()
​
    lo.add_argument('--reference', '-R', default=None, type=Path,
                    help='The reference to find similarity to. If none provided'
                         ', use the coordinates of the model file.')
    lo.add_argument('--cutoff', '-c', default=0.1, type=float,
                    help='The highest RMSD that should be tolerated as '
                         'corresponding to a matching frame.')
    lo.add_argument('--out-structure-prefix', '-o', default=None, type=str,
                    help='If provided, write any discovered frames to PDBs with'
                         ' this prefix. The file(s) will be numbered in the '
                         'order that they are found.')
    lo.add_argument('--first-only', '-f', default=False,
                    action=argparse.BooleanOptionalAction,
                    help='Stop looking after the first hit.')
    lo.add_argument('--refsel', default='all', type=str,
                    help='If provided, subset the reference input by this selection.')
​
​
    args = lo.parse_args()
    header = lo.header()
​
    system = loos.createSystem(args.model)
    subset = loos.selectAtoms(system, args.refsel)
    if args.reference:
        ref = loos.selectAtoms(loos.createSystem(args.reference),
                               args.refsel)
    else:
        ref = subset.copy()
​

    diff_inds = list_difference_by_name(ref, subset)
    if len(diff_inds) > 0:
        print('The Atomic Groups selected as reference and subset are different based on their length, or the one-to-one association of their atom names. The following atoms were found to be different in this way:')
        print(*diff_inds)
        exit(2)

    found_frame = False
​
    # Stop at first hit
    if args.first_only:
        if args.out_structure_prefix:
            for traj_ix, traj_name in enumerate(args.traj):
                traj = pyloos.Trajectory(traj_name, system, skip=args.skip,
                                         stride=args.stride,
                                         subset=args.selection)
                for fra_ix, frame in enumerate(traj):
                    # perform superposition
                    ref.alignOnto(subset)
                    # compute RMSD
                    rmsd = ref.rmsd(subset)
                    if rmsd < args.cutoff:
                        print(traj_name, traj_ix, fra_ix, rmsd)
                        pdb = loos.PDB.fromAtomicGroup(system)
                        pdb.remarks.add('trajectory: {}'.format(traj_name))
                        pdb.remarks.add('frame: {}'.format(fra_ix))
                        pdb.remarks.add('rmsd: {}'.format(rmsd))
                        with open(args.out_structure_prefix + "-0.pdb", 'w') as f:
                            f.write(str(pdb))
                        exit(0)
        else:
            for traj_ix, traj_name in enumerate(args.traj):
                traj = pyloos.Trajectory(traj_name, system, skip=args.skip,
                                         stride=args.stride, subset=args.selection)
                for fra_ix, frame in enumerate(traj):
                    # perform superposition
                    ref.alignOnto(subset)
                    # compute RMSD
                    rmsd = ref.rmsd(subset)
                    boolprint(rmsd < args.cutoff, traj_name, traj_ix, fra_ix, rmsd)
                    exit(0)
    else:
        if args.out_structure_prefix:
            for traj_ix, traj_name in enumerate(args.traj):
                traj = pyloos.Trajectory(traj_name, system, skip=args.skip,
                                         stride=args.stride,
                                         subset=args.selection)
                for fra_ix, frame in enumerate(traj):
                    # perform superposition
                    ref.alignOnto(subset)
                    # compute RMSD
                    rmsd = ref.rmsd(subset)
                    if rmsd < args.cutoff:
                        found_frame = True
                        print(traj_name, traj_ix, fra_ix, rmsd)
                        pdb = loos.PDB.fromAtomicGroup(system)
                        pdb.remarks.add('trajectory: {}'.format(traj_name))
                        pdb.remarks.add('frame: {}'.format(fra_ix))
                        pdb.remarks.add('rmsd: {}'.format(rmsd))
                        with open(args.out_structure_prefix + "-0.pdb", 'w') as f:
                            f.write(str(pdb))
        else:
            for traj_ix, traj_name in enumerate(args.traj):
                traj = pyloos.Trajectory(traj_name, system, skip=args.skip,
                                         stride=args.stride, subset=args.selection)
                for fra_ix, frame in enumerate(traj):
                    # perform superposition
                    ref.alignOnto(subset)
                    # compute RMSD
                    rmsd = ref.rmsd(subset)
                    if rmsd < args.cutoff:
                        found_frame = True
                        print(traj_name, traj_ix, fra_ix, rmsd)
​
    if found_frame:
        exit(0)
    else:
        exit(1)