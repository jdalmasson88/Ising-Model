import random
import math
import os
import argparse
import time
import glob
import numpy as np
import cPickle as pkl
from matplotlib import pylab as plt

# to run the code call python cluster_ising.py followed by
# the number of of dimension of the lattice and the number
# of sites on the edge. For example:
# python cluster_ising.py 2 100 to create a 100x100 lattice
parser = argparse.ArgumentParser()

parser.add_argument(
    'n_dim', help='number of dimensions of the hypercubic lattice')
parser.add_argument('length', help='number of sites per edge')

args = parser.parse_args()
n_dim = int(args.n_dim)
n_dim_file = 'pickle' + args.n_dim + 'd_multip.p'
length = int(args.length)
N = pow(length, n_dim)

print(n_dim_file, length)


def lattice(n_dim, Len):
    H = [pow(Len, i) for i in xrange(0, n_dim + 2)]
    return {i: np.ravel([((i // H[j + 1]) * H[j + 1] + (i + H[j]) % H[j + 1], (i // H[j + 1]) * H[j + 1] + (i - H[j]) % H[j + 1]) for j in xrange(0, n_dim)]) for i in xrange(H[n_dim])}


# call to the function creating the lattice and the nn structure
# loop going through the different temperatures
# for a given lattice dimension
n_lattice = lattice(n_dim, length)

T = [np.linspace(1.5, 3, 26), np.linspace(4.03, 5, 10),
     np.linspace(5.5, 7.5, 21), np.linspace(7, 9.5, 26)]
av_M = []

for i, e in enumerate(T[n_dim - 2]):
    M = []
    for smp in xrange(0, 4):
        start_time = time.time()
        # acceptance probability for the cluster expansion
        p = 1.0 - math.exp(-2.0 / e)
        # number of Monte Carlo steps
        nsteps = 10000
        # +/-1 is randomly assigned to every site
        S = [random.choice([1, -1]) for k in xrange(N)]
        # beginning of the Wolff algorithm (see the article for reference)
        for step in xrange(nsteps):
            k = random.randint(0, N - 1)
            Pocket, Cluster = [k], [k]
            while Pocket != []:
                j = random.choice(Pocket)
                for l in n_lattice[j]:
                    if S[l] == S[j] and l not in Cluster and random.uniform(0.0, 1.0) < p:
                        Pocket.append(l)
                        Cluster.append(l)
                Pocket.remove(j)
            S = np.asarray(S) * (-1)
            stop_loop = time.time()
            # this condition is commented if the first method is followed
            # if stop_loop-start_time > 2000:
            #   print('stopping %i step'%step)
            #   break
        # magnetization array created
        M.append(abs(np.mean(S)))
    # the output of the simulation is saved in av_M (commented for
    # first method)
    av_M.append([e, np.mean(M), np.std(M), S])
    # commented for second method
    # av_M.append([e,S])

# the output is saved in a pickle file
with open(n_dim_file, 'wb') as wfp:
    pkl.dump(av_M, wfp)
