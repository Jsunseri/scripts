#!/usr/bin/python

import numpy as np
import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

def readin_dat(dat, dir, traj, fname='mmpbsa.dat'):
    # read data in from the MMPBSA.py output
    with open(fname) as file:
        i = 0
        for line in file:
            if line.startswith('DELTA TOTAL'):
                if i == 0:
                    dat.append({'System' : dir, 'Energy' : line.split()[2],
                        'Traj' : traj, 'Solvent' : 'g'})
                    i = 1
                else:
                    dat.append({'System' : dir, 'Energy' : line.split()[2],
                        'Traj' : traj, 'Solvent' : 'p'})

# Generate a data frame

_, dirs, _ = os.walk('.').next()

dat = []
means = []
trajs = ['1traj', '3traj']
solvent = ['p', 'g']

# TODO: generalize. Currently writing this for profilin, but this should be
# general-purpose

for dir in dirs:
    for i in range(1,51):
        for traj in trajs:
            os.chdir('/mnt/dvorak/net/pulsar/home/koes/jss97/profilin/' + dir + '/run' + str(i) +
                    '_shorter/' + traj + '_enthalpy')
            readin_dat(dat, dir, traj)

os.chdir('/home/jocelynsunseri/Documents/profilin/')
dframe = pd.DataFrame(dat)

# Bootstrap each subset and compute the standard deviation of the resulting
# sample means; also plot the resulting distributions for reference

for dir in dirs:
    for traj in trajs:
        for sol in solvent:
            sample = dframe.loc[(dframe['System'] == dir) & (dframe['Traj'] ==
                traj) & (dframe['Solvent'] ==
                    sol)].as_matrix(columns=['Energy'])
            sample = sample.ravel()
            sample = sample.astype(np.float32, copy=False)
            for i in range(0,10000):
                boot = np.random.choice(sample.ravel(), sample.size)
                means.append({'System' : dir, 'Traj' : traj, 'Solvent' : sol,
                    'Mean' : np.mean(boot)})

dmeans = pd.DataFrame(means)

for dir in dirs:
    for traj in trajs:
        for sol in solvent:
            mean_sample = dmeans.loc[(dmeans['System'] == dir) &
                    (dmeans['Traj'] == traj) & (dmeans['Solvent'] ==
                        sol)].as_matrix(columns=['Mean'])
            mean_sample = mean_sample.ravel()
            mean_sample = mean_sample.astype(np.float32, copy=False)
            # fig = plt.figure()
            # sns.distplot(mean_sample, kde=False, rug=False)
            # fig.savefig(dir + '_' + traj + '_' + sol +
            # '_means.png',bbox_inches='tight')

            print dir, traj, sol, np.mean(mean_sample), np.std(mean_sample)
