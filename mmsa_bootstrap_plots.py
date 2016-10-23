#!/usr/bin/python

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import random
import sys

plt.style.use('seaborn-white')
palette = sns.color_palette()

data = []
jit = {}
jit["1traj"] = -1
jit["3traj"] = 1

for line in open('boot.out'):
    line_args = line.split()
    data.append({'System' : line_args[0].strip(),
        'Traj' : line_args[1].strip(), 
        'Solvent' : line_args[2].strip(), 'Energy' :
        line_args[3].strip(), 'Error' : line_args[4].strip()})

df = pd.DataFrame(data)
systems = pd.unique(df.System.ravel())

trajs = pd.unique(df.Traj.ravel())
solvents = pd.unique(df.Solvent.ravel())
x = np.arange(1, len(systems)+1)

fig,ax = plt.subplots()
for solvent in solvents:
    for traj in trajs:
        mc = palette[jit.keys().index(traj)]
        if solvent == "g":
            fc = 'white'
        else:
            fc = mc
        vals = df.loc[(df['Traj'] == traj) & (df['Solvent'] == solvent)]
        plt.errorbar(x+random.uniform(jit[traj]*.05,jit[traj]*.25),
            vals.as_matrix(columns=['Energy']).astype(np.float32, copy=False),
            yerr=vals.as_matrix(columns=['Error']).astype(np.float32,
                copy=False), marker='o',
            ecolor=mc, markerfacecolor=fc, ls='None',
            markeredgecolor=mc, markeredgewidth=0.5, label='MM%sBSA %s' %
            (solvent.upper(),traj))

ax.set_xticks(x)
ax.set_xticklabels(systems,rotation=90)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
ax.set_ylabel('Estimated Enthalpy of Binding')
plt.title(sys.argv[1])
fig.savefig(sys.argv[2],bbox_inches='tight')
