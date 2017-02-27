#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import argparse
from sklearn.metrics import roc_curve, roc_auc_score
import seaborn as sns

plt.style.use('seaborn-white')
vina_color='#CCBB44'

def multi_rocplot(infiles):
    """
    Read in true classes and scores from as many files as you want. Assumes the
    first two columns contain the true class and the classifier score. 
    """
    count = 0
    for file in infiles:
        if count == 0:
            contents = np.genfromtxt(file, usecols=(0,1), delimiter=" ")
            data = np.empty((contents.shape[0],2,len(infiles)))
            data[:,:,count] = contents
        else:
            data[:,:,count] = np.genfromtxt(file, usecols=(0,1), delimiter=" ")
        count += 1
    return data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generates a combined roc \
    plot from input files containing cols in (score, label) order.')
    parser.add_argument('-i', '--infiles', nargs='+', required=True,
            help='Input files containing classes and classifier output scores')
    parser.add_argument('-o', '--outname', default='multi_mean_rocplot',
    help='Optionally specify the output filename for the plot')
    args = parser.parse_args()

    # CSAR, DUD-E, 2:1
    cnn_colors = ['#332288', '#4477aa', '#88ccee']
    sns.set_palette([vina_color] + cnn_colors)
    zorders = [3, 5, 2, 4, 1]

    data = multi_rocplot(args.infiles)
    fig,ax = plt.subplots()
    aucs = []
    zorder_count = 0
    for i in range(data.shape[2]):
        fpr,tpr,_ = roc_curve(data[:,1,i],data[:,0,i])
        auc = roc_auc_score(data[:,1,i],data[:,0,i])
        aucs.append(auc)
        ax.plot(fpr,tpr,label=args.infiles[i].split('.')[0].split('_')[0]+'(AUC=%.3f)' %auc,
                linewidth=4, solid_joinstyle='miter', zorder=zorders[zorder_count])
        zorder_count += 1

    box = ax.get_position()
    handles,labels = ax.get_legend_handles_labels()
    indices = [i[0] for i in sorted(enumerate(aucs), key=lambda x:x[1],
        reverse=True)]
    new_handles = [handles[i] for i in indices]
    new_labels = [labels[i] for i in indices]
    ax.legend(new_handles, new_labels, loc='lower right', frameon=True,
            fontsize=12)
    ax.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--',
            zorder=zorders[zorder_count])
    ax.set_xlabel('False Positive Rate', fontsize=22)
    ax.set_ylabel('True Positive Rate', fontsize=22)
    ax.set_aspect('equal')
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xticks(ax.get_xticks().tolist()[1:])
    ax.set_yticks(ax.get_yticks().tolist()[1:])
    fig.savefig(args.outname+'.pdf', bbox_inches='tight')
