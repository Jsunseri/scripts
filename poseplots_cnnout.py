#!/usr/bin/python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import seaborn as sns

plt.style.use('seaborn-white')
sns.set_palette(sns.color_palette("Set1", n_colors=12, desat=.5))

def toprank_rmsd(data, fnames, targets, combined, outname, target_labels):
    rmsds = []
    best_possible = []
    for file in fnames:
        for target in targets:
            idx = data.loc[(data["file"] == file) & (data["rec"] == target),
            "score"].idxmax()
            rmsd = data.loc[(data["file"] == combined) & (data["lig"] ==
                data.loc[idx]["lig"]), "rmsd"]
            rmsds.append({'file' : file, 'target' : target, 'rmsd' : rmsd})
    rmsds = pd.DataFrame(rmsds)
    for target in targets:
        best_possible.append(data.loc[(data["file"] == combined) & (data["rec"] ==
            target),"rmsd"].min())
    fig,ax = plt.subplots()
    x = np.arange(1, len(targets)+1)
    ax.plot(x, best_possible, marker='o', alpha=0.5, ls='None', label='Best\
            RMSD (mean=%.2f)' % np.array(best_possible).mean())
    for file in fnames:
        ax.plot(x,rmsds.loc[rmsds["file"] == file, "rmsd"], marker='o',
                alpha=0.5, ls='None', label=file+'(mean=%.2f)' %
                (rmsds.loc[rmsds["file"] ==
                    file,"rmsd"]).as_matrix().astype(np.float32).mean())
    ax.legend(loc='best', bbox_to_anchor=(1,0.5))
    ax.set_xlabel('Target')
    ax.set_ylabel('RMSD of top pose')
    ax.set_xticks(x)
    ax.set_xticklabels(target_labels, rotation=90)
    fig.savefig(outname+'.png', bbox_inches='tight')

def cumul_best_mean(data, fnames, targets, combined, outname, maxrank):
    best_by_rank = np.empty((maxrank, len(targets), len(fnames)))
    for file in fnames:
        for target in targets:
            sorted = data.loc[(data["file"] == file) & (data["rec"] ==
                target)].sort_values(["score"], ascending=False)
            for i in range(maxrank):
                best_by_rank[i,np.where(targets==target),fnames.index(file)] =\
                sorted['rmsd'][:i+1].min()
    fig,ax = plt.subplots()
    x = np.arange(1, maxrank+1)
    for file in fnames:
        if file == "Vina":
            ax.plot(x, best_by_rank[:,:,fnames.index(file)].mean(1),
                    label=file, color="black")
        else:
            ax.plot(x, best_by_rank[:,:,fnames.index(file)].mean(1), label=file)
    ax.legend(loc='best', bbox_to_anchor=(1,0.5))
    ax.set_xlabel('Rank')
    ax.set_ylabel('Mean Cumulative Best RMSD')
    ax.set_xticks(x)
    fig.savefig(outname+'_cumulbest.pdf', bbox_inches='tight')

def per_target_top1(data, fnames, targets, outname):
    overall_best = []
    for target in targets:
        print data.iloc[data.loc[data["rec"] == target,"rmsd"].argmin()]["rec"]
    # indices = overall_best[:,1].argsort()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create plots summarizing\
    pose ranking for a classifier. At a minimum, need a file with scores and\
    RMSDs; add extra files that contain scores for the same poses to summarize\
    performance for several classifiers")
    parser.add_argument('-c','--combined',required=True, help='File that\
    contains combined score and RMSD info')
    parser.add_argument('-a','--additional', nargs='+', required=False,
            help='all the other files you want to include in the summary\
            plots; these don\'t need RMSD info')
    parser.add_argument('-o','--outname', default='top_pose', help='Output\
    file name; defaults to top_pose.png')
    args = parser.parse_args()

data = pd.read_csv(args.combined, sep=" ", header=None, names=["score",
"label", "rec", "lig", "rmsd"])
data["file"] = args.combined.split('.')[0]
data["score"] = data["score"]
fnames = [args.combined.split('.')[0]]
targets = data.rec.unique()
target_labels = [x.split('/')[-2] for x in targets]

if args.additional:
    fnames = fnames + [x.split('.')[0] for x in args.additional]
    for file in args.additional:
        df = pd.read_csv(file, sep=" ", header=None, names=["score", "label",
            "rec", "lig", "rmsd"])
        df["file"] = file.split('.')[0]
        data = pd.concat([data,df], ignore_index=True)

# cumul_best_mean(data, fnames, targets, args.combined.split('.')[0],
        # args.outname, 10)
# toprank_rmsd(data, fnames, targets, args.combined.split('.')[0], args.outname,
        # target_labels)
per_target_top1(data, fnames, targets, args.outname)
