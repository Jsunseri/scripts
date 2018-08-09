#! /usr/bin/env python
import os,sys,math
import collections
from sklearn.metrics import roc_curve, roc_auc_score
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.lines as mlines

plt.style.use('seaborn-white')
SMALL_SIZE=14
MEDIUM_SIZE=16
BIGGER_SIZE=24
SUBFIG_SIZE=30

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)

paper_palettes = {}
paper_palettes['Vina'] = '#000000' #the original CNN paper used ccbb44
paper_palettes['CSAR'] = '#332288'
paper_palettes['DUD-E'] = '#4477aa'
paper_palettes['2:1'] = '#88ccee'
paper_palettes['CNN Affinity Rescore'] = '#6da4c0'
paper_palettes['CNN Affinity Refine'] = '#332288'
paper_palettes['CNN Scoring Rescore'] = '#ffd91c'
paper_palettes['CNN Scoring Refine'] = '#877b25'
paper_palettes['Experiment'] = '#498540'
swarm_markers = ['^', '>', '*', 's', 'o', '<']
blue_star = mlines.Line2D([], [], color='black', marker='^', linestyle='None',
                                  markersize=10, label='Blue stars')
red_square = mlines.Line2D([], [], color='red', marker='s', linestyle='None',
                                  markersize=10, label='Red squares')
purple_triangle = mlines.Line2D([], [], color='purple', marker='^',
        linestyle='None',
                                  markersize=10, label='Purple triangles')

def calc_auc(target_and_method, target_predictions):
	y_true=[]
	y_score=[]
	for i,item in enumerate(target_predictions):
	    try:
                label = float(item[0])
                score = float(item[1])
	        y_true.append(label)
	        y_score.append(score)
	    except Exception as e:
	        print 'Error: %d %f %s\n'%(label, score, target_and_method[i])
	        continue
        fpr,tpr,_ = roc_curve(y_true,y_score)
	return roc_auc_score(y_true,y_score),fpr,tpr

def mean_auc(data, methods, targets, args):
        #use this palette if the methods don't correspond to methods used in
        #any of the old papers, which have associated colors
        backup_palette = sns.color_palette("Set1", n_colors=len(methods), desat=.5)

        #overall_auc tracks and reports the average AUC across all targets
        #bytarget tracks and reports the AUC for each target for each method
        overall_auc = {}
        bytarget = {}
        for method in methods:
	    bytarget[method] = open('%s_%s'%(method,args.outprefix),'w')
        if args.make_boxplot:
            boxplot_dat = []
        total_plots = len(targets)
        if total_plots > 9:
            mpl.rcParams['xtick.labelsize'] = 8
            mpl.rcParams['ytick.labelsize'] = 8
        grid_width = int(math.ceil(math.sqrt(total_plots)))
        grid_length = int(math.ceil(float(total_plots)/grid_width))
        fig,fig_ax = plt.subplots(figsize=(8,8))

        #if there is only one output plot, print AUCs on the legend; otherwise
        #make a single shared legend for all plots without AUCs
        legend_dict = {}
        num_lines = {}
	for t,l in data.iteritems():
		auc,fpr,tpr = calc_auc(t,l)
                bytarget[t[1]].write('%s %.3f\n' %(t[0],auc))
                if args.make_boxplot:
                    boxplot_dat.append({'Method' : t[1], 'AUC' : auc, 'Target'
                        : t[0]})
                plot_num = targets.index(t[0])
                if t[1] not in overall_auc.keys():
                    overall_auc[t[1]] = 0
                overall_auc[t[1]] += float(auc)
                if plot_num not in legend_dict.keys(): 
                    num_lines[plot_num] = 0 
                    ax = plt.subplot2grid((grid_length,grid_width),
                            (plot_num / grid_width, plot_num % grid_width))
                    ax.set_aspect('equal')
                    for tick in ax.get_xticklabels():
                        tick.set_rotation(-90)
                    legend_dict[plot_num] = ax
                method = t[1]
                color = paper_palettes[method] if method in paper_palettes else backup_palette[methods.index(method)]
                label = '%s, AUC=%0.2f' % (t[1], auc) if total_plots == 1 else t[1]
                legend_dict[plot_num].plot(fpr, tpr, color=color, label=label,
                        zorder=2) 
                num_lines[plot_num] += 1
                if int(plot_num) / grid_width == grid_length-1:
                    legend_dict[plot_num].set_xlabel('False Positive Rate')
                if plot_num % grid_width == 0:
                    legend_dict[plot_num].set_ylabel('True Positive Rate')
                if num_lines[plot_num] == len(methods):
                    legend_dict[plot_num].set_title('Target %s' % (t[0]))
                    if total_plots == 1:
                        legend_dict[plot_num].legend(loc='center left',
                                bbox_to_anchor=(1,0.5))

        #if we have multiple subplots, make a shared legend; constrain CNN
        #paper order to follow line order
        handles,labels = legend_dict[0].get_legend_handles_labels()
        if args.color_scheme == 'cnn':
            shortlabels = [x.split(',')[0] for x in labels]
            indices = []
            cnnpaper_order = ['DUD-E', 'Vina', '2:1', 'CSAR']
            for label in cnnpaper_order:
                indices.append(shortlabels.index(label))
            handles = [handles[i] for i in indices]
            labels = [labels[i] for i in indices]
            box = legend_dict[total_plots-2].get_position()
            fig.legend(handles, labels, loc=(box.x0+0.465, box.y0-0.035), frameon=True)
        elif args.color_scheme == 'd3r':
            indices = []
            cnnpaper_order = ['CNN Affinity Rescore', 'CNN Affinity Refine',
                    'CNN Scoring Rescore', 'CNN Scoring Refine', 'Vina']
            for label in cnnpaper_order:
                indices.append(labels.index(label))
            handles = [handles[i] for i in indices]
            labels = [labels[i] for i in indices]
            fig.legend(handles, labels, loc='lower center', ncol=3, frameon=True)
        else:
            box = legend_dict[total_plots-1].get_position()
            if grid_length == grid_width:
                #put the legend underneath
                xloc = box.x0 + 1
                yloc = box.y0 + 5
            else:
                xloc = box.x0 + 0.465
                yloc = box.y0 - 0.035
            legend = fig.legend(handles, labels, loc='best', frameon=True)

        #add in line showing random performance
        for plot_num in legend_dict:
            legend_dict[plot_num].plot([0, 1], [0, 1], color='gray', lw=2,
               linestyle='--', zorder=1)
        if args.color_scheme == 'cnn':
            fig.subplots_adjust(wspace=0.05, hspace=0.5)
        else:
            #TODO: currently doing this via manual tuning - can it be
            #automated? if you are using this, you may need to fiddle with
            #these numbers and be aware that if you try to use tight_layout it
            #seems to override whatever you do here
            fig.subplots_adjust(hspace=-0.4, wspace=0.45)
        for method in overall_auc.keys():
	    overall_auc[method] /= total_plots
	    # outfile.write('%s, AUC %.2f\n' 
                    # % (method.split('_')[0],overall_auc[method]))
        for method in methods:
            bytarget[method].close()
        fig.savefig(args.plotname+'.pdf', bbox_inches='tight')

        #now do boxplot
        fig,ax = plt.subplots()
        if args.make_boxplot:
            if args.color_scheme == 'd3r' or args.color_scheme == 'cnn':
                palette = paper_palettes
            else:
                palette = backup_palette
            boxplot_df = pd.DataFrame(boxplot_dat)
            #if we're doing the d3r paper figs, David wants a different marker
            #style for each target because having per-target ROC plots isn't
            #enough for him...
            if args.color_scheme == 'd3r':
                cnnpaper_order = ['CNN Affinity Rescore', 'CNN Affinity Refine',
                        'CNN Scoring Rescore', 'CNN Scoring Refine', 'Vina']
                leghands = []
                for marker_id,target in enumerate(targets):
                    marker = swarm_markers[marker_id]
                    if marker == '+' or marker == '*':
                        mew = 10
                        size = 15
                    else:
                        mew = 0.5
                        size = 10
                    sns.stripplot(x='Method', y='AUC',
                            data=boxplot_df[boxplot_df['Target']==target],
                            split=True, edgecolor='black', size=size, linewidth=0,
                            linewidths=mew, jitter = True,
                            palette=palette, marker=marker,
                            order=cnnpaper_order, ax=ax)
                    leghands.append(mlines.Line2D([], [], color='black',
                        fillstyle='none', marker=marker, linestyle='None',
                        mew=1,
                        markersize=size, label=target))
                # ax.legend(handles=leghands, bbox_to_anchor=(1,-0.2),
                        # frameon=True)
                ax.legend(handles=leghands, loc='lower left', ncol=2, 
                        frameon=True)
            else:
                sns.swarmplot(x='Method', y='AUC',
                        data=boxplot_df, split=True, edgecolor='black', size=7,
                        linewidth=0, palette = palette, ax=ax)
            sns.boxplot(x='Method', y='AUC', data=boxplot_df,
                    color='white', order=cnnpaper_order, ax=ax)
            ax.set_ylabel('AUCs')
            ax.set_xlabel('')
            ax.set(ylim=(0,1.1))
            xlims = ax.get_xlim()
            ax.plot([xlims[0],xlims[1]],[0.5, 0.5], linestyle='--', color='gray',
                    lw=2, zorder=1, alpha=0.5)
            for tick in ax.get_xticklabels():
                tick.set_rotation(45)
            if args.color_scheme == 'd3r':
                labels = ax.get_xticklabels()
                labels = [label.get_text().split() for label in labels]
                labels = ['%s %s\n%s'%(x[0],x[1],x[2]) if not x[0] == 'Vina' else x[0] for x in labels]
                ax.set_xticklabels(labels)
            fig.savefig(args.plotname+'_boxplot.pdf', bbox_inches='tight')
	return overall_auc

if __name__=='__main__':
	parser = ArgumentParser(description='Calculate AUC by target for multiple methods')
	parser.add_argument('-p','--predictions', nargs='*', default=[], 
                help='files of predictions, formatted LABELS PREDICTIONS TARGET METHOD where METHOD is stylized as desired for output, with spaces replaced with underscores')
	parser.add_argument('-o','--outprefix',type=str,default='bytarget',help='prefix for all output files')
        parser.add_argument('-n','--plotname', default='roc_curve_bytarget',
                help='Output filename for by-target ROC curves')
        parser.add_argument('-make_boxplot', action='store_true',
                default=False, help='Make a boxplot of the by-target AUC values \
associated with each method')
        parser.add_argument('-color_scheme', required=False, choices=['cnn', 'd3r'], 
            help='Specify color scheme, choosing from the one used in the CNN \
paper or the one used in the 2018 D3R paper; if used, the \
predictions files must have names indicating the correct methods to \
use those color schemes')
	args= parser.parse_args()
	
        data = {}
        methods,targets = [],[]
        if args.color_scheme == 'cnn':
            paper_palettes['Vina'] = '#CCBB44'
	for file in args.predictions:
	    for line in open(file,'r'):
                contents = line.split()
                target = contents[2].replace('_', ' ')
                method = contents[3].replace('_', ' ')
                if target not in targets:
                    targets.append(target)
                if method not in methods:
                    methods.append(method)
                this_key = (target,method)
                if this_key not in data:
                    data[this_key] = []
                data[this_key].append((contents[0],contents[1]))
	auc = mean_auc(data, methods, targets, args)
        print auc
