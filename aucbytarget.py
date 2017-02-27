#! /usr/bin/env python
import os,sys,math
import collections
from sklearn.metrics import roc_curve, roc_auc_score
from argparse import ArgumentParser
import numpy as np
import bootstrap
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

plt.style.use('seaborn-white')
vina_color='#CCBB44'
# CSAR, DUD-E, 2:1
cnn_colors = ['#332288', '#4477aa', '#88ccee']
sns.set_palette([vina_color] + cnn_colors)

def calc_auc(target_predictions,use_bootstrap):
	if use_bootstrap: 
		_,auc,_= bootstrap.bootstrap(target_predictions,1000)
		return auc
	y_true=[]
	y_score=[]
	for line in target_predictions:
		try:
			values= line.split(" ")
			if float(values[0]) == 1 or float(values[0]) ==0:
				y_true.append(float(values[0]))
				y_score.append(float(values[1]))
			elif float(values[1]) ==1 or float(values[1])==0:
				y_true.append(float(values[1]))
				y_score.append(float(values[0]))
		except Exception as e:
			print 'Error: %s'%line
			continue
        fpr,tpr,_ = roc_curve(y_true,y_score)
	return roc_auc_score(y_true,y_score),fpr,tpr

def mean_auc(ordered_dict,methods,target_list,args):
        overall_auc = {}
        bytarget_dict = {}
        for method in methods:
	    bytarget_dict[method] = open('%s_%s'%(method,args.outprefix),'w')
        if args.make_boxplot:
            boxplot_dat = []
        total_plots = len(target_list)
        if total_plots > 9:
            mpl.rcParams['xtick.labelsize'] = 8
            mpl.rcParams['ytick.labelsize'] = 8
        # sns.set_palette(sns.color_palette("Set1", n_colors=len(methods), desat=.5))
        grid_width = int(math.ceil(math.sqrt(total_plots)))
        grid_length = int(math.ceil(float(total_plots)/grid_width))
        fig,fig_ax = plt.subplots(figsize=(8,8))
        legend_dict = {}
        num_lines = {}
	for t,l in ordered_dict.items():
		auc,fpr,tpr = calc_auc(l,args.bootstrap)
                bytarget_dict[t[1]].write('%s %.3f\n'%(t[0],auc))
                if args.make_boxplot:
                    boxplot_dat.append({'Method' : t[1].split('_')[0], 'AUC' : auc})
                plot_num = target_list.index(t[0])
                if t[1] not in overall_auc.keys():
                    overall_auc[t[1]] = 0
                overall_auc[t[1]] += float(auc)
                if plot_num not in legend_dict.keys(): 
                    num_lines[plot_num] = 0 
                    ax = plt.subplot2grid((grid_length,grid_width),
                            (plot_num / grid_width, plot_num % grid_width))
                    ax.set_aspect('equal')
                    legend_dict[plot_num] = ax
                legend_dict[plot_num].plot(fpr,tpr,#color=sns.color_palette()[methods.index(t[1])-1],
                        label=t[1].split('_')[0], zorder=2) #'%s, AUC=%0.2f' % (t[1].split('_')[0],auc))
                num_lines[plot_num] += 1
                if int(plot_num) / grid_width == grid_length-1:
                    legend_dict[plot_num].set_xlabel('False Positive Rate')
                if plot_num % grid_width == 0:
                    legend_dict[plot_num].set_ylabel('True Positive Rate')
                if num_lines[plot_num] == len(methods):
                    legend_dict[plot_num].set_title('Target %s' % (t[0]))
                    # legend_dict[plot_num].legend(loc='center left', bbox_to_anchor=(1,0.5), 
                            # fontsize='xx-small')
		if args.create_files:
			pfile= open('%s_predictions'%t[1], 'w')
			for line in l:
				pfile.write('%s\n'%line)
			print '%s_predictions created'%t[1]
			pfile.close()
        handle_list = []
        label_list = []
        for plot in legend_dict.keys():
            handles,labels = legend_dict[plot].get_legend_handles_labels()
            handle_list = handle_list + handles
            label_list = label_list + labels
        #TODO: remove the next few lines after CNN paper figs are done...
        indices = [2,0,3,1]
        final_handles = [handle_list[i] for i in indices]
        final_labels = [label_list[i] for i in indices]
        box = legend_dict[total_plots-2].get_position()
        # fig.legend(handle_list[0:len(args.predictions)],
                # label_list[0:len(args.predictions)],
                # loc=(box.x0-0.24,box.y0-0.105),
                # fontsize=14, frameon=True, ncol=4)
        fig.legend(final_handles, final_labels, loc=(box.x0+0.465, box.y0-0.035),
                fontsize=14, frameon=True)
        for plot_num in legend_dict:
            legend_dict[plot_num].plot([0, 1], [0, 1], color='gray', lw=2,
               linestyle='--', zorder=1)
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.05,hspace=0.5)
        for method in overall_auc.keys():
	    overall_auc[method] /= total_plots
	    # outfile.write('%s, AUC %.2f\n' 
                    # % (method.split('_')[0],overall_auc[method]))
        for method in methods:
            bytarget_dict[method].close()
        fig.savefig(args.plotname+'.pdf', bbox_inches='tight')
        fig,ax = plt.subplots()
        if args.make_boxplot:
            boxplot_df = pd.DataFrame(boxplot_dat)
            ax = sns.swarmplot(x='Method', y='AUC',
                    data=boxplot_df, split=True, edgecolor='black', linewidth=1)
            ax = sns.boxplot(x='Method', y='AUC', data=boxplot_df)
            handles,labels = ax.get_legend_handles_labels()
            ax.legend(handles[0:len(args.predictions)],
                    labels[0:len(args.predictions)], loc='best', bbox_to_anchor=(1,0.5))
            ax.set_ylabel('AUCs')
            fig.savefig(args.plotname+'_boxplot.pdf', bbox_inches='tight')
	return overall_auc

if __name__=='__main__':
	parser = ArgumentParser(description='Calculate AUC by targets')
	parser.add_argument('-p','--predictions', nargs='*',default=[],help='files of predictions')
	parser.add_argument('-o','--outprefix',type=str,default='bytarget',help='prefix for all output files')
	parser.add_argument('-create_files', action='store_true',default=False, help='create files for each target')
	parser.add_argument('-b','--bootstrap', action='store_true',default=False,help='use bootstrapping when calculating AUCs')
        parser.add_argument('-n','--plotname', default='roc_curve_bytarget',
                help='Output filename for by-target ROC curves')
        parser.add_argument('-make_boxplot', action='store_true',
                default=False, help='Make a boxplot of the AUC values \
                associated with each method')
	args= parser.parse_args()
	
	 
	ordered_dict = collections.OrderedDict()
        methods = []
        target_list = []
	for file in args.predictions:
		for line in open(file,'r'):
			target = line.split(" ")[2]
                        temp_target = os.path.basename(os.path.dirname(target))
                        #for swapping the decoys associated with the target
                        #shared between MUV and ChEMBL
                        # if "chembl" in file.lower() and temp_target == "466":
                            # temp_target = "11631"
                        # elif "muv" in file.lower() and temp_target == "11631":
                            # temp_target = "466"
        		target = (temp_target,
                                os.path.splitext(os.path.basename(file))[0])
        		if target in ordered_dict: ordered_dict[target].append(line)
        		else:
                		#print target 
                		ordered_dict[target]=[line]
                        if target[0] not in target_list: target_list.append(target[0])
                        if target[1] not in methods: methods.append(target[1])
	auc = mean_auc(ordered_dict,methods,target_list,args)
	print auc
	
