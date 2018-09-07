#!/usr/bin/env python

#generate model files and PBS scripts (using the train.py script) to do a
#random hyperparameter search for the subcube LSTM.
#sample: - LSTM hidden layer dimension
#        - number of input conv layers (and number of filters per layer)
#        - learning rate
#        - momentum
#        - weight initialization
#        - batch size
#        - subgrid dim
#        - resolution
#        - activation function
#        - loss? (if we're doing affinity prediction) (TBD)

from make_caffemodel import make_model
import numpy as np
import argparse
import math

parser = argparse.ArgumentParser(description="Generate model prototxt and job \
script for LSTM random hyperparameter optimization.")

parser.add_argument('-n', '--num_models', default=100, type=int, help='Number of files to \
generate (default=100).')
parser.add_argument('-q', '--queue', default='dept_gpu', help='Queue you want \
to submit the job to (default=dept_gpu).')
parser.add_argument('-o', '--output_prefix', default='lstm', help='Prefix for \
.model file and job script (default="lstm").')
parser.add_argument('-aff', '--affinity', default=False, action='store_true',
        help='Indicate that input data has affinity.')
parser.add_argument('-lab', '--label', default=False, action='store_true',
        help='Indicate that input data has a class label.')
parser.add_argument('-lmap', '--ligmap', default='', help='Location of \
desired nonstandard ligmap file.')
parser.add_argument('-rmap', '--recmap', default='', help='Location of \
desired nonstandard recmap file.')
parser.add_argument('-stratrec', '--stratify_receptor', default=False,
        action='store_true', help='Stratify by receptor when training.')
parser.add_argument('-p', '--prefix', default='', help='.types file prefix.')
parser.add_argument('--trainargs', required=False, nargs='+', help='Any \
additional arguments you want to pass through to the train.py script.')
args = parser.parse_args()

layerspec = {}
layerspec['MolGridData'] = {}
layerspec['MolGridData']['affinity'] = False
layerspec['MolGridData']['label'] = True

#TODO: loss?
hidden_dim_ = 2 * np.logspace(1, 3)
n_cnn_layers_ = range(6)
n_cnn_filters_ = [32, 64]
learning_rate_ = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1] 
momentum_ = 1 - np.array([1e-2, 1e-1, 1])
weight_init_ = ['constant', 'uniform', 'gaussian', 'xavier', 'positive_unitball']
constant_weight = [0.5, 0.75, 1]
gaussian_weight = np.linspace(0, 1)
xavier_weight = ['FAN_IN', 'FAN_OUT', 'AVERAGE']
batch_size_ = [10, 32, 64, 128]
subgrid_dim_ = [6, 8, 10, 12]
resolution_ = [0.25, 0.5, 0.75]
activation_function_ = ['ELU', 'ReLU']

for i in range(args.num_models):
    outname = args.output_prefix + '_' + str(i)
    #make .model file
    layerspec = {}
    layerspec['MolGridData'] = {}
    layerspec['MolGridData']['affinity'] = False
    layerspec['MolGridData']['label'] = False
    if args.affinity:
        layerspec['MolGridData']['affinity'] = True
    if args.label:
        layerspec['MolGridData']['label'] = True

    layerspec['MolGridData']['batch_size'] = np.random.choice(batch_size_)
    subgrid_dim = np.random.choice(subgrid_dim_)
    layerspec['MolGridData']['subgrid_dim'] = subgrid_dim
    resolution = np.random.choice(resolution_)
    layerspec['MolGridData']['resolution'] = resolution
    layerspec['MolGridData']['recmap'] = args.recmap
    layerspec['MolGridData']['ligmap'] = args.ligmap
    layerspec['MolGridData']['stratify_receptor'] = args.stratify_receptor
    #close to 23.5, but choose something evenly divisible even though it
    #shouldn't actually be necessary
    dimension = 23.5
    factor = round((dimension - subgrid_dim) / (subgrid_dim + resolution))
    layerspec['MolGridData']['dimension'] = factor * (subgrid_dim + resolution) + subgrid_dim

    layerspec['LSTM'] = {}
    layerspec['LSTM']['hidden_dim']= [np.random.choice(hidden_dim_)]
    layerspec['LSTM']['n_lstm_layers'] = 1 #let's not bother with varying 
                                           #this yet due to memory issues
    layerspec['LSTM']['lstm_weight_fill'] = {}
    weight_init = np.random.choice(weight_init_)
    layerspec['LSTM']['lstm_weight_fill']['type'] = weight_init
    #depending on weight init type, we need to make some more choices
    if weight_init == 'constant':
        layerspec['LSTM']['lstm_weight_fill']['value'] = np.random.choice(constant_weight)
    if weight_init == 'uniform':
        layerspec['LSTM']['lstm_weight_fill']['min'] = 0.25
        layerspec['LSTM']['lstm_weight_fill']['max'] = 0.75
    if weight_init == 'gaussian':
        layerspec['LSTM']['lstm_weight_fill']['std'] = np.random.choice(gaussian_weight)
    if weight_init == 'xavier':
        layerspec['LSTM']['lstm_weight_fill']['variance_norm'] = np.random.choice(xavier_weight)
    layerspec['LSTM']['lstm_weight_fill'] = [layerspec['LSTM']['lstm_weight_fill']]

    layerspec['Convolution'] = {}
    layerspec['Convolution']['n_cnn_layers'] = np.random.choice(n_cnn_layers_)
    layerspec['Convolution']['n_cnn_filters'] = [np.random.choice(n_cnn_filters_)]
    layerspec['Convolution']['stride'] = [1]
    layerspec['Convolution']['pad'] = [1]
    layerspec['Convolution']['kernel_size'] = [3]
    layerspec['Convolution']['activation'] = np.random.choice(activation_function_)
    layerspec['Convolution']['cnn_weight_fill'] = [{'type' : 'xavier'}]

    layerspec['general'] = {}
    layerspec['general']['output'] = outname + '.model'
    make_model(layerspec)

    #make job script
    learning_rate = np.random.choice(learning_rate_)
    momentum = np.random.choice(momentum_)
    with open(outname + '.pbs', 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#PBS -N %s\n' %outname)
        f.write('#PBS -j oe\n')
        f.write('#PBS -l nodes=1:ppn=1:gpus=1:exclusive_process\n')
        f.write('#PBS -l walltime=28:00:00:00\n')
        f.write('#PBS -q %s\n' %args.queue)
        f.write('\n')
        f.write('echo Running on `hostname`\n')
        f.write('echo workdir $PBS_O_WORKDIR\n')
        f.write('cd $PBS_O_WORKDIR\n')
        f.write('\n')
        f.write('SCRDIR=/scr/jsunseri/$PBS_JOBID\n')
        f.write('if [[ ! -e $SCRDIR ]]; then\n')
        f.write('	mkdir -p $SCRDIR\n')
        f.write('fi\n')
        f.write('\n')
        f.write('echo scratch drive ${SCRDIR}\n')
        f.write('cp %s.model ${SCRDIR}\n' %outname)
        f.write('cd ${SCRDIR}\n')
        f.write('\n')
        #ick, why isn't there a "current cuda version" symlink? I should make one for myself
        f.write('export LD_LIBRARY_PATH=/usr/local/hdf5/lib:/usr/local/cuda-9.0/lib64:/net/pulsar/home/koes/jss97/git/gnina/build/caffe/lib:/net/pulsar/home/koes/dkoes/local/lib:$LD_LIBRARY_PATH\n')
        f.write('export PATH=/net/pulsar/home/koes/jss97/git/gninascripts:$PATH\n')
        f.write('export PYTHONPATH=/net/pulsar/home/koes/jss97/git/gnina/caffe/python:$PYTHONPATH\n')
        f.write('\n')
        f.write('trap "cp * $PBS_O_WORKDIR" EXIT\n')
        f.write('train.py -m %s.model -p %s --base_lr %s --momentum %s -o %s %s &> %s.out\n' 
                %(outname, args.prefix, learning_rate, momentum, 
                    outname, ' '.join(args.trainargs), outname))
        f.write('exit')
