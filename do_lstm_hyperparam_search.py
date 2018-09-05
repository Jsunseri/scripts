#!/usr/bin/env python

#generate model files and PBS scripts (using the train.py script) to do a
#random hyperparameter search for the subcube LSTM.
#sample: - LSTM hidden layer dimension
#        - number of input conv layers (and number of filters per layer)
#        - learning rate
#        - momentum
#        - weight initialization
#        - loss (if we're doing affinity prediction)
#        - batch size
#        - subgrid dim
#        - resolution
#        - activation function

from make_caffemodel import make_model
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Generate model prototxt and job \
script for LSTM random hyperparameter optimization.")

parser.add_argument('-n', '--num_models', default=10, help='Number of files to \
generate.')
parser.add_argument('-aff', '--affinity', default=False, action='store_true',
        help='Indicate whether input data has affinity.')
parser.add_argument('-lab', '--label', default=False, action='store_true',
        help='Indicate whether input data has a class label.')
parser.add_argument('-lmap', '--ligmap', default='', help='Location of \
desired ligmap file.')
parser.add_argument('-rmap', '--recmap', default='', help='Location of \
desired recmape file.')
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
learning_rate_ = np.logspace(-6, -1)
momentum_ = 1 - np.logspace(-2, 0)
weight_init_ = ['constant', 'uniform', 'gaussian', 'xavier', 'positive_unitball']
constant_weight = [0.5, 0.75, 1]
gaussian_weight = np.linspace(0, 1)
xavier_weight = ['FAN_IN', 'FAN_OUT', 'AVERAGE']
batch_size_ = [10, 32, 64, 128]
subgrid_dim_ = [6, 8, 10, 12]
resolution_ = [0.25, 0.5, 0.75]
activation_function_ = ['ELU', 'ReLU']

for i in range(args.num_models):
    #make .model file
    layerspec = {}
    layerspec['MolGridData'] = {}
    if args.affinity:
        layerspec['MolGridData']['affinity'] = True
    if args.label:
        layerspec['MolGridData']['label'] = True

    layerspec['MolGridData']['batch_size'] = np.random.choice(batch_size_)
    layerspec['MolGridData']['subgrid_dim'] = np.random.choice(subgrid_dim_)
    layerspec['MolGridData']['resolution'] = np.random.choice(resolution_)
    #at _least_ 23.5, but choose something evenly divisible even though it
    #shouldn't actually be necessary
    dimension = 23.5
    factor = round((dimension - subgrid_dim) / (subgrid_dim + resolution))
    layerspec['MolGridData']['dimension'] = factor * (subgrid_dim + resolution) + subgrid_dim

    layerspec['LSTM'] = {}
    layerspec['LSTM']['hidden_dim']= np.random.choice(hidden_dim_)
    layerspec['LSTM']['n_lstm_layers'] = 1 #let's not bother with varying 
                                           #this yet due to memory issues
    layerspec['LSTM']['weight_filler'] = {}
    weight_init = np.random.choice(weight_init_)
    layerspec['LSTM']['weight_filler']['type'] = weight_init
    #depending on weight init type, we need to make some more choices
    if weight_init == 'constant':
        layerspec['LSTM']['weight_filler']['value'] = np.random.choice(constant_weight)
    if weight_init == 'uniform':
        layerspec['LSTM']['weight_filler']['min'] = 0.25
        layerspec['LSTM']['weight_filler']['max'] = 0.75
    if weight_init == 'gaussian':
        layerspec['LSTM']['weight_filler']['std'] = np.random.choice(gaussian_weight)
    if weight_init == 'xavier':
        layerspec['LSTM']['weight_filler']['variance_norm'] = np.random.choice(xavier_weight)

    layerspec['CNN'] = {}
    layerspec['CNN']['n_cnn_layers'] = np.random.choice(n_cnn_layers_)
    layerspec['CNN']['n_cnn_filters'] = np.random.choice(n_cnn_filters_)

    #make job script
    learning_rate = np.random.choice(learning_rate_)
    momentum = np.random.choice(momentum_)
    activation_function = np.random.choice(activation_function_)
