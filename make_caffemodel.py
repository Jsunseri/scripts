#!/usr/bin/env python

import argparse
import json
from layer import layer

activation_types = ['ELU', 'ReLU', 'PReLU', 'Sigmoid', 'TanH', 'Power', 'Exp',
        'Log', 'BNLL', 'Threshold', 'Bias', 'Scale']

def add_layer_of_type(type, layerparam, layers):
    '''
    sets up layerparam, top, and bottom for specified layer type
    '''
    params = {}
    #check whether we have specific params for this type
    try:
        params = layerparam[type]
    except:
        pass

    if type == 'MolGridData':
        assert params['affinity'] or params['label'], \
        "Must have at least one of affinity and label."
        name = "data"
        top = ['data']
        aff = params.pop('affinity')
        if aff: 
            top.append('affinity')
            params['has_affinity'] = True
        lab = params.pop('label')
        if lab: 
            top.append('label')
        stratrec = params.pop("stratify_receptor")
        params['root_folder'] = "DATA_ROOT"
        if not params['ligmap']:
            params.pop('ligmap')
        if not params['recmap']:
            params.pop('recmap')
        for phase in ["TEST", "TRAIN"]:
            params['source'] = '%sFILE' %phase
            if phase == "TEST":
                params["shuffle"] = False
                params["balanced"] = False
                params["stratify_receptor"] = False
            else:
                params["shuffle"] = True
                params["balanced"] = True
                params["random_rotation"] = True
                params["random_translate"] = 2
                params["stratify_receptor"] = stratrec
            layers.append(layer(type, layerparam, top, [], name, phase))
        if lab:
            params['label'] = True
        if aff:
            params['predaff'] = True

    if type == 'Convolution' or type == "Pooling" or type == "LSTM":
        bottom = layers[-1].get_top()
        if len(bottom) > 1:
            assert "data" in bottom, "Input layer must create a data blob."
            bottom = ["data"]
        if type == "LSTM":
            bottom.append("seqcont")
        layers.append(layer(type, layerparam, [], bottom))

    if type in activation_types:
        bottom = layers[-1].get_top()
        top = bottom
        layers.append(layer(type, layerparam, top, bottom))

    if type == "InnerProduct":
        name = "output_fc"
        if params['num_output'] == 1:
            name += '_aff'
        i = len(layers)-1
        while i >=0:
            last_type = layers[i].get_type()
            if last_type in activation_types or last_type == "Convolution" or last_type == "LSTM":
                bottom = layers[i].get_top()
                break
            i -= 1
        layers.append(layer(type, layerparam, [], bottom, name))

    if type == "Softmax":
        #must be predicted label
        bottom = ['output_fc']
        name = 'output'
        layers.append(layer(type, layerparam, [], bottom, name))

    if type == 'SoftmaxWithLoss':
        name = 'loss'
        bottom = ['output_fc', 'label']
        top = ['loss']
        layers.append(layer(type, layerparam, top, bottom, name))

    if type == "AffinityLoss":
        name = "rmsd"
        bottom = ['output_fc_aff', 'affinity']
        layers.append(layer(type, layerparam, [], bottom, name))

    if type == "Reshape":
        name = params.pop('name')
        phase = ''
        try:
            phase = params.pop('phase')
        except:
            pass
        bottom = [params.pop('bottom')]
        layers.append(layer(type, layerparam, [], bottom, name, phase))

def make_model(layerspec):
    '''
    layerspec is a dict mapping layer types to param dicts; from those, generate a list
    of layer objects and write a model protobuf
    '''
    layers = []
    pooling_param = {'pool': 'MAX', 'kernel_size': 2, 'stride': 2}
    has_lstm = layerspec['LSTM']['n_lstm_layers'] > 0
    #start with MolGridDataLayer
    #make two layer objects, one each for train and test phase
    layerparam = layerspec['MolGridData']
    add_layer_of_type('MolGridData', {'MolGridData': layerparam}, layers)

    #now check whether we have conv layers and how many
    layerparam = layerspec['Convolution']
    if int(layerparam['n_cnn_layers']) > 0:
        nlayers = int(layerparam.pop('n_cnn_layers'))
        nfilters = layerparam.pop('n_cnn_filters')
        nfilters = [int(elem) for elem in nfilters] #seems unnecessary, figure out why it isn't
        if len(nfilters) < nlayers:
            if len(nfilters) > 1:
                print "Number of CNN filter sizes provided is smaller than \
the number of CNN layers; inferring additional filter dimensions \
from the first.\n"
            for i in range(1, nlayers):
                nfilters.append(int(nfilters[i-1]) * 2)
        stride = layerparam.pop('stride')
        if len(stride) < nlayers:
            stride = stride * nlayers
        pad = layerparam.pop('pad')
        if len(pad) < nlayers:
            pad = pad * nlayers
        kernel_size = layerparam.pop('kernel_size')
        if len(kernel_size) < nlayers:
            kernel_size = kernel_size * nlayers
        activation = layerparam.pop('activation')
        if len(activation) < nlayers:
            if len(activation) > 1:
                print "Warning: fewer activation classes provided than conv layers; \
replicating first specified for all layers.\n"
            activation = activation * nlayers
        cnn_weight_fill = layerparam.pop('cnn_weight_fill')
        if len(cnn_weight_fill) < nlayers:
            if len(cnn_weight_fill) > 1:
                print "Warning: fewer weight fill classes provided than conv layers; \
replicating first specified for all layers.\n"
            cnn_weight_fill = cnn_weight_fill * nlayers

        for i in range(nlayers):
            layerparam['stride'] = stride[i]
            layerparam['pad'] = pad[i]
            layerparam['kernel_size'] = kernel_size[i]
            layerparam['num_output'] = nfilters[i]
            if has_lstm:
                layerparam['axis'] = 2
            layerparam['weight_filler'] = cnn_weight_fill[i]
            add_layer_of_type('Convolution', {"Convolution": layerparam}, layers)
            add_layer_of_type(activation[i], {}, layers)
            add_layer_of_type('Pooling', {"Pooling" : pooling_param}, layers)

    #now check whether we have lstm layers and how many
    layerparam = layerspec['LSTM']
    if int(layerparam['n_lstm_layers']) > 0:
        assert layers[0].get_type() == 'MolGridData'
        layers[0].add_top("seqcont")
        nlayers = int(layerparam.pop('n_lstm_layers'))
        hidden_dim = layerparam.pop('hidden_dim')
        hidden_dim = [int(elem) for elem in hidden_dim] #seems unnecessary, figure out why it isn't
        if len(hidden_dim) < nlayers:
            if len(hidden_dim) > 1:
                print "Number of LSTM hidden dimension sizes provided is less than \
the number of LSTM layers; inferring additional dimensions from the \
first.\n"
            hidden_dim = hidden_dim * nlayers
        layerparam['num_output'] = hidden_dim
        lstm_weight_fill = layerparam.pop('lstm_weight_fill')
        if len(lstm_weight_fill) < nlayers:
            if len(lstm_weight_fill) > 1:
                print "Warning: fewer weight fill classes provided than LSTM layers; \
replicating first specified for all layers.\n"
            lstm_weight_fill = [lstm_weight_fill for lstm_weight_fill in
                    range(nlayers)]
        for i in range(nlayers):
            layerparam['num_output'] = hidden_dim[i]
            layerparam['weight_filler'] = lstm_weight_fill[i]
            add_layer_of_type('LSTM', {"LSTM": layerparam}, layers)

    #now output layer(s)
    outparam = layerspec['MolGridData']
    layerparam = {}
    reshape_param = {}
    loss_param = {}
    for output in [('label', 2), ('predaff', 1)]:
        #inner product, then prediction, loss, and propagated label for TEST phase
        if output[0] not in outparam:
            continue
        layerparam['num_output'] = output[1]
        layerparam['weight_filler'] = {'type': 'xavier'} #don't get a choice...for now
        if has_lstm:
            layerparam['axis'] = 2
        add_layer_of_type('InnerProduct', {"InnerProduct": layerparam}, layers)
        if output[0] == "label":
            layerparam = {}
            if has_lstm:
                layerparam['axis'] = 2
            loss_param["ignore_label"] = -1 #TODO: allow customization?
            add_layer_of_type('Softmax', {"Softmax": layerparam}, layers)
            add_layer_of_type('SoftmaxWithLoss', {"Softmax": layerparam,
                "Loss": loss_param}, layers)
            reshape_param['Reshape'] = {'shape': {'dim': (0,-1)},
                    'name': 'labelout', 'bottom': 'label', 'phase': 'TEST'}
        else:
            reshape_param['Reshape'] = {'shape': {'dim': (0,-1)},
                    'name': 'predaff', 'bottom': 'output_fc_aff'}
            add_layer_of_type('Reshape', reshape_param, layers)
            affinity_loss_param = {'affinity_loss': {'scale': 0.1}}
            add_layer_of_type('AffinityLoss', affinity_loss_param, layers)
            reshape_param['Reshape'] = {'shape': {'dim': (0,-1)},
                    'name': 'affout', 'bottom': 'affout', 'phase': 'TEST'}
        add_layer_of_type('Reshape', reshape_param, layers)

    ofile = open(layerspec['general']['output'], 'w')
    for layer in layers:
        layer.write(ofile)
    ofile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Sets up model prototxt for \
    specific hyperparameter options.")

    parser.add_argument('-o', '--output', default='gnina.model', help='Ouput \
    model name, including .model extension.')
    #MolGridDataLayer options
    molgrid = parser.add_argument_group('MolGridData')

    molgrid.add_argument('-aff', '--affinity', default=False, action =
            'store_true',
            help='Indicate whether input data has affinity.')

    molgrid.add_argument('-lab', '--label', default=False, action= 'store_true', 
            help='Indicate whether input data has label.')

    molgrid.add_argument('-bs', '--batch_size', default=32, help='MolGrid batch \
    size.')

    molgrid.add_argument('-d', '--dimension', default=23.5, help='Total grid \
    dimension.')

    molgrid.add_argument('-r', '--resolution', default=0.5, help='Grid \
    resolution.')

    molgrid.add_argument('-sd', '--subgrid_dim', default=0, help='Subgrid \
    dimension for cubic recurrence. Default is 0 (no recurrence).')
    
    molgrid.add_argument('--ligmap', default='', help='Custom ligmap.')

    molgrid.add_argument('--recmap', default='', help='Custom recmap.')

    molgrid.add_argument('-stratrec', '--stratify_receptor', default=False,
            action='store_true', help='Stratify by receptor when producing \
    examples during training.')

    #CNN layer options
    cnn = parser.add_argument_group('Convolution')

    cnn.add_argument('-cnn', '--n_cnn_layers', default=0, help='Number \
    of input convolutional layers; default=0.')

    cnn.add_argument('-ncf', '--n_cnn_filters', nargs='*', default=[32], help='Number of \
    CNN filters; can be a list for each CNN layer present, or a base \
    value that will be multipled by a factor of 2 for subsequent \
    layers. If none provided, base=32.')

    cnn.add_argument('-cnn_wf', '--cnn_weight_fill', nargs='*',
            default=['xavier'],
            help='Conv layer weight filler type; default=xavier.')

    cnn.add_argument('-kp', '--pad', nargs='*', default=[1], help='Kernel pad for \
    conv layers; if one value provided it is duplicated for all layers, \
    otherwise give a list equal to the number of layers. Default is \
    1.')

    cnn.add_argument('-ksz', '--kernel_size', nargs='*', default=[3], help='Kernel size for \
    conv layers; if one value provided it is duplicated for all layers, \
    otherwise give a list equal to the number of layers. Default is \
    3.')

    cnn.add_argument('-kstr', '--stride', nargs='*', default=[1], help='Kernel stride for \
    conv layers; if one value provided it is duplicated for all layers, \
    otherwise give a list equal to the number of layers. Default is \
    1.')

    cnn.add_argument('-act', '--activation', nargs='*', default=['ELU'],
            help='CNN activation function; default=ELU.')

    #LSTM layer options
    lstm = parser.add_argument_group('LSTM')

    lstm.add_argument('-lstm', '--n_lstm_layers', default=0, help='Number of \
    input LSTM layers; default=0.')

    lstm.add_argument('-hd', '--hidden_dim', nargs='*', default=[1000], help='Hidden unit \
    dimension; default=1000.')

    lstm.add_argument('-lstm_wf', '--lstm_weight_fill', nargs='*',
            default={'type': 'xavier'}, help='LSTM layer weight filler type; \
    default=xavier.', type=json.loads) 

    args = parser.parse_args()
    #parse into layer-specific dicts
    grouped_args = {}
    for group in parser._action_groups[2:]:
        group_options = [action.dest for action in group._group_actions]
        grouped_args[group.title] = {name: value for (name, value) in
                args._get_kwargs() if name in group_options}

    general_group = parser._action_groups[1]
    general_options = [action.dest for action in general_group._group_actions]
    grouped_args['general'] = {name: value for (name, value) in
            args._get_kwargs() if name in general_options}

    make_model(grouped_args)
