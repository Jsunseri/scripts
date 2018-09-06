#!/usr/bin/python
from copy import deepcopy
import traceback

class layer:
    '''
    Generates caffe layers to automate model generation.
    
    Attributes:
        name: layer identifier, specific to individual instantiation
        type: kind of layer, common to all instantiations of this layer type
        top: list of names of output blobs
        bottom: list of names of input blobs
        phase: optional inclusion phase
        params: dict of layer-specific param dicts
    '''
    layer_count = {}

    def __init__(self, type = '', params={}, top=[], bottom=[], name = '', phase=''):
        self.name = name
        self.type = type
        if type == 'MolGridData':
            #reset count when we start a new model
            for key in self.layer_count.keys():
                self.layer_count.pop(key)
        if self.type not in self.layer_count:
            self.layer_count[self.type] = 0
        self.layer_count[self.type] += 1
        if not self.name:
            self.name = self.type.lower() + str(self.layer_count[self.type])
        self.top = top
        if not self.top:
            self.top = [self.name]
        self.bottom = bottom
        self.phase = phase
        self.params = deepcopy(params)

    def add_top(self, top):
        self.top.append(top)

    def get_top(self):
        return deepcopy(self.top)

    def add_bottom(self, bottom):
        self.bottom.append(bottom)

    def get_bottom(self):
        return deepcopy(self.bottom)

    def add_param(self, param):
        self.params.update(param)

    def add_phase(self, phase):
        self.phase.append(phase)

    def get_type(self):
        return self.type

    def format_paramname(self, name):
        #TODO: remove? hacky...
        if name == 'MolGridData':
            return 'molgrid_data_param'
        elif name == 'AffinityLoss':
            return 'affinity_loss_param'
        elif name == 'InnerProduct':
            return 'inner_product_param'
        elif name == 'LSTM':
            return 'recurrent_param'
        elif name == 'weight_filler' or name == 'shape':
            return name
        else:
            return name.lower() + '_param'

    def write(self, ofile):
        '''
        Write self in protobuf format to ofile. 
        '''
        ofile.write('layer {\n')
        ofile.write('  name: "%s"\n' %self.name)
        ofile.write('  type: "%s"\n' %self.type)
        for b in self.bottom:
            ofile.write('  bottom: "%s"\n' %b)
        for t in self.top:
            ofile.write('  top: "%s"\n' %t)
        if self.phase:
            ofile.write('  include {\n')
            ofile.write('    phase: %s\n' %self.phase)
            ofile.write('  }\n')
        self.write_params(self.params, ofile)
        ofile.write('}\n')

    def write_params(self, params, ofile, nspaces=2):
        '''
        Special formatting function for writing layer params. 
        params is a dict mapping strings to values that may include other dicts
        '''
        for param,pargs in params.iteritems():
            ofile.write('{}{} {{\n'.format(nspaces * ' ', self.format_paramname(param)))
            for paramname,value in pargs.iteritems():
                if isinstance(value, dict):
                    self.write_params({paramname: value}, ofile, nspaces+2)
                else:
                    if paramname == 'dim':
                        for dim in value:
                            ofile.write('{}{}: {}\n'.format((nspaces +2)* ' ', 
                                paramname, dim))
                    elif paramname == 'pool':
                        ofile.write('{}{}: {}\n'.format((nspaces +2)* ' ', 
                            paramname, value))
                    else:
                        if isinstance(value, str):
                            if paramname == "variance_norm":
                                ofile.write('{}{}: {}\n'.format((nspaces +2)* ' ', 
                                    paramname, value))
                            else:
                                ofile.write('{}{}: "{}"\n'.format((nspaces +2)* ' ', 
                                    paramname, value))
                        elif isinstance(value, bool):
                            ofile.write('{}{}: {}\n'.format((nspaces +2)* ' ', 
                                paramname, str(value).lower()))
                        else:
                            ofile.write('{}{}: {}\n'.format((nspaces +2)* ' ', 
                                paramname, value))
            ofile.write('{}}}\n'.format(nspaces * ' '))
