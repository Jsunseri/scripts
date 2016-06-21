#!/usr/bin/python

# Generate conformers for Paramfit input. Infers parameters to be sampled using
# an input frcmod file

import sys, os, argparse, random
import subprocess as sp

def read_from_frcmod(frcmod, params, all_params, exclude, mask):
    '''
    Set up lists of bonded parameters to be fit based on frcmod
    '''
    with open(frcmod) as file:
        for line in file:
            contents = line.strip()
            if contents:
                maybe_param = filter(None, contents.split("  "))[0].lower()
                if maybe_param in all_params:
                    index = all_params.index(maybe_param)
                    mask = [0 for item in mask]
                    mask[index] = 1
                elif maybe_param in exclude:
                    mask = [0 for item in mask]
                    continue
                else:
                    if sum(mask):
                        idx = [i for i,e in enumerate(mask) if e][0]
                        params[all_params[idx]].append(maybe_param)
    count=0
    for param in all_params:
        count += len(params[param])
    return count

def all_params_frcmod(params,all_params,exclude,mask):
    '''
    Generates an frcmod with all parameters and reads that in
    '''
    if not os.path.isfile(args.prefix+'.mol2'):
        print 'Error: prefix does not correspond to a mol2.'
        sys.exit()
    frcmod = 'all_params.frcmod'
    cmd = ['parmchk','-i',args.prefix+'.mol2','-f','mol2','-o',frcmod,'-a','Y']
    p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    out,err = p.communicate()
    if err:
        # If we can't generate parameters, there's nothing else we can do
        print err
        sys.exit()
    count = read_from_frcmod(frcmod,params,all_params,exclude,mask)
    if count == 0:
        print 'Failed to generate parameters from mol2 with provided prefix. Check system and try again.'
        sys.exit()

def make_directory(base):
    i = 0
    while 1:
        try:
            if i == 0:
                dirname = base
            else:
                dirname = base + '_' + str(i)
            os.mkdir(dirname)
            break
        except OSError:
            if not os.path.isdir(dirname):
                raise
            i += 1
            pass 
    return dirname

# TODO: maybe make it possible to specify a particular bond for
# parametrization? or maybe not...seems like it might be a bad idea to do that
parser = argparse.ArgumentParser(description='Generate Paramfit inputs')
parser.add_argument('--frcmod', default='LIG.frcmod')
parser.add_argument('--lib', default='LIG.lib')
parser.add_argument('--prefix', default='LIG')
parser.add_argument('--fit_all',default='N')
parser.add_argument('-e', '--exhaustiveness', default=2)
args = parser.parse_args()

# Create directory for the starting conformers to be stored in; supports
# multiple such directories existing with distinct names
dirname = make_directory('starting_conformers')

params = {}
# exclude impropers entirely because they aren't bonded?
all_params = ['bond','angle','dihe']
param_range = [(0.5,3.5),(10,180),(0,180)]
param_vars = [2,3,4]
exclude = ['improper','mass','nonbon']
mask = [0] * len(all_params)
for param in all_params:
    params[param] = []

if args.fit_all=='N':
    frcmod = args.frcmod
    count = read_from_frcmod(frcmod,params,all_params,exclude,mask)
    # If frcmod doesn't contain any bonded parameters to be fit, generate a
    # temporary frcmod if possible and warn user that all parameters will be fit
    if count == 0:
        print "frcmod contains no parameters to be fit. Fitting all parameters."
        all_params_frcmod(all_params,exclude,mask)
else:
    all_params_frcmod(all_params,exclude,mask)

# Use lib to generate some relevant data structures: list of atom types, in the
# "atom number" order; dictionary of bonds that are double bonds or more, to
# fix the Gaussian input file later since paramfit assumes everything is a
# single bond; a list of adjacency lists
atoms_flag = 0
connect_flag = 0
# Array of atoms consists of tuples of (name,type) pairs
atoms = []
higher_bonds = {}

with open(args.lib) as file:
    for line in file:
        contents = line.strip()
        if contents.startswith('!'):
            if contents.split(' ')[0].split('.')[-1] == 'atoms' and not atoms_flag:
                atoms_flag += 1
            elif contents.split(' ')[0].split('.')[-1] == 'connectivity' and not connect_flag:
                connect_flag += 1
                adjacency_list = [[] for i in range(len(atoms))]
            elif atoms_flag == 1:
                atoms_flag += 1
            elif connect_flag == 1:
                connect_flag += 1
        else:
            if atoms_flag == 1:
                atoms.append((contents.split(' ')[0].strip('"'),contents.split(' ')[1].strip('"')))
            elif connect_flag == 1:
                fields = contents.split(' ')
                fields = [int(i) for i in fields]
                adjacency_list[fields[0]-1].append(fields[1]-1)
                adjacency_list[fields[1]-1].append(fields[0]-1)
                if fields[2] > 1:
                    if fields[2] not in higher_bonds:
                        higher_bonds[fields[2]] = []
                    higher_bonds[fields[2]].append((fields[0],fields[1]))

# Generate all atom subsequences of 4, 3, or 2 atoms using the adjacency lists
# TODO: probably should omit generating length and angle parameters for bonds
# inside of rings...currently hoping that doing this will result in high 
# energies that lead to structures being tossed out, but we'll see.
# The following is a bit ugly, but I think it's mostly due to using string
# formatting to pad out individual characters to the correct number of spaces
# to match the format of the frcmod files.
subsequences = {}
for param in all_params:
    subsequences[param] = []

for index,sublist in enumerate(adjacency_list):
    for child in sublist:
        forward = '{0:<2}'.format(atoms[index][1]) + '-' + atoms[child][1]
        reverse = (atoms[child][0],atoms[index][0])

        if forward in params['bond'] and reverse not in subsequences['bond']:
            subsequences['bond'].append((atoms[index][0], atoms[child][0]))

        for gchild in adjacency_list[child]:
            forward = ('{0:<2}'.format(atoms[index][1]) + '-' +
                    '{0:<2}'.format(atoms[child][1]) + '-' + atoms[gchild][1])
            reverse = (atoms[gchild][0],atoms[child][0],atoms[index][0])

            if (forward in params['angle'] and reverse not in subsequences['angle'] and atoms[index][0] != atoms[gchild][0]):
                subsequences['angle'].append((atoms[index][0],atoms[child][0],atoms[gchild][0]))

            for ggchild in adjacency_list[gchild]:
                forward = ('{0:<2}'.format(atoms[index][1]) + '-' +
                        '{0:<2}'.format(atoms[child][1]) + '-' +
                        '{0:<2}'.format(atoms[gchild][1]) + '-' +
                        atoms[ggchild][1])
                reverse = (atoms[ggchild][0],atoms[gchild][0],atoms[child][0],atoms[index][0])

                if (forward in params['dihe'] and reverse not in
                        subsequences['dihe'] and atoms[index][0] !=
                        atoms[gchild][0] and atoms[child][0] != atoms[ggchild][0]):
                    subsequences['dihe'].append((atoms[index][0],atoms[child][0],atoms[gchild][0],atoms[ggchild][0]))

# TODO: maybe make alternative exhaustiveness values correspond to different
# sampling options
if args.exhaustiveness != 2:
    print 'Exhaustiveness values other than 2 not currently supported. Continuing with exhaustiveness of 2.'

# Store all the values generated in values.dat for reference
vfile = open('values.dat','w')
gentraj = open('gen_mdcrd.traj','w')

idx = 0
print 'Generating conformers that sample parameters of interest.'
for i,param in enumerate(all_params):
    for sequence in subsequences[param]:
        for j in range(0,200):
            # Generate random param value in permitted range
            value = random.uniform(param_range[i][0],param_range[i][1]) 
            
            # Generate leap input file for imposing the new parameter value
            leap = open('param_leap.in','w')
            leap.write('source leaprc.gaff\n')
            leap.write('loadamberparams %s\n' % (frcmod))
            leap.write('loadoff %s\n' % (args.lib))
            leap.write('lig=loadpdb %s.pdb\n' % (args.prefix))
            if param_vars[i] == 4:
                leap.write('impose lig {1}, { {"%s" "%s" "%s" "%s" %3.4f} }\n' %
                    (sequence[0],sequence[1],sequence[2],sequence[3],value))
                vfile.write('(%s,%s,%s,%s) %3.4f\n' %
                        (sequence[0],sequence[1],sequence[2],sequence[3],value))
            elif param_vars[i] == 3:
                leap.write('impose lig {1}, { {"%s" "%s" "%s" %3.4f} }\n' %
                    (sequence[0],sequence[1],sequence[2],value))
                vfile.write('(%s,%s,%s) %3.4f\n' %
                        (sequence[0],sequence[1],sequence[2],value))
            elif param_vars[i] == 2:
                leap.write('impose lig {1}, { {"%s" "%s" %3.4f} }\n' %
                    (sequence[0],sequence[1],value))
                vfile.write('(%s,%s) %3.4f\n' % (sequence[0],sequence[1],value))
            leap.write('saveamberparm lig tmp.prmtop "%s/LIG_struc_%d.mdcrd"\n' %
                    (dirname, idx))
            leap.write('quit')
            leap.close()
    
            # Run leap; if error occurs with a given conformation just keep going
            # if possible (TODO: error handling?)
            cmd = ['tleap', '-f', 'param_leap.in']
            p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE)
            out,err = p.communicate()
            if err:
                print err
            
            gentraj.write('trajin %s/LIG_struc_%d.mdcrd\n' % (dirname, idx))
            idx += 1

vfile.close()

# Collect all trajectories into one file
gentraj.write('trajout LIG_all.mdcrd\n')
gentraj.write('go')
gentraj.close()
p = sp.Popen('cpptraj -p tmp.prmtop < gen_mdcrd.traj',stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
out,err = p.communicate()
if err:
    print err

p = sp.Popen(['rm', '-f', 'gen_mdcrd.traj'],stdout=sp.PIPE,stderr=sp.PIPE)
out,err = p.communicate()
if err:
    print err
# Use paramfit to prune conformers just to those beneath some energy threshold
dummy = open('dummy.dat','w')
jobfile = open('job.in','w')
valid_dat = open('LIG_valid.dat','w')
gentraj = open('gen_mdcrd.traj','w')
energy_threshold = 2000
print 'Pruning conformers to only those under %3.4f kcal/mol.' % (energy_threshold)

dummy.write('0.000\n')
dummy.close()
jobfile.write('ALGORITHM=NONE\nNSTRUCTURES=1\nCOORDINATE_FORMAT=RESTART\n')
jobfile.close()

good_structs=0
for i in range(0,idx):
    restart = dirname + '/LIG_struc_'+str(i)+'.mdcrd'
    param_command = "paramfit -i job.in -p tmp.prmtop -c " + restart + " -q dummy.dat"
    cmd = param_command + " | grep 'Calculated energy with initial parameters' | awk '{print $10}'"
    p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
    out,err = p.communicate()
    if err:
        print err
    if float(out) < energy_threshold:
        valid_dat.write('LIG_struc_%d.mdcrd\n' % (i))
        gentraj.write('trajin %s/LIG_struc_%d.mdcrd\n' % (dirname,i))
        good_structs += 1

gentraj.write('trajout LIG_valid.mdcrd\ngo')
valid_dat.close()
gentraj.close()
p = sp.Popen('cpptraj -p tmp.prmtop < gen_mdcrd.traj', stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
out,err = p.communicate()
if err:
    print err
sp.call(['rm', '-f', 'job.in'])
sp.call(['rm', '-f', 'dummy.dat'])
sp.call(['rm', '-f', 'gen_mdcrd.traj'])
print "Total good structures=%d." % (good_structs)

# Now generate quantum input files using paramfit
header = open('gaussian.header','w')
header.write('%RWF=job.rwf\n%Chk=job.chk\n%NoSave\n%Mem=2GB\n#SP b3lyp/6-31G* geom=connectivity\n')
# NProcShared=4\n')
header.close()

q_in = make_directory('quantum_inputs')

qjob = open('job_quantum.in','w')
qjob.write('NSTRUCTURES=' + str(good_structs) + '\n')
qjob.write('COORDINATE_FORMAT=TRAJECTORY\n')
qjob.write('RUNTYPE=CREATE_INPUT\n')
qjob.write('QMFILEFORMAT=GAUSSIAN\n')
qjob.write('QMHEADER=gaussian.header\n')
qjob.write('QMFILEOUTSTART=' + q_in + '/quantum_LIG_\n')
qjob.write('QMFILEOUTEND=.com\n')
qjob.close()

cmd = 'paramfit -i job_quantum.in -p tmp.prmtop -c LIG_valid.mdcrd > paramfit_quantum.out'
p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
out,err = p.communicate()
if err:
    print err

# Edit quantum input files so that bond order information is correct.

# I stored the bond order information above so I could do this file parsing
# manually. Currently I am generating a temporary gaussian input from an sdf in
# the directory, then copying the connectivity section and pasting it onto the
# bottom of the inputs. Probably this is a bad idea - requires the file to
# exist or I have to create it, assumes various files in the directory and
# might be tricked by files with the same name but different atom orderings,
# etc. Also, investigate whether Gaussian has an order value for aromatics.

p = sp.Popen(['obabel',args.prefix+'.sdf','-O','gaussian.com','-xb'],stdout=sp.PIPE,stderr=sp.PIPE)
out,err = p.communicate()
if err:
    print err

space_count = 0
connect_list = []
with open('gaussian.com','r') as file:
    for line in file:
        if not line.strip():
            space_count += 1    
        if space_count == 3:
            connect_list.append(line)

os.system('rm -f gaussian.com')

for file in os.listdir(q_in):
    space_count = 0
    with open('tmp.com','w') as new:
        with open(q_in + '/' + file,'r') as old:
            for line in old:
                if not line.strip():
                    space_count += 1
                if space_count < 3:
                    new.write(line)
        for line in connect_list:
            new.write(line)
        new.write('\n')
    cmd = 'mv tmp.com ' + q_in + '/' + file
    os.system(cmd)

print 'Generating PBS job control files in quantum output directory.'
q_out = make_directory('quantum_outputs')

# Generate PBS job file and exit
for input in os.listdir(q_in):
    out_number = os.path.splitext(input)[0].split('_')[-1]
    with open(q_out + '/run' + out_number + '.pbs','w') as file:
        file.write('#!/bin/bash\n')
        file.write('#PBS -N ' + os.getcwd().split('/')[-2] + '_' + out_number +
                '\n')
        file.write('#PBS -j oe\n')
        file.write('#PBS -l nodes=1:ppn=1\n')
        file.write('#PBS -l walltime=28:00:00:00\n')
        file.write('#PBS -q dept_24_core\n')
        file.write('\n')
        file.write('echo Running on `hostname`\n')
        file.write('echo workdir $PBS_O_WORKDIR\n')
        file.write('echo ld_library_path $LD_LIBRARY_PATH\n')
        file.write('\n')
        file.write('cd $PBS_O_WORKDIR\n')
        file.write('\n')
        file.write('SCRDIR=/scr/$PBS_JOBID\n')
        file.write('\n')
        file.write('if [[ ! -e $SCRDIR ]]; then\n')
        file.write('\tmkdir $SCRDIR\n')
        file.write('fi\n')
        file.write('\n')
        file.write('echo scratch drive ${SCRDIR}\n')
        file.write('cp $PBS_O_WORKDIR/../' + q_in + '/' + input + ' ${SCRDIR}\n')
        file.write('cd ${SCRDIR}\n')
        file.write('module load g09B.01\n')
        file.write('\n')
        file.write('trap "mv quantum_LIG_' + out_number + '.out $PBS_O_WORKDIR" EXIT\n')
        file.write('g09 <' + input + ' >quantum_LIG_' + out_number + '.out\n')
        file.write('cp quantum_LIG_' + out_number + '.out $PBS_O_WORKDIR\n')
        file.write('exit\n')
