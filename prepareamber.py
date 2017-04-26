#!/usr/bin/python

import argparse
import parsepdb as pdb
import os, shutil, glob
from plumbum.cmd import sed, grep, cut, uniq, wc
try:
    from plumbum.cmd import obabel
except ImportError:
    raise ImportError('Check that obabel is on your path')
try:
    from plumbum.cmd import antechamber, parmchk, pdb4amber, tleap, pmemd_cuda
except ImportError:
    raise ImportError('Check that AMBER binaries are on your path')

def get_cmd(input_str):
    cmd_dict = {'.frcmod' : 'loadamberparams', '.lib' : 'loadoff', '.off' : 
            'loadoff', '.prep' : 'loadAmberPrep', '.mol2' : 'loadmol2', '.pdb' :
            'loadpdb'}
    ext = os.path.splitext(input_str)[-1]
    return cmd_dict[ext]

def get_base(fname):
    '''
    Strip path and extension from filename in order to generate other filenames
    '''
    return os.path.splitext(os.path.basename(fname))[0]

def get_fname(fname):
    '''
    Generate a filename that doesn't overwrite anything in the current
    directory
    '''
    i = 0
    base,ext = fname.split('_')
    if ext: ext = '_' + ext
    while 1:
        if i == 0:
            fname = base + ext
        else:
            fname = base + str(i) + ext
        if os.path.isfile(fname):
            i += 1
        else:
            return fname

def get_molname(molname):
    '''
    Generate a unique molname, important if there are multiple mols
    '''
    molname = molname[:3].upper()
    while os.path.isfile(molname + '.lib'):
        num = ''.join(m for m in molname if m.isdigit())
        num = int(num)+1 if num else 1
        numstr = str(num)
        molname = ''.join(m for m in molname if m.isalpha())
        molname = molname[:3-len(numstr)] + numstr
    return molname

def get_units(lib):
    '''
    Get the names of the units defined by an AMBER library file
    '''
    units = []
    ext = os.path.splitext(lib)[-1]
    with open(lib,'r') as f:
        if ext != 'prep':
            copy = False
            for line in f:
                if line.startswith('!!index array str'):
                    copy = True
                elif line.startswith('!'):
                    break
                elif copy:
                    units.append(line.strip('"'))
        else:
            i = 0
            for line in f:
                if i == 4:
                    units.append(line.split()[0])
                    break
                i += 1
    return units

def get_charge(mol2):
    '''
    Get the net charge on a molecule from a structure file. Currently only
    works with mol2 format
    '''
    with open(mol2, 'r') as f:
        copy = False
        charge = 0
        for line in f:
            if line.startswith('@'):
                record = line.split('>')[-1]
                if record == 'ATOM':
                    copy = True 
                elif record != 'ATOM':
                    if copy == True:
                        break
                    else:
                        continue
            elif (copy == True):
                charge += float(line.split()[-1])
    return charge 

def make_amber_parm(fname, base, ff, wat_dist = 0, libs=[]):
    '''
    Generate AMBER parameters with tleap
    '''
    inpcrd = base + '.inpcrd'
    prmtop = base = '.prmtop'
    with open(base + '.tleap', 'w') as leap_input:
        leap_input.write('source ' + ff + '\n' + 
                'source leaprc.gaff\n' + 
                base + '=' + get_cmd(fname) + ' ' + fname + '\n')
        for lib in libs:
            leap_input.write(get_cmd(lib) + ' ' + lib + '\n')
        if wat_dist:
            leap_input.write('loadamberparams frcmod.ionsjc_tip3p\n' + 
                    'solvateoct prot TIP3PBOX ' + wat_dist + '\n' + 
                    'addions prot Na+ 0\n' + 
                    'addions prot Cl- 0\n')
        leap_input.write('saveamberparm ' + base + ' ' + prmtop + ' ' + inpcrd + '\n')
        leap_input.write('quit\n')
    leap['-f', base + '.tleap']()

def do_amber_min_constraint(base):
    '''
    Do AMBER minimization with protein constraint
    '''
    numres = (grep['ATOM', base + '.pdb'] | cut['-b', '23-26'] | uniq | wc)()
    with open(base + '_min1.in', 'w') as min_input:
        min_input.write(base + ': initial minimization solvent + ions\n' + 
            ' &cntrl\n' + 
            '  imin   = 1,\n' + 
            '  maxcyc = 1000,\n' +
            '  ncyc   = 500,\n' +
            '  ntb    = 1,\n' +
            '  ntr    = 1,\n' +
            '  cut    = 10.0\n' +
            ' /\n' +
            'Hold the protein fixed\n' +
            '500.0\n' +
            'FIND\n' +
            '* * S *\n' +
            '* * B *\n' +
            '* * 3 *\n' +
            '* * E *\n' +
            '* * M *\n' +
            'SEARCH\n' +
            'RES 1 ' + str(numres) + '\n' +
            'END\n' +
            'END\n')
    pmemd_cuda['-O', '-i', base+'_min1.in', '-o', base+'_min1.out', '-p',
            base+'.prmtop', '-c', base+'.inpcrd', '-r', base+'_min1.rst',
            '-ref', base+'.inpcrd']()

def do_amber_min(base):
    '''
    Do unconstrained AMBER minimization 
    '''
    #TODO: assumes order and naming convention
    with open(base + '_min2.in', 'w') as min_input:
        min_input.write(base + ':  minimization whole system\n' + 
            ' &cntrl\n' + 
            '  imin   = 1,\n' + 
            '  maxcyc = 5000,\n' + 
            '  ncyc   = 2000,\n' + 
            '  ntb    = 1,\n' + 
            '  ntr    = 0,\n' + 
            '  cut    = 10.0\n' + 
            ' /\n')
    pmemd_cuda['-O', '-i', base+'_min2.in', '-o', base+'_min2.out', '-p',
            base+'.prmtop', '-c', base+'_min1.rst', '-r', base+'_min2.rst']()

def do_amber_warmup(base, temperature):
    '''
    Do AMBER MD to gradually increase system to target temp
    '''
    #TODO: assumes order and naming convention
    numres = (grep['ATOM', base + '.pdb'] | cut['-b', '23-26'] | uniq | wc)()
    with open(base + '_md1.in','w') as md_input:
        md_input.write(' &cntrl\n' + 
              '  imin   = 0,\n' + 
              '  irest  = 0,\n' + 
              '  ntx    = 1,\n' + 
              '  ntb    = 1,\n' + 
              '  cut    = 10.0,\n' + 
              '  ntr    = 1,\n' + 
              '  ntc    = 2,\n' + 
              '  ntf    = 2,\n' + 
              '  tempi  = 0.0,\n' + 
              '  temp0  = $temp,\n' + 
              '  ntt    = 3,\n' + 
              '  gamma_ln = 1.0,\n' + 
              '  nstlim = 50000, dt = 0.002, ntxo = 2,\n' + 
              '  ntpr = 1000, ntwx = 1000, ntwr = 10000,\n' + 
              '  ioutfm = 1\n' + 
              ' /\n' + 
            'Keep prot fixed with weak restraints\n' + 
            '10.0\n' + 
            'FIND\n' + 
            '* * S *\n' + 
            '* * B *\n' + 
            '* * 3 *\n'  
            '* * E *\n' + 
            '* * M *\n' + 
            'SEARCH\n' + 
            'RES 1 ' + numres + '\n' + 
            'END\n' + 
            'END\n')
    pmemd_cuda['-O', '-i', base+'_md1.in', '-o', base+'_md1.out', '-p',
            base+'.prmtop', '-c', base+'_min2.rst', '-r', base+'_md1.rst',
            '-ref', base+'_min2.rst', '-x', base+'_md1.nc']()

def do_amber_constant_pressure(base):
    '''
    Do AMBER MD to equilibrate system at constant pressure
    '''
    #TODO: assumes order and naming convention
    with open(base + '_md2.in','w') as md_input:
        md_input.write(base + ': 100ps MD\n' + 
             ' &cntrl\n' + 
             '  imin = 0, irest = 1, ntx = 7,\n' + 
             '  ntb = 2, pres0 = 1.0, ntp = 1,\n' + 
             '  taup = 2.0,\n' + 
             '  cut = 10.0, ntr = 0,\n' + 
             '  ntc = 2, ntf = 2,\n' + 
             '  tempi = $temp, temp0 = $temp,\n' + 
             '  ntt = 3, gamma_ln = 1.0,\n' + 
             '  nstlim = 50000, dt = 0.002, ntxo = 2,\n' + 
             '  ntpr = 5000, ntwx = 5000, ntwr = 500000,\n' + 
             '  ioutfm = 1\n' + 
             ' /\n')
    pmemd_cuda['-O', '-i', base+'_md2.in', '-o', base+'_md2.out', '-p',
            base+'.prmtop', '-c', base+'_md1.rst', '-r', base+'_md2.rst',
            '-x', base+'_md2.nc']()

def make_amber_production_input(base, length, keep_velocities):
    '''
    Make input files for production run AMBER MD; length is in nanoseconds
    '''
    nstlim = float(length) / .000002
    irest = int(keep_velocities)
    ntx = 7 if keep_velocities else 1
    with open(base + '_md3.in','w') as md_input:
        md_input.write(base + ': ' + length + 'ns MD\n' + 
            ' &cntrl\n' + 
            '  imin = 0, irest = '+irest+', ntx = '+ntx+',\n' + 
            '  ntb = 2, pres0 = 1.0, ntp = 1,\n' + 
            '  taup = 2.0,\n' + 
            '  cut = 10.0, ntr = 0,\n' + 
            '  ntc = 2, ntf = 2,\n' + 
            '  tempi = $temp, temp0 = $temp,\n' + 
            '  ntt = 3, gamma_ln = 1.0,\n' + 
            '  nstlim = '+nstlim+', dt = 0.002, ntxo = 2,\n' + 
            '  ntpr = 5000, ntwx = 5000, ntwr = 500000,\n' + 
            '  ioutfm = 1\n' + 
             '/\n')

def do_amber_preproduction(base, args):
    '''
    Do minimization with constraints, minimization without constraints, initial
    MD as temperature is raised to target temp, second MD where system is
    equilibrated at constant pressure, and generate input files for production
    run MD but don't run it (becuz it's PREproduction, see?)
    '''
    do_amber_min_constraint(base)
    do_amber_min(base)
    do_amber_warmup(base, args.temperature)
    do_amber_constant_pressure(base)
    make_amber_production_input(base, args.prod_length, args.keep_velocities)

def do_amber_production(base):
    '''
    Does AMBER production run MD locally
    '''
    pmemd_cuda['-O', '-i', base+'_md3.in', '-o', base+'_md3.out', '-p',
            base+'.prmtop', '-c', base+'_md2.rst', '-r', base+'_md3.rst',
            '-x', base+'_md3.nc']()

def do_antechamber(fname, net_charge, ff, molbase = ''):
    '''
    Run antechamber and get correctly named versions of the following: mol2
    with bcc charges, frcmod, lib, prmtop, inpcrd
    '''
    if not molbase: molbase = get_base(fname)
    ext = os.path.splitext(fname)[-1]
    mol2 = molbase + '.mol2'
    antechamber['-i', fname, '-fi', ext, '-o', mol2, '-fo', 'mol2', '-c',
            'bcc', '-nc', str(net_charge), '-s', '2']()
    sed['-i',"'s/\<MOL\>/%s/g' %s" % (molbase, mol2)]()
    frcmod = molbase + '.frcmod'
    parmchk['-i', mol2, '-f', 'mol2', '-o', frcmod]()
    make_amber_parm(mol2, molbase, ff)

def set_matches(fname, libs, reslist, orphaned_res):
    '''
    Find whether any units defined by a lib are required; if they are, update
    the liblist to include that lib and remove the units it defines from
    orphaned_res
    '''
    units = get_units(fname)
    matches = set(units).intersection(reslist)
    if matches:
        #avoid redefining anything, but warn the user about the
        #duplication...this doesn't avoid _all_ redefinitions,
        #though, FIXME?
        if not matches.intersection(orphaned_res):
            print 'Unit %s found in %s defined previously, \
            skipping to avoid redefinition\n' % (' '.join(match for match in
                matches), fname)
        else:
            libs.add(fname)
            frcmod = get_base(fname) + '.frcmod'
            if os.path.isfile(frcmod):
                lib.add(frcmod)
            orphaned_res -= matches

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generates pre-production files for AMBER \
    MD. Can handle a receptor or ligand by themselves (checks for ligand library \
    files in the current directory and generates them if they don't exist) or \
    sets up the complex if given both a receptor and ligand.")

    parser.add_argument('-s', '--structures', nargs='+', required=True, help='Structures \
    for which you want to run a simulation. N.B. if more than one is provided \
    they will be simulated together.')

    parser.add_argument('-n', '--complex_name', default='complex',
    help='Optionally provide a filename for the created complex; only \
    meaningful if you provide more than one structure. Default is "complex."')

    parser.add_argument('-p', '--prefix', nargs='+', required=False, help="Optionally specify \
    a prefix for nonstandard residue library files; this can include their path if they aren't \
    in the current directory. If the relevant files exist we assume you want \
    to use them, otherwise we assume this is where you want them to go and \
    derive the residue name accordingly.")

    parser.add_argument('-w', '--water_dist', default=12, help='Water box \
    distance; defaults to 12.')

    parser.add_argument('-ff', '--force_field', default='oldff/leaprc.ff14SB',
            help='Force field; defaults to ff14SB.')

    parser.add_argument('-t', '--temperature', default=300, help='Simulation \
    temperature; defaults to 300K.')

    parser.add_argument('-l', '--prod_length', default=100, help='Length of \
    production-run MD. This is used to generate the input files for that MD, \
    but note that by default this MD is not run. Units are nanoseconds.')

    parser.add_argument('-k', '--keep_velocities', default=False,
            action='store_true', help='Keep velocities from preproduction run \
                    when starting the production MD. Default is false.')

    parser.add_argument('-r', '--run_prod_md', default=False,
    action='store_true', help='Run production MD locally when all preprocessing \
    is finished. Default is false, because you might want to run it on a \
    cluster.')

    parser.add_argument('-noh', '--no_touch_hyd', dest='noh', default=False,
            action='store_true', help="Don't remove any hydrogens.")

    parser.add_argument('-wat', '--water_cutoff', required=False, help="If \
    specified, keep waters that are the provided distance from nonprotein \
    molecules; if -1, keep all waters in the input")

    parser.add_argument('-parm', '--parm_only', action='store_true', default =
    False, help="Only generate the necessary ligand parameters, don't do the \
    preproduction MDs")

    args = parser.parse_args()

    #do we have nonstandard residues?
    mol_data = {}
    nonstandard_res = {}
    for structure in args.structures:
        assert os.path.isfile(structure),'%s does not exist\n' % structure
        mol_res = {}
        mol_data[structure] = pdb.parsepdb(structure)
        mol_res[structure] = set(mol_data[structure].resname)
        nonstandard_res[structure] = list(mol_res[structure] - standard_res)

    #if so, do we have the necessary library files? 
    #check for prep, lib, and off; just add the frcmod if there is one
    libs = set([])
    for struct,reslist in nonstandard_res.items():
        #track which units you don't have libs for
        orphaned_res = set(reslist)
        #try any user-provided locations first
        for prefix in args.prefix:
            for ext in ['.lib','.off','.prep']:
                fname = prefix + ext
                if os.path.isfile(fname):
                    set_matches(fname, libs, reslist, orphaned_res)
        #if necessary, check the current directory too
        if orphaned_res:
            local_libs = [name for name in glob.glob('*.lib') +
                    glob.glob('*.off') + glob.glob('*.prep')]
            for lib in local_libs:
                set_matches(fname, libs, reslist, orphaned_res)
   
        is_protein = pdb.is_protein(mol_data[struct])
        #for now, require that the ligand be provided separately from the protein -
        #that way we don't need to worry about differentiating between modified
        #residues (or other things we don't want to strip out of the protein) and
        #small molecules we can parametrize with antechamber
        assert(not (is_protein and orphaned_res),
        "Undefined units in protein - check for modified residues, ions, or \
        cofactors\n")

        if is_protein and not args.noh: mol_data[struct] = pdb.strip_hydrogen(mol_data[struct])

        #TODO? we can actually do this for the user
        assert(len(orphaned_res)==1, "%s has multiple ligands; break them into \
        separate files to process with antechamber\n" % struct)

        molname = get_molname(args.prefix[:3].upper() if args.prefix else
                'LIG')
        #if we're handling a ligand and don't have library files, we will need at 
        #least the pdb-formatted data and a mol2 from which we can derive gasteiger 
        #charges for antechamber; make these with babel and find the net charge
        if orphaned_res:
            mol_data[struct] = rename_atoms(mol_data[struct])
            tempname = get_fname(get_base(struct) + '_temp.pdb')
            ligname = get_fname(get_base(struct) + '_amber.pdb')
            mol2 = os.path.splitext(ligname)[0] + '.mol2'
            if not has_hydrogen(mol_data[struct]):
                pdb.writepdb(mol_data[struct], tempname)
                obabel[tempname, '-O', ligname, '-h']()
                os.remove(tempname)
            else:
                pdb.writepdb(mol_data[struct], ligname)
            obabel[ligname, '-O', mol2]()
            net_charge = get_charge(mol2)
            do_antechamber(ligname, net_charge, args.force_field, molname)
            libs.add(molname + '.lib')
            libs.add(molname + '.frcmod')

    #ok, now we can be pretty sure we know what to do and that we are able to do it
    complex = args.complex_name + '.pdb'
    pdb.writepdb(mol_data.values(), complex)
    base = os.path.splitext(complex)[-1]
    make_amber_parm(complex, base, args.ff, args.water_dist, libs)
    if not args.parm_only: do_amber_preproduction(complex, args)
    if args.run_prod_md: do_amber_production(base)
