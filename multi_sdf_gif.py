#!/usr/bin/python

import argparse
import sys,os,shutil,math
import subprocess as sp

def make_directory(base):
    """ 
    Make directory with error checking that 
    avoids overwriting an existing directory.
    """
    	
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

def extract_data(line):
    """
    Extract useful information from gdb output
    """
    tmp = line.split('=',1)[1]
    tmp = tmp.split('data = ')[1:]
    tmp = [x.split('}') for x in tmp]
    tmp = [x[0] for x in tmp]
    tmp = [x.lstrip('{') for x in tmp]
    return tmp

def roundoff_error(exact, approximate):
    if (exact == 0.0 or approximate == 0.0):
        return abs(exact + approximate)
    else:
        return abs(approximate/exact - 1.0)

def almosteq(float1, float2, epsilon=.0001):
    return (roundoff_error(float1, float2) < epsilon)

#TODO: add option to choose between visualizing energies or forces if gdb
# output is provided? Ultimately maybe want energies for coloring and forces as
# arrows
parser = argparse.ArgumentParser(description='Generates animated GIF from \
multi-state sdf file; if you include correctly formatted GDB output with \
per-atom forces, it can also generate a heatmap of those on the GIF')
parser.add_argument('-l', '--ligand', required=True, help='Ligand sdf to \
generate GIF from; has to have more than one state/molecule if you want \
something animated')
parser.add_argument('-r', '--receptor', required=True, help='Receptor, which \
you might want to include to make the picture prettier')
parser.add_argument('-g', '--gdb_output', required=False, help='Output from \
running gdb such that each ligand molecule has an associated gdb output entry \
that includes both the per-atom forces (to generate heatmap coloring) and the \
atomic coordinates (to provide a mapping between the atom ordering of the sdf \
and the forces array)')
parser.add_argument('-o', '--outname', required=True, help='Output name for \
gif')
args = parser.parse_args()

#see how many ligands (and therefore frames) there are in the input
count = 0
with open(args.ligand,'r') as f:
    for line in f:
        if line.startswith('$$$$'):
            count += 1

#if you have gdb output, generate an array of the magnitudes of the forces for
# each iteration to update the b values in pymol later
if args.gdb_output:
    forces = []
    seen_coords = 0
    copyvar = {}
    copyvar['forces'] = 0
    copyvar['coords'] = 0
    with open(args.gdb_output,'r') as f:
        for line in f:
            if line.startswith('Forces'):
                copyvar['forces'] = 1
            elif line.startswith('Coords'):
                copyvar['coords'] = 1
            elif copyvar['forces'] == 1:
                tforces = extract_data(line)
                # generate magnitude of force for each atom
                for i,elem in enumerate(tforces):
                    force = map(float, elem.split(','))
                    # mag = math.sqrt((float(force[0]))**2 +
                            # (float(force[1]))**2 +
                            # (float(force[2]))**2) 
                    idx,_ = max(enumerate(map(lambda x: x**2, force[:-1])),
                            key=lambda y: y[1])
                    maxf = force[idx]
                    tforces[i] = maxf
                forces.append(tforces)
                copyvar['forces'] = 0
            elif not seen_coords and copyvar['coords'] == 1:
                # reorder forces to match sdf atom ordering
                # first read in gdb coords
                tcoords = extract_data(line)
                for i,elem in enumerate(tcoords):
                    coord = elem.split(',')
                    tcoords[i] = (round(float(coord[0]),4),
                            round(float(coord[1]),4),
                            round(float(coord[2]),4))
                # then parse starting coords from sdf
                line_count = 1
                sdf_coords = []
                # just make natoms bigger than 4 to start, since line 4 is
                # where we will set it to the correct value
                natoms = 100
                with open(args.ligand, 'r') as f:
                    for line in f:
                        if line_count == 4:
                            natoms = int(line.split()[0])
                        elif line_count > 4 and line_count < (natoms + 5):
                            line = line.split()
                            # will the new H-removed code output Hs? might be
                            # able to remove this check
                            if line[3] != 'H':
                                sdf_coords.append((float(line[0]),float(line[1]),float(line[2])))
                        elif line_count > (natoms + 5):
                            break
                        line_count+=1
                if len(tcoords) != len(sdf_coords):
                    print "Error: mismatch between number of coords in gdb output and sdf.\n"
                    sys.exit()
                # now find the indices that will reorder the gdb output to
                # match the sdf
                indices = []
                for j,coord in enumerate(tcoords):
                    match = 0
                    for diff in [.0001,.001,.01]:
                        #relax equality check if you can't find a match
                        for i,elem in enumerate(sdf_coords):
                            if almosteq(coord[0],elem[0],diff) and \
                            almosteq(coord[1],elem[1],diff) and \
                            almosteq(coord[2],elem[2],diff):
                                indices.append(i)
                                match += 1
                        if match > 0:
                            break
                    if match > 1:
                        print "Error: either two atoms in sdf overlap, or atom %d in tcoords has no good match.\n" % j
                        sys.exit()
                seen_coords = 1

if len(forces) != count:
    print "Error: mismatch between number of times forces were output and the number of ligands provided.\n"
    sys.exit()

# reorder the forces to match the sdf atom ordering
for i,force in enumerate(forces):
    forces[i] = [force[j] for j in indices]

# find max and min of forces overall for b-factor coloring 
maxforce = max(map(max,forces))
minforce = min(map(min,forces))

#make a temp directory to store all the pngs
dirname = make_directory('pngs')
shutil.copy(args.ligand,dirname)
shutil.copy(args.receptor,dirname)
os.chdir(dirname)

#generate the pymol input file for making pngs
#N.B. if you're picky about position, you can tune this manually by opening
#pymol, positioning the scene as desired, doing cmd.get_view(), and then adding
#cmd.set_view([$(18 floats from previous step)]) below
pymol_fname = 'pymol_in.py'
ligname = args.ligand.split('.')[0]
with open(pymol_fname,'w') as f:
    f.write('from pymol import cmd\n\n')
    f.write('cmd.load("%s")\n' %args.receptor)
    f.write('cmd.load("%s")\n' %args.ligand)
    f.write('cmd.bg_color("white")\n')
    f.write('cmd.show("sticks")\n')
    f.write('cmd.hide("(hydro)")\n')
    f.write('cmd.show("surface", cmd.get_names()[0])\n')
    f.write('cmd.show("cartoon", cmd.get_names()[0])\n')
    f.write('cmd.set("stick_radius", .15, cmd.get_names()[0])\n')
    f.write('cmd.set("transparency",".2")\n')
    f.write('cmd.util.cbaw(cmd.get_names()[0])\n')
    #example of specializing the view for mhc.pdb and bug.sdf
    f.write('cmd.set_view([0.547671377658844, 0.032215431332588196,\
            -0.8360742926597595, 0.4791826009750366, 0.8070738911628723,\
            0.3449873924255371, 0.6858854293823242, -0.5895708799362183,\
            0.42657366394996643, 0.0, 0.0, -45.422935485839844,\
            -25.35818862915039, 6.457786560058594, 44.009033203125,\
            35.811790466308594, 55.034080505371094, -20.0])\n')
    for i in range(1,count+1):
        f.write('cmd.set("state",%d)\n' %i)
        if args.gdb_output:
            f.write('cmd.select("sele","%s",state=%d)\n' %(ligname,i))
            f.write('cmd.alter("sele", "b=0.0")\n')
            f.write('tlist=[%s]\n' % ','.join(map(str,forces[i-1])))
            f.write('cmd.alter("sele", "b=tlist.pop(0)")\n')
            f.write('cmd.spectrum("b", "green_yellow_red", "sele", %.4f, \
                    %.4f)\n' %(minforce,maxforce))
        f.write('cmd.ray()\n')
        f.write('cmd.png("conf_%d.png")\n' %i)

#run pymol
cmd = 'pymol -cq ' + pymol_fname
p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
out,err = p.communicate()
if err:
        print err

#generate gif with convert
cmd = 'convert -delay 20 -loop 0 $(ls -v *.png) ' + args.outname + '.gif'
p = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True)
out,err = p.communicate()
if err:
        print err

#clean up the detritis
shutil.copy(args.outname + '.gif', '..')
os.chdir('..')
shutil.rmtree(dirname)
