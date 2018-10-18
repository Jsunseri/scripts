#!/usr/bin/python
# Autogenerate PBS run file for MD
import sys, os

i = sys.argv[1]

# Pass the job name as the second argument
# on the command line, or it will be 
# inferred from the directory structure.
if sys.argv > 2:
    dir_name = sys.argv[2]
else:
    dir_name = os.getcwd()
    dir_name = dir_name.strip('"')
    dir_name = dir_name.split('/')[-3]

run_script = open("run.pbs","w")
run_script.write('#!/bin/bash\n#PBS -N %s\n' % (dir_name))
run_script.write('#PBS -j oe\n')
run_script.write('#PBS -l nodes=1:ppn=1:gpus=1\n')
run_script.write('#PBS -l walltime=28:00:00:00\n')
run_script.write('#PBS -q dept_gpu_8GB\n\n')
run_script.write("i=%s\n\n" % (i))
run_script.write("echo Running on `hostname`\n")
run_script.write("echo workdir $PBS_O_WORKDIR\n")
run_script.write("echo ld_library_path $LD_LIBRARY_PATH\n\n")
run_script.write("cd $PBS_O_WORKDIR\n\n")
run_script.write("SCRDIR=/scr/jsunseri/$PBS_JOBID\n\n")
run_script.write("if [[ ! -e $SCRDIR ]]; then\n")
run_script.write("\tmkdir -p $SCRDIR\n")
run_script.write("fi\n\n")
run_script.write("echo scratch drive ${SCRDIR}\n\n")
run_script.write("cp $PBS_O_WORKDIR/${i}*.in ${SCRDIR}\n")
run_script.write("cp $PBS_O_WORKDIR/${i}*.prmtop ${SCRDIR}\n")
run_script.write("cp $PBS_O_WORKDIR/${i}*_md2.rst ${SCRDIR}\n")
run_script.write("cp $PBS_O_WORKDIR/${i}*.inpcrd ${SCRDIR}\n\n")
run_script.write("cd ${SCRDIR}\n\n")
run_script.write("export AMBERHOME=/usr/local/amber18\n")
run_script.write("export PATH=/usr/local/amber18/bin:$PATH\n")
run_script.write("export \
LD_LIBRARY_PATH=/usr/local/cuda-9.0/lib64:/usr/local/amber18/lib:/usr/include/x86_64-linux-gnu:$LD_LIBRARY_PATH\n\n")
run_script.write('trap "mv *md3.nc $PBS_O_WORKDIR" EXIT\n\n')
run_script.write("pmemd.cuda -O -i ${i}_md3.in -o $PBS_O_WORKDIR/${i}_md3.out \
-p ${i}.prmtop -c ${i}_md2.rst -r ${i}_md3.rst -x ${i}_md3.nc -inf \
$PBS_O_WORKDIR/mdinfo_${i}\n\n")
run_script.write('exit')
run_script.close()
#os.system("find . -type d -exec bash -c "cd '{}' && cocrystal_process.py *.pdb && ligandparams3.py LIG.pdb receptor.pdb" \;")
