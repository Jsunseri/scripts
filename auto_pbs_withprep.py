#!/usr/bin/python
#Autogenerate PBS run file for MD
import sys, os

dir_name = os.getcwd()
dir_name = dir_name.strip('"')
values = dir_name.split('/')
dir_name = values[-2]
model_name = values[-1]

run_script = open("run.pbs","w")
run_script.write('#!/bin/bash\n#PBS -N ')
run_script.write(dir_name)
run_script.write('\n#PBS -j oe\n#PBS -l nodes=1:ppn=1:gpus=1\n#PBS -l walltime=24:00:00\n#PBS -q any_gpu\n\ni=com\n\necho Running on `hostname`\necho workdir $PBS_O_WORKDIR\necho ld_library_path $LD_LIBRARY_PATH\n\ncd $PBS_O_WORKDIR\n\nSCRDIR=/scr/$PBS_JOBID\n\nif [[ ! -e $SCRDIR ]]; then\n\tmkdir $SCRDIR\nfi\n\necho scratch drive ${SCRDIR}\n\ncp $PBS_O_WORKDIR/* ${SCRDIR}\n\ncp /net/pulsar/home/koes/jss97/scripts/prepareamber_short.pl ${SCRDIR}\n\ncd ${SCRDIR}\n\nexport AMBERHOME=/usr/local/amber14\nexport PATH=/usr/local/amber14/bin:$PATH\n\ntrap "mv *md3.nc $PBS_O_WORKDIR/com_md3_' + model_name + '.nc" EXIT\n\n./prepareamber_short.pl ${i}.pdb LIG\n\npmemd.cuda -O -i ${i}_md3.in -o $PBS_O_WORKDIR/${i}_md3.out -p ${i}.prmtop -c ${i}_md2.rst -r ${i}_md3.rst -x ${i}_md3.nc -inf $PBS_O_WORKDIR/mdinfo\n\nexit')
run_script.close()
#os.system("find . -type d -exec bash -c "cd '{}' && cocrystal_process.py *.pdb && ligandparams3.py LIG.pdb receptor.pdb" \;")
