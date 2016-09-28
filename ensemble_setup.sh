#!/bin/bash

for arg
do
    cd "$arg"
    for compound in LIG receptor LIGreceptor
    do
        sed -i -e 's/50000000/2000000/g' -e 's/100ns/4ns/g' ${compound}_md3.in
        sed -i -e 's/ntx = 7/ntx = 1/g' -e 's/irest = 1/irest = 0/g' ${compound}_md3.in
    done
    for ((i=1;i<=50;i++))
    do
        mkdir run${i}
        cp *.{in,prmtop,rst,inpcrd} run${i}
        for compound in LIG receptor LIGreceptor
        do
            auto_pbs.py ${compound} $(basename $arg)_run${i}_${compound}
            mv run.pbs run${i}/run_${compound}.pbs
        done
    done
done
