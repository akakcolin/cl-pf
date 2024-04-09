#!/bin/sh

# extract energies from output files (log_test) in folders run_* 
for i in run_*; do
cd $i && grep -wns "free  energy   TOTEN" OUTCAR | awk '{print $6}' > ../tmp1 && pwd | cut -d "/" -f9 | cut -d "_" -f2 > ../tmp2 && paste -d' ' ../tmp2 ../tmp1 > ../energies_$i
cd ../
done

# collect energies to one file
cat $(ls */energies_* | sort -n -t _ -k 2) > f_energies

# delete unnecesary files
rm tmp* energies_run_*
