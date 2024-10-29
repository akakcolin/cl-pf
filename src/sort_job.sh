rm out_sorted.extxyz
for i in  qm_c0_12.extxyz qm_c0_18.extxyz qm_c0_20.extxyz qm_c0_22.extxyz qm_c0_23.extxyz qm_c0_24.extxyz qm_c0_5.extxyz
do
#grep energy $i | sort -r | awk '{print  $18}' >sorted_energy_id_1
grep energy $i | sort -r | awk '{print $11}' >sorted_energy
python gather_xyz.py $i sorted_energy
done
mv out_sorted.extxyz all.extxyz
grep energy all.extxyz | sort -r | awk '{print $11}' >sorted_energy
python gather_xyz.py all.extxyz sorted_energy

