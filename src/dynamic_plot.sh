(
echo " plot 'profile_1-4.dat' ,  'profile_14.dat'  w l title 'test'";
while :; do sleep 5; echo replot; done
) | gnuplot -persist
