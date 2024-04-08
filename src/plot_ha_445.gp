#!/usr/local/bin/gnuplot --persist
###########################################################################################################################
###########################################################################################################################
####                                                                                                                   ####
####   Gnuplot Script To plot Thermal_properties Output File From phonopy                                              ####
####                                                                                                                   ####
###########################################################################################################################

if  (ARGC != 3){print "\n       Arguments Error... ";
print "===================================================================================================================="
print "  Usage: phonpy-proplot --gnuplot thermal_propertiey.yaml > A_thermal_properties.gp"
print "  Usage: paste A_thermal_propertity-dft.gp"
print "  Usage: paste A_thermal_propertiy.gp B_thermal_properties.gp C_thermal_properties.gp > A_thermal_properties "
print "  Usage: gnuplot -c pthermal_propertity.gp name 8 filename                                                            " 
print "====================================================================================================================\n"
exit
}

############################################################################################################
############################################################################################################
file2plot=ARG1
file2plotd3=ARG2
file2plotlj=ARG3
print "\n================================================="
print "File name                        : ", file2plot
print "File name                        : ", file2plotd3
print "File name                        : ", file2plotlj
print "-------------------------------------------------\n"
############################################################################################################

set terminal epslatex standalone size 18cm,12cm color colortext
set output sprintf("relative_free_energy_dft_d3_lj_%s.tex",ARG1)
# Rank1
latt1= -163537.380890521
# Rank84
latt2= -163537.32574901
d3latt1=-6989.632
d3latt2=-6990.1891 
ljlatt1= -6952.0363
ljlatt2=-6952.4198
ev2kj=96.4853905398362
symm1=2
symm2=2
symm3=4
set style fill   solid 1.00 border lt -1
set grid nopolar
#set grid noxtics nomxtics ytics nomytics noztics nomztics nortics nomrtics \
# nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   lt 0 linecolor 0 linewidth 0.500,  lt 0 linecolor 0 linewidth 0.500
set style increment default

set style func linespoints

set key outside
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set yrange [ * : * ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback


set xlabel 'Temperature(K)'
set ylabel 'Relative free energy (kJ/mol)'
#set label 1 at 50, 1 sprintf(ARG1)

#set label 1 sprintf(ARG1) at 50, 2 font ",12 bold" textcolor rgb "black" 
#set label 2 '$k_x$' at graph 1.07, first 0 front center


set style line 11 lt 1 lc rgb '#0072bd' # blue

set style line 12 lt 1 lc rgb '#d95319' # orange
set style line 13 lt 1 lc rgb '#edb120' # yellow
set style line 14 lt 1 lc rgb '#7e2f8e' # purple
set style line 15 lt 1 lc rgb '#77ac30' # green
set style line 16 lt 1 lc rgb '#4dbeee' # light-blue
set style line 17 lt 1 lc rgb '#a2142f' # red

set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2 pointtype 6 pointsize 1.5
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 2 pointtype 6 pointsize 1.5
set style line 3 linecolor rgb '#a2142f' linetype 1 linewidth 2 pointtype 6 pointsize 1.5
set style line 4 linecolor rgb '#77ac30' linetype 1 linewidth 2 pointtype 6 pointsize 1.5
set style line 5 linecolor rgb '#d95319' linetype 1 linewidth 2 pointtype 6 pointsize 1.5


plot sprintf("%s", file2plot) using 1:($2-$2) w lp ls 5 pt 10  t "445",\
     sprintf("%s", file2plot) using 1:($7/4-$2/4)+(latt2-latt1)*ev2kj/4 w lp ls 1 pt 8 t "PBE-MBD 447",\
     sprintf("%s", file2plotd3) using 1:($7/4-$2/4)+(d3latt2-d3latt1)*ev2kj/4 w lp ls 3 pt 4 t "DFTB-D3 447",\
     sprintf("%s", file2plotlj) using 1:($7/4-$2/4)+(ljlatt2-ljlatt1)*ev2kj/4 w lp ls 4 pt 6 t "DFTB-LJ 447",\
     

print "          ...Done\n" 
print "\nWriting Output file ---> ",sprintf("relative_free_energy_dft_d3_lj_%s.tex",ARG1),"\n"
