#!/usr/local/bin/gnuplot --persist
###########################################################################################################################
###########################################################################################################################
####                                                                                                                   ####
####   Gnuplot Script To plot Thermal_properties Output File From phonopy                                              ####
####                                                                                                                   ####
###########################################################################################################################

if  (ARGC != 1){print "\n       Arguments Error... ";
print "===================================================================================================================="
print "====================================================================================================================\n"
exit
}

############################################################################################################
############################################################################################################
file2plot=ARG1
print "\n================================================="
print "File name                        : ", file2plot
print "-------------------------------------------------\n"
############################################################################################################

set terminal epslatex standalone size 18cm,12cm color colortext
set output sprintf("opt_energy_force_%s.tex",ARG1)
set grid nopolar
set grid xtics nomxtics noytics nomytics noztics nomztics nortics nomrtics \
 nox2tics nomx2tics y2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   lt 0 linecolor 0 linewidth 0.500,  lt 0 linecolor 0 linewidth 0.500
set key fixed center top vertical Right noreverse enhanced autotitle nobox
set style data lines
set xtics border out scale 1,0.5 mirror norotate  autojustify
set xtics  norangelimit logscale autofreq 
set ytics border out scale 1,0.5 nomirror norotate  autojustify
set ytics  norangelimit logscale autofreq 
set ztics border out scale 1,0.5 nomirror norotate  autojustify
set cbtics border out scale 1,0.5 mirror norotate  autojustify
set rtics axis out scale 1,0.5 nomirror norotate  autojustify

set xlabel 'VASP-ENCUT'
set xrange [ * : *] noreverse nowriteback

set ylabel 'Relative Energy (kj/mol)'
set yrange [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

plot  sprintf("%s", file2plot)  using 1:($2-$4)*96.4869 w lp title "relative energy" ,\

print "          ...Done\n" 
print "\nWriting Output file ---> ",sprintf("opt_energy_force_%s.tex",ARG1),"\n"
