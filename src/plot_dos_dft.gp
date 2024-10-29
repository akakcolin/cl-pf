#set terminal png size 1200,900 linewidth 3
set terminal epslatex standalone size 18cm,12cm color colortext linewidth 2
#set output sprintf("relative_free_energy_dft_%s.tex",ARG1)
set output 'phonoDOS.tex'
set style fill   solid 1.00 border lt -1
#set grid nopolar

#set grid noxtics nomxtics ytics nomytics noztics nomztics nortics nomrtics \
# nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
#set grid layerdefault   lt 0 linecolor 0 linewidth 0.500,  lt 0 linecolor 0 linewidth 0.500
set style increment default

set style func linespoints

set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set yrange [ * : * ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback

set linestyle 1 linecolor rgb 'blue'
set linestyle 2 linecolor rgb 'red'
set multiplot layout 3,1
# first plot top
set bmargin at screen 0.65
set tmargin at screen 0.95
set xtics format ' '
#set ytics 0,2,8

plot 'xxiii_dos_dft.dat' with lines linestyle 1 title "PBE-MBD"
# Second plot middle
set bmargin at screen 0.35
set tmargin at screen 0.65
plot 'xxiii_dos_d3.dat' with lines linestyle 1 title "DFTB-D3"
#plot 'benzene.dat' with lines linestyle 1 notitle,'fermi.dat' with lines linestyle 2 notitle

# Third plot
set bmargin at screen 0.05
set tmargin at screen 0.35
set xtics format '%g'
plot 'xxiii_dos_lj.dat' with lines linestyle 1 title "DFTB-LJ"
unset multiplot
set xlabel "Frequency (THz)"
