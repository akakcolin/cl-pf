set terminal epslatex standalone lw 2 color 11
set output 'fa_dftb.tex'
set style fill   solid 1.00 border lt -1
set xlabel "Favip DFTB energy(KJ/mol) " 
set ylabel 'DFTB energy (KJ/mol)'
set grid nopolar
set grid noxtics nomxtics ytics nomytics noztics nomztics nortics nomrtics \
 nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   lt 0 linecolor 0 linewidth 0.500,  lt 0 linecolor 0 linewidth 0.500
set style increment default
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set yrange [ * : * ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback
f(x) = m*x + b
f1(x) = m1*x + b1
fit f(x) 'fa.dat' using 1:3 via m, b
fit f1(x) 'fa.dat' using 1:2 via m1, b1
plot 'fa.dat' using 1:3 with points pt 7 title "cellopt2", 'fa.dat' using 1:2 with points pt 7 title "unit", f(x) title sprintf("Cellopt2 Fit y=%.2fx+%.2f", m,b), f1(x) title sprintf("Unit Fix y=%.2fx+%.2f", m1, b1)

