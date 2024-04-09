set terminal epslatex standalone size 18cm,12cm color colortext linewidth 2
set style fill   solid 1.00 border lt -1
#set grid nopolar
set grid
ev2kj=96.4853905398362
set style data points
#set key right

set out "dftb_mlff_model2.tex"
set title "MLFF vs DFTB"

set style data points
#set key right
set xlabel 'DFTB (eV)'
set ylabel 'MLFF (eV)'


set style line 11 lt 1 lc rgb '#0072bd' # blue

set style line 12 lt 1 lc rgb '#d95319' # orange
set style line 13 lt 1 lc rgb '#edb120' # yellow
set style line 14 lt 1 lc rgb '#7e2f8e' # purple
set style line 15 lt 1 lc rgb '#77ac30' # green
set style line 16 lt 1 lc rgb '#4dbeee' # light-blue
set style line 17 lt 1 lc rgb '#a2142f' # red

set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 1 pointtype 6 pointsize 1.
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 1 pointtype 6 pointsize 1.
set style line 3 linecolor rgb '#a2142f' linetype 1 linewidth 1 pointtype 6 pointsize 1.
set style line 4 linecolor rgb '#77ac30' linetype 1 linewidth 1 pointtype 6 pointsize 1.
set style line 5 linecolor rgb '#d95319' linetype 1 linewidth 1 pointtype 6 pointsize 1.

pw1=-1568.307466
tb1=-1567.6470375079548
pw1=-1568.307466	
tb1=-1559.6326264649324
pw1=-1232.3453079
tb1=-1213.160618108204
tb1=-1229.116156350934


#plot 'ee' using ($1/12-pw1/12):($2/12-tb1/12) w p  notitle, x notitle
plot 'ee4' using ($1/64-pw1/64):($2/64-tb1/64) w p  notitle, x notitle

