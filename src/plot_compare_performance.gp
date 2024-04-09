set terminal epslatex standalone size 18cm,12cm color colortext linewidth 2
set style fill   solid 1.00 border lt -1
set grid nopolar
#set xtics rotate by 45 right


set out "performance.tex"
set title ""
set key left
set xlabel 'Number of Atoms'
set ylabel 'Single Step Times(s)'

#set key outside right center
set yrange[0:400]
set style data histogram

set autoscale fix

set style data histogram
set style histogram cluster gap 1
set style fill solid
set boxwidth 0.9
set boxwidth 1 relative
set style fill solid 1.0 border -1

plot 'ddd' using 2:xtic(1) t "4cpu" ,'ddd' using 3 t "4cpu+3090" ,'ddd' using 4 t "4cpu+3090(ours)"

