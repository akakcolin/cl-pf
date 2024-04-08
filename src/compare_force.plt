set terminal png 
set output 'cg_c2.png'
plot 'dumpto1/force.log' with lines, 'dumpto2/force.log' with lines
