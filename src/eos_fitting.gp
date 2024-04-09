E(V) = E0 + 9.0/8.0*B0*V0*((V0/V)**(2.0/3.0)-1)**2 + 9.0/16.0*B0*(B0P-4)* \
       V0*((V0/V)**(2.0/3.0)-1.0)**3.0 + R4*((V0/V)**(2.0/3.0)-1.0)**4.0

B0P   = 2.60567 ; B0 = 0.111821
V0     = 167.834 ; E0 = 0.1            ; R4  = -73.033

set terminal png
set output "eos2.png"
set grid
set xlabel 'Volume({\305}^3)
set ylabel 'E_{total} (eV)'
fit E(x) 'ev2.dat'  u 1:($2*13.6058) via B0P,B0,V0,E0,R4
plot 'ev2.dat' u 1:2 w  p pt 7 ps 1 

#plot E(x) w l,'e-v.dat' u 1:2 w  p pt 7 ps 1 
#
