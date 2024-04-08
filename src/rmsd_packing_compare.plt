set terminal epslatex standalone
set output 'Favip_cg_vasp_rmsd.tex'
set grid lt -1
set grid ytics
set grid xtics
set title "Structure rmsd value map (using Mercury)"
set ylabel "Favipiravir struct ID (VASP-OPT) " 
set xlabel "Favipiravir struct ID (cg) " 
set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)
#set palette rgbformulae -7, 2, -7
set cblabel "rmsd value (num matched mol > 13)" 
plot "cg_vasp_rmsd.dat"  with points pt 7 palette  notitle

