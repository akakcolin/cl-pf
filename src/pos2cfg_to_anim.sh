for i in  0 1 2 3 4 5 6 7 8 9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 
do
cp pos_$i.vasp POSCAR
atomsk POSCAR $i.cfg
echo "$i.cfg pic/$i.png" >> scr_anim
done
 #convert -delay 10 -loop 0 pic/*.png animation.gif

