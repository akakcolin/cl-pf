set terminal png
set output "qhq_1324.png"
set xrange [17:0]
set yrange [17:0]
set xtics 1.0
set ytics 1.0
set grid
#set palette model RGB rgbformulae 7,5,15 negative
#set palette model RGB defined (0 "white", 1 "yellow", 2 "dark-yellow", 2 "red", 3 "dark-red" )
#
#set palette defined ( 0 '#000090',1 '#000fff', 2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000')
#set palette defined (0 1 1 1, 1 0 0 0)
#set palette defined ( 0 "#e0e0e0", 1 "#00ff00", 2 "#ff0000", 3 "#101010" ) 
#set palette model RGB defined ( 0 "#e0e0e0", 1 "#00ff00", 2 "#ff0000" ) 
#set palette model RGB defined ( 0 "white", 1 "#e0e0e0", 2 "#000090", 3 "#ff0000" ) 
set palette model RGB defined ( 0 "white", 1 "#000090", 2 "#ff0000" ) 
set cbrange [0:0.001]

plot 'qhq.txt' using 2:1:3 with image
