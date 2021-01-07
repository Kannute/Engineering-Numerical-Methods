set term png
set grid

set output "zad1.png"
set xlabel "t"
set ylabel "u(t), z(t)"
set title "Niejawna metoda trapezów z iteracją Picarda"

plot "zad1.dat" u 1:2 w l lw 2 t "u(t)", "" u 1:3 w l lw 2 t "z(t)"

set output "zad2.png"
set xlabel "t"
set ylabel "u(t), z(t)"
set title "Niejawna metoda trapezów z iteracją Newtona"

plot "zad2.dat" u 1:2 w l lw 2 t "u(t)", "" u 1:3 w l lw 2 t "z(t)"

set output "zad3.png"
set xlabel "t"
set ylabel "u(t), z(t)"
set title "Niejawna metoda RK2 max it = 3"

plot "zad3.dat" u 1:2 w l lw 2 t "u(t)", "" u 1:3 w l lw 2 t "z(t)"