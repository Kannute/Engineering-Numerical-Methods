set term png
set grid

set title "RK2 x(t)"
set output "RK2_1.png"
set xlabel "t"
set ylabel "x(t)"


plot "zad1.dat" i 0 u 1:3 w l lw 1 t "TOL 10^{-2}", "" i 1 u 1:3 w l lw 1 t "TOL 10^{-5}"


set title "RK2 v(t)"
set output "RK2_2.png"
set xlabel "t"
set ylabel "v(t)"


plot "zad1.dat" i 0 u 1:4 w l lw 1 t "TOL 10^{-2}", "" i 1 u 1:4 w l lw 1 t "TOL 10^{-5}"


set title "RK2 dt(t)"
set output "RK2_3.png"
set xlabel "t"
set ylabel "dt(t)"


plot "zad1.dat" i 0 u 1:2 w l lw 1 t "TOL 10^{-2}", "" i 1 u 1:2 w l lw 1 t "TOL 10^{-5}"

set title "RK2 v(x)"
set output "RK2_4.png"
set xlabel "x"
set ylabel "v(x)"


plot "zad1.dat" i 0 u 3:4 w l lw 1 t "TOL 10^{-2}", "" i 1 u 3:4 w l lw 1 t "TOL 10^{-5}"


set title "MT x(t)"
set output "MT_1.png"
set xlabel "t"
set ylabel "x(t)"


plot "zad2.dat" i 0 u 1:3 w l lw 1 t "TOL 10^{-2}", "" i 1 u 1:3 w l lw 1 t "TOL 10^{-5}"


set title "MT v(t)"
set output "MT_2.png"
set xlabel "t"
set ylabel "v(t)"


plot "zad2.dat" i 0 u 1:4 w l lw 1 t "TOL 10^{-2}", "" i 1 u 1:4 w l lw 1 t "TOL 10^{-5}"


set title "MT dt(t)"
set output "MT_3.png"
set xlabel "t"
set ylabel "dt(t)"


plot "zad2.dat" i 0 u 1:2 w l lw 1 t "TOL 10^{-2}", "" i 1 u 1:2 w l lw 1 t "TOL 10^{-5}"


set title "MT v(x)"
set output "MT_4.png"
set xlabel "x"
set ylabel "v(x)"


plot "zad2.dat" i 0 u 3:4 w l lw 1 t "TOL 10^{-2}", "" i 1 u 3:4 w l lw 1 t "TOL 10^{-5}"