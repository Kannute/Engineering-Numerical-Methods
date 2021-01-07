set term png
set grid

#zad1
set output "1_1.png"
set xlabel "t"
set ylabel "y"
set title "RK2"
plot "zad1_1.dat" i 0 using 1:2 w l lw 2 t "dt= 0.01",\
    "zad1_1.dat" i 1 using 1:2 w p t "dt= 0.1", \
    "zad1_1.dat" i 2 using 1:2 w p t "dt= 1.0", \
    "zad1_1.dat" i 1 using 1:2 w p t "dok"

set output "1_2.png"
set xlabel "t"
set ylabel "|yn-ya|"
set title "RK2"
plot "zad1_2.dat"i 0 using 1:2 w l lw 2 t "dt = 0.01", \
    "zad1_2.dat"i 1 using 1:2 w l t "dt = 0.1", \
    "zad1_2.dat" i 2 using 1:2 w l t "dt = 1.0"

#zad2


set output "2_1.png"
set xlabel "t"
set ylabel "y"
set title "RK2 - trapez"
plot "zad2_1.dat" i 0 using 1:2 w l lw 2 t "dt= 0.01",\
    "zad2_1.dat" i 1 using 1:2 w p t "dt= 0.1", \
    "zad2_1.dat" i 2 using 1:2 w p t "dt= 1.0", \
    "zad2_1.dat" i 1 using 1:2 w p t "dok"

set output "2_2.png"
set xlabel "t"
set ylabel "|yn-ya|"
set title "RK2"
plot "zad2_2.dat"i 0 using 1:2 w l lw 2 t "dt = 0.01", \
    "zad2_2.dat"i 1 using 1:2 w l t "dt = 0.1", \
    "zad2_2.dat" i 2 using 1:2 w l t "dt = 1.0"

#zad3_1

set output "3_1.png"
set xlabel "t"
set ylabel "y"
set title "RK4"
plot "zad3_1.dat" i 0 using 1:2 w l lw 2 t "dt = 0.01", \
    "zad3_1.dat" i 1 using 1:2 w p t "dt = 0.1", \
    "zad3_1.dat" i 2 using 1:2 w p t "dt = 1.0", \
    "zad3_1.dat" i 1 using 1:2 w p t "dok"


set output "3_2.png"
set xlabel "t"
set ylabel "|yn-ya|"
set title "RK4"
plot "zad3_2.dat"i 0 using 1:2 w l lw 2 t "dt = 0.01", \
    "zad3_2.dat"i 1 using 1:2 w l t "dt = 0.1", \
    "zad3_2.dat" i 2 using 1:2 w l t "dt = 1.0"


#zad4_1

set output "4_1.png"
set xlabel "t"
set ylabel "Q"
set title "z4 RLC Q"
set xrange [0.0:0.25]
set yrange [-0.002:0.0035]
plot "zad4_1.dat" i 0 using 1:2 w l lw 2 t "0.5 w_0", \
      "zad4_1.dat" i 1 using 1:2 w l lw 2 t "0.8 w_0", \
      "zad4_1.dat" i 2 using 1:2 w l lw 2 t "1.0 w_0", \
      "zad4_1.dat" i 3 using 1:2 w l lw 2 t "1.2 w_0"


set output "4_2.png"
set xlabel "t"
set ylabel "I"
set title "z4 RLC I"
set xrange [0.0:0.25]
set yrange [-0.15:0.125]
plot "zad4_2.dat" i 0 using 1:2 w l lw 2 t "0.5 w_0", \
      "zad4_2.dat" i 1 using 1:2 w l lw 2 t "0.8 w_0", \
      "zad4_2.dat" i 2 using 1:2 w l lw 2 t "1.0 w_0", \
      "zad4_2.dat" i 3 using 1:2 w l lw 2 t "1.2 w_0"
