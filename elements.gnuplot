set datafile separator ','

set key autotitle columnhead

set ylabel "Relative Abuundance"
set xlabel "Element"
set xrange[0.9:6.1]
set logscale y 10
set yrange[1e-15:1]

set xtics ("Yp" 1, "H2/H" 2, "He3/H" 3, "Li7/H" 4, "Li6/H" 5, "He7/H" 6)
set xtics nomirror
set ytics nomirror

set style fill transparent solid 1 #noborder
set style circle radius 0.005

#plot 'stand_cosmo.csv' using 1:3:4:2 with errorbars 
plot 'alter_vs.csv' using 1:3:4:2 with errorbars

#plot 'stand_cosmo.csv' using 1:2 with circles, '' using 1:3 with circles, '' using 1:4 with circles, '' using 1:5 with circles, '' using 1:6 with circles, '' using 1:7 with circles
