set datafile separator ','

set key autotitle columnhead

set ylabel "rho_{vs} (MeV)"
set xlabel "Temperature (MeV)"
set logscale x 10
set xrange[15:0.6]
set logscale y 10
set yrange[1e-40:1e5]

set xtics nomirror
set ytics nomirror
set key right top

set style line 1 lc rgb "#e6402e" lw 2 #red for rho_vs
set style line 2 lc rgb "#e6872e" lw 2 #orange for (d/dTnu)(rho_nu)
set style line 3 lc rgb "#2e56e6" lw 2 #blue for rho_nu_vs
set style line 4 lc rgb "#2ee6e6" lw 2 #light blue for (d/dTnu)(rho_nu_vs)

plot 'evolution_vs.csv' using 3:10 with lines ls 1
