set datafile separator ','

set key autotitle columnhead

set ylabel "\\rho_{\\nu} (GeV)"
set xlabel "Temperature (MeV)"
set logscale x 10
set xrange[15:0.0007]
set logscale y 10
set yrange[1e-25:1e-5]

set xtics nomirror
set ytics nomirror
set key right top

set style line 1 lc rgb "#e6402e" lw 2 #red for rho_nu
set style line 2 lc rgb "#e6872e" lw 2 #orange for (d/dTnu)(rho_nu)
set style line 3 lc rgb "#2e56e6" lw 2 #blue for rho_nu_vs
set style line 4 lc rgb "#2ee6e6" lw 2 #light blue for (d/dTnu)(rho_nu_vs)

plot 'evolution.csv' using 3:7 with lines ls 1, 'evolution.csv' using 3:8 with lines ls 2, 'evolution_vs.csv' using 3:7 with lines ls 3, 'evolution_vs.csv' using 3:8 with lines ls 4
#plot 'evolution.csv' using 3:7 with lines ls 1, 'evolution_vs.csv' using 3:7 with lines ls 3
