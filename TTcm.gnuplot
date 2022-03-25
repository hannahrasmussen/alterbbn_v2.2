set datafile separator ','

set key autotitle columnhead

set ylabel "Plasma Temp (MeV)"
set xlabel "Comoving Temp(MeV)"
set logscale x 10
set xrange[15:0.0007]
set logscale y 10
set yrange[1e-3:15]

set xtics nomirror
set ytics nomirror
set key right top

set style line 1 lc rgb "#848587" lw 2 #grey for Tcm vs Tcm
set style line 2 lc rgb "#e6d32e" lw 2 #gold for T vs Tcm (std cosmo)
set style line 3 lc rgb "#782ee6" lw 2 #purple for T vs Tcm (our model)
set style line 4 lc rgb "#e6872e" lw 2 #orange for T vs Tcm (phi model)

plot 'evolution.csv' using 4:3 with lines ls 2, 'evolution_vs.csv' using 4:3 with lines ls 3, 'evolution_phi.csv' using 4:3 with lines ls 4, 'evolution_vs.csv' using 4:4 with lines ls 1

#plot 'evolution.csv' using 4:3 with lines ls 2, 'evolution.csv' using 4:4 with lines ls 1
