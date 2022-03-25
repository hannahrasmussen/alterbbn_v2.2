set datafile separator ','

set key autotitle columnhead

set ylabel "Relative Abundance"
#set xlabel "Element"
set xlabel "Temperature (MeV)"
set logscale x 10
set xrange[2.33:0.00086]
set logscale y 10
set yrange[1e-15:1]

#set xtics ("Y(n)" 1, "Y(p)" 2, "2H" 3, "3H" 4, "3He" 5, "4He" 6, "6Li" 7, "7Li" 8, "7Be" 9)
set xtics nomirror
set ytics nomirror
set key left center

#set style fill transparent solid 1 #noborder
#set style circle radius 0.005

set style line 1 lc rgb "#e6402e" lw 2 #red for Y(n)
set style line 2 lc rgb "#e6872e" lw 2 #orange for Y(p) aka H
set style line 3 lc rgb "#e6d32e" lw 2 #gold for Y(2H)
set style line 4 lc rgb "#37e62e" lw 2 #lime green for Y(3H)
set style line 5 lc rgb "#09ba47" lw 2 #green for Y(3He)
set style line 6 lc rgb "#2ee6e6" lw 2 #light blue for Y(4He)
set style line 7 lc rgb "#2e56e6" lw 2 #blue for Y(6Li)
set style line 8 lc rgb "#782ee6" lw 2 #purple for Y(7Li)
set style line 9 lc rgb "#f763c1" lw 2 #pink for Y(7Be)

#set terminal pngcairo size 800,600 enhanced font 'Segoe UI,10'
#set output 'stand_cosmo.png'

#plot 'stand_cosmo.csv' using 1:3:4:2 with errorbars
#plot 'alter_vs.csv' using 1:3:4:2 with errorbars

#plot 'stand_cosmo.csv' using 1:2 with circles, '' using 1:3 with circles, '' using 1:4 with circles, '' using 1:5 with circles, '' using 1:6 with circles, '' using 1:7 with circles

plot 'evolution_vs.csv' using 3:12 with lines ls 1, '' using 3:13 with lines ls 2, '' using 3:14 with lines ls 3, '' using 3:15 with lines ls 4, '' using 3:16 with lines ls 5, '' using 3:17 with lines ls 6, '' using 3:18 with lines ls 7, '' using 3:19 with lines ls 8, '' using 3:20 with lines ls 9
#plot 'evolution.csv' using 3:11 with lines ls 1, '' using 3:12 with lines ls 2, '' using 3:13 with lines ls 3, '' using 3:14 with lines ls 4, '' using 3:15 with lines ls 5, '' using 3:16 with lines ls 6, '' using 3:17 with lines ls 7, '' using 3:18 with lines ls 8, '' using 3:19 with lines ls 9
