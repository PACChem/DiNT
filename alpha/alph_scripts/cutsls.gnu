set encoding iso_8859_1
set term postscript eps enhanced color "Times-Roman, 12"
set output "test.eps"

set xlabel "{/Times-Italic R}, {\305}" font "Times-Roman,16"
set ylabel "{/Times-Italic V}_I, cm^-^1" font "Times-Roman,16"

set label "(a)" at 1.8,2520 font "Times-Roman,16"
set label "(b)" at 1.8,-1020 font "Times-Roman,16"

set lmargin 10.

set format x "%5.1f" 
set format y "%5.0f"

set xtics font "Times-Roman,16"
set ytics font "Times-Roman,16"

set size 0.5, 1.
set origin 0., 0.

set multiplot
set xrange [3.:6.5]
set yrange [*:100]
set size 0.5, .5
set origin 0., 0.
set style data points
unset key
plot 'fit.dat' u 1:3 lw 1 lc -1 dt 1 smooth csplines title "Cut 1", \
     'fit.dat' u 1:5 lw 1 lc 1 dt 1 smooth csplines title "Cut 2", \
     'fit.dat' u 1:7 lw 1 lc 2 dt 1 smooth csplines title "Cut 3", \
     'fit.dat' u 1:9 lw 1 lc 3 dt 1 smooth csplines title "Cut 4", \
     'fit.dat' u 1:11 lw 1 lc 4 dt 1 smooth csplines title "Cut 5", \
     'fit.dat' u 1:13 lw 1 lc 5 dt 1 smooth csplines title "Cut 6", \
     'fit.dat' u 1:2 lw 2 lt 0 lc -1 smooth csplines title "", \
     'fit.dat' u 1:4 lw 2 lt 0 lc 1 smooth csplines title "", \
     'fit.dat' u 1:6 lw 2 lt 0 lc 2 smooth csplines title "", \
     'fit.dat' u 1:8 lw 2 lt 0 lc 3 smooth csplines title "", \
     'fit.dat' u 1:10 lw 2 lt 0 lc 4 smooth csplines title "", \
     'fit.dat' u 1:12 lw 2 lt 0 lc 5 smooth csplines title "", \
     'fit.dat' u 1:2 ps .75 lt 7 lc -1 title "", \
     'fit.dat' u 1:4 ps .75 lt 7 lc 1 title "", \
     'fit.dat' u 1:6 ps .75 lt 7 lc 2 title "", \
     'fit.dat' u 1:8 ps .75 lt 7 lc 3 title "", \
     'fit.dat' u 1:10 ps .75 lt 7 lc 4 title "", \
     'fit.dat' u 1:12 ps .75 lt 7 lc 5 title ""

set xrange [2.5:4.5]
set yrange [*:2500]
set size 0.5, .5
set origin 0., .5
#set key top right maxrows 6 width -7
set key top right 
unset xlabel
plot 'fit.dat' u 1:3 lw 1 lc -1 dt 1 smooth csplines title "Cut 1", \
     'fit.dat' u 1:5 lw 1 lc 1 dt 1 smooth csplines title "Cut 2", \
     'fit.dat' u 1:7 lw 1 lc 2 dt 1 smooth csplines title "Cut 3", \
     'fit.dat' u 1:9 lw 1 lc 3 dt 1 smooth csplines title "Cut 4", \
     'fit.dat' u 1:11 lw 1 lc 4 dt 1 smooth csplines title "Cut 5", \
     'fit.dat' u 1:13 lw 1 lc 5 dt 1 smooth csplines title "Cut 6", \
     'fit.dat' u 1:2 lw 2 lt 0 lc -1 smooth csplines title "", \
     'fit.dat' u 1:4 lw 2 lt 0 lc 1 smooth csplines title "", \
     'fit.dat' u 1:6 lw 2 lt 0 lc 2 smooth csplines title "", \
     'fit.dat' u 1:8 lw 2 lt 0 lc 3 smooth csplines title "", \
     'fit.dat' u 1:10 lw 2 lt 0 lc 4 smooth csplines title "", \
     'fit.dat' u 1:12 lw 2 lt 0 lc 5 smooth csplines title "", \
     'fit.dat' u 1:2 ps .75 lt 7 lc -1 title "", \
     'fit.dat' u 1:4 ps .75 lt 7 lc 1 title "", \
     'fit.dat' u 1:6 ps .75 lt 7 lc 2 title "", \
     'fit.dat' u 1:8 ps .75 lt 7 lc 3 title "", \
     'fit.dat' u 1:10 ps .75 lt 7 lc 4 title "", \
     'fit.dat' u 1:12 ps .75 lt 7 lc 5 title ""

