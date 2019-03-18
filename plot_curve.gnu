set style line 1 linecolor rgb "#DD0000" lw 3  pt 5 ps 1
set style line 2 linecolor rgb "#009900" lw 3  pt 7 ps 1
set style line 3 linecolor rgb "#0000BB" lw 3  pt 9 ps 1

set term pos eps enh color solid "Helvetica" 18
set encoding iso_8859_1
set out "melting_curves.eps"
set title 'DNA melting temperature'
set xlabel 'T (K)'
set ylabel '{/Helvetica-Oblique f}'  %f=1 fully hbonded; f=0 no hbond between strand
set xrange [250:400]
set yrange [0:1]
pl \
 'bre_melting_curve.out' u 1:2 t 'Breslauer' w lp ls 1, \
 'san_melting_curve.out' u 1:2 t 'SantaLucia' w lp ls 2, \
 'sug_melting_curve.out' u 1:2 t 'Sugimoto' w lp ls 3
exit
