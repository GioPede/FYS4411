# Definizioni Generali di Stile

set mxtics 10
set mytics 10
set key top left
set grid
set bars small
set terminal pngcairo enhanced font 'Arial, 20' background "#ffffff" fontscale 1.0 dashed size 1280,1280

set style line 1  linetype 1 linecolor rgb "red"  linewidth 2 
set style line 3  linetype 1 linecolor rgb "yellow"  linewidth 2 
set style line 4  linetype 1 linecolor rgb "green"  linewidth 2 
set style line 6  linetype 1 linecolor rgb "blue"  linewidth 2 
set style line 5  linetype 1 linecolor rgb "cyan"  linewidth 2 
set style line 2  linetype 1 linecolor rgb "orange"  linewidth 2 
set style line 7  linetype 1 linecolor rgb "black"  linewidth 2

# Parametri Specifici Del Grafico
set output 'plot_jastrow_ellipse.png'
#set xformat "%.0f"
#set yformat "%.0f"
set xrange [0.1:0.75]
#set yrange [0:2]
#set logscale x 
#set logscale y
set xlabel "{/Symbol a}"
set ylabel "E_{Interacting}/E_{Non-Interacting}"
set title "Ellyptical Trap Interaction Contribution"

# Plot del modulo
plot  '10.csv' u ($1):($2):($3) ls 1 w yerrorbars ti sprintf("10 particles"), \
	  '50.csv' u ($1):($2):($3) ls 6 w yerrorbars ti sprintf("50 particles"),\
	  '100.csv' u ($1):($2):($3) ls 4 w yerrorbars ti sprintf("100 particle")
