reset
set datafile separator ','
set term pdf enhanced lw 1.0
set output './resultados/Python/anticorposPython.pdf'
set xlabel '{/*1.4 tempo (dias)}'
set ylabel '{/*1.4 Concentration}'
set style data lines
set key left top
set key horizontal
set grid
set xtics 10
val =1.5
set arrow from graph 0, first val to graph 1, first val nohead dashtype 2 linewidth 2
plot 'output_python.csv' u 1:($14) title "IgM+IgG", './dados/anticorposPorcoInoculado.csv' u 1:2 w p lw 1 title 'dados experimentais'
#plot 'output_python.csv' u 1:($14+$15) title "IgM+IgG", './dados/anticorposPorcoInoculado.csv' u 1:2 w p lw 1 title 'dados experimentais'
