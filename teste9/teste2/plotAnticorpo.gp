reset
set datafile separator ','
set term pdf enhanced lw 1.0
set output './resultados/anticorpos.pdf'
set xlabel '{/*1.4 tempo (dias)}'
set ylabel '{/*1.4 Concentration}'
set style data lines
set key left top
set key horizontal
set grid
plot 'output.csv' u 1:13, './dados/anticorposPorcoInoculado.csv' u 1:2 w p lw 1#u 1:2:3:4 w yerrorbars 
#plot 'output.csv' u 1:($12+$14), './dados/anticorposPorcoInoculado.csv' u 1:2 w p lw 1
