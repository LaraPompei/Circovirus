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
plot 'output_python.csv' u 1:($13+$14) title "IgM+IgG", './dados/anticorposPorcoInoculadoConvertido.csv' u 1:2 w p lw 1 title 'dados experimentais'
