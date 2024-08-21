reset
set datafile separator ','
set term pdf enhanced lw 1.0 
set output './resultados/apresentadora.pdf'
set xlabel '{/*1.4 tempo (dias)}'
set ylabel '{/*1.4 Concentration}'
set style data lines
set key left top
set key horizontal
set grid
#plot 'output_python.csv' u 1:3 title "Apresentadora naive"
#plot 'output_python.csv' u 1:4 title "Apresentadora madura"
#plot 'output_python.csv' u 1:5 title "T helper naive"
#plot 'output_python.csv' u 1:6 title "T helper efetora"
#plot 'output_python.csv' u 1:7 title "T killer naive"
#plot 'output_python.csv' u 1:8 title "T killer efetora"
#plot 'output_python.csv' u 1:9 title "celula B"
#plot 'output_python.csv' u 1:10 title "Plasma short"
#plot 'output_python.csv' u 1:11 title "Plasma long"
#plot 'output_python.csv' u 1:12 title "Celula B de memoria"
#plot 'output_python.csv' u 1:13 title "IgM"
#plot 'output_python.csv' u 1:14 title "IgG"
#plot 'output_python.csv' u 1:15 title "sistema imune inato" 

plot 'output.csv' u 1:3 title "Apresentadora naive"
plot 'output.csv' u 1:4 title "Apresentadora madura"
plot 'output.csv' u 1:5 title "T helper naive"
plot 'output.csv' u 1:6 title "T helper efetora"
plot 'output.csv' u 1:7 title "T killer naive"
plot 'output.csv' u 1:8 title "T killer efetora"
plot 'output.csv' u 1:9 title "celula B"
plot 'output.csv' u 1:10 title "Plasma short"
plot 'output.csv' u 1:11 title "Plasma long"
plot 'output.csv' u 1:12 title "Celula B de memoria"
plot 'output.csv' u 1:13 title "IgM"
plot 'output.csv' u 1:14 title "IgG"
plot 'output.csv' u 1:15 title "sistema imune inato"
