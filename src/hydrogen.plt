
plot 'epsilon.txt' index 0 title 'Calculated' with points pointtype 1 pointsize 1, \
	'epsilon.txt' index 1 title 'Extrapolated' with points pointtype 4 pointsize 1, \
	'epsilon.txt' index 2 title 'Interpolated' with lines 
set grid
set title '-Epsilon vs 1/Rho0'
set xlabel '1/Rho0'
set ylabel '-Epsilon'
set terminal gif large size 1024,760
set output 'Epsilon.gif'
replot


plot 'hydrogen10.txt' index 0 using 1:2 title 'u10' with lines linetype 0, \
	'hydrogen10.txt' index 0 using 1:3 title 'R10' with lines linetype 1, \
	'hydrogen10.txt' index 0 using 1:4 title 'P10' with lines linetype 2, \
	'hydrogen10.txt' index 1 title 'R10exact' with points pointtype 1 pointsize 1
set grid
set title 'Hydrogen state 10 (ground state)'
set xlabel 'rho'
set ylabel 'Dimensionless u, R, P'
set terminal gif large size 1024,760
set output 'Hydrogen10.gif'
replot


plot 'hydrogen20.txt' index 0 using 1:2 title 'u20' with lines linetype 0, \
	'hydrogen20.txt' index 0 using 1:3 title 'R20' with lines linetype 1, \
	'hydrogen20.txt' index 0 using 1:4 title 'P20' with lines linetype 2, \
	'hydrogen20.txt' index 1 title 'R20exact' with points pointtype 1 pointsize 1
set grid
set title 'Hydrogen state 20 (first excited state)'
set xlabel 'rho'
set ylabel 'Dimensionless u, R, P'
set terminal gif large size 1024,760
set output 'Hydrogen20.gif'
replot


plot 'hydrogen21.txt' index 0 using 1:2 title 'u21' with lines linetype 0, \
	'hydrogen21.txt' index 0 using 1:3 title 'R21' with lines linetype 1, \
	'hydrogen21.txt' index 0 using 1:4 title 'P21' with lines linetype 2, \
	'hydrogen21.txt' index 1 title 'R21exact' with points pointtype 1 pointsize 1
set grid
set title 'Hydrogen state 21'
set xlabel 'rho'
set ylabel 'Dimensionless u, R, P'
set terminal gif large size 1024,760
set output 'Hydrogen21.gif'
replot


plot 'hydrogen30.txt' index 0 using 1:2 title 'u30' with lines linetype 0, \
	'hydrogen30.txt' index 0 using 1:3 title 'R30' with lines linetype 1, \
	'hydrogen30.txt' index 0 using 1:4 title 'P30' with lines linetype 2, \
	'hydrogen30.txt' index 1 title 'R30exact' with points pointtype 1 pointsize 1
set grid
set title 'Hydrogen state 30'
set xlabel 'rho'
set ylabel 'Dimensionless u, R, P'
set terminal gif large size 1024,760
set output 'Hydrogen30.gif'
replot



plot 'hydrogen31.txt' index 0 using 1:2 title 'u31' with lines linetype 0, \
	'hydrogen31.txt' index 0 using 1:3 title 'R31' with lines linetype 1, \
	'hydrogen31.txt' index 0 using 1:4 title 'P31' with lines linetype 2, \
	'hydrogen31.txt' index 1 title 'R31exact' with points pointtype 1 pointsize 1
set grid
set title 'Hydrogen state 31'
set xlabel 'rho'
set ylabel 'Dimensionless u, R, P'
set terminal gif large size 1024,760
set output 'Hydrogen31.gif'
replot



plot 'hydrogen32.txt' index 0 using 1:2 title 'u32' with lines linetype 0, \
	'hydrogen32.txt' index 0 using 1:3 title 'R32' with lines linetype 1, \
	'hydrogen32.txt' index 0 using 1:4 title 'P32' with lines linetype 2, \
	'hydrogen32.txt' index 1 title 'R32exact' with points pointtype 1 pointsize 1
set grid
set title 'Hydrogen state 32'
set xlabel 'rho'
set ylabel 'Dimensionless u, R, P'
set terminal gif large size 1024,760
set output 'Hydrogen32.gif'
replot


