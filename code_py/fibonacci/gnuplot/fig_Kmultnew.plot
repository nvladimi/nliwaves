set term pdf fsize 7 enhanced color dashed lw 1 size 4, 2


set datafile commentschars '%#'


set style line 1 lt 1 pt  8 ps 1.2   lc rgb "#D00040"         lw 3
set style line 2 lt 1 pt  4 ps 1.1   lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 10 ps 1.3   lc rgb "orange"          lw 2
set style line 4 lt 1 pt  4 ps 1.2   lc rgb "#60B000"         lw 3
set style line 5 lt 1 pt  6 ps 1.2   lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1.2   lc rgb  "blue"           lw 3
set style line 7 lt 1 pt  4 ps 1     lc rgb "#800080"         lw 2


set style line 11 dt 3 pt 8 ps 1    lc rgb "#D00040"         lw 2
set style line 22 dt 3 pt 4 ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 3 pt 10 ps 1   lc rgb "orange"          lw 2
set style line 44 dt 3 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 45 dt 3 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 66 dt 3 pt 12 ps 1   lc rgb  "blue"           lw 2
set style line 77 dt 3 pt 4 ps 1    lc rgb "#800080"         lw 2

set style line 10   dt 1 pt 7 ps 0.5    lc rgb "black"       lw 4
set style line 100  dt 3 pt 6 ps 1    lc rgb "black"       lw 4
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'fig_Kmultnew.pdf'

set multiplot

set xlabel font "Times-Italic"
set key font "Times-Italic"
#set label font "Times-Italic"

LM1=0.08;  RM1=0.48
LM2=0.58;  RM2=0.98

BM1=0.16;  TM1=0.96

#==========================================================================

set tmargin at screen TM1
set bmargin at screen BM1

set logscale y

set xlabel '{/Symbol s}_i' offset 0,0.6
set ylabel 'probability'   offset 0, 0

set ytics 100
set format y "10^{%L}"

set key r spacing 1.2 samplen 2

set xtics 4


#==========================================================================

set lmargin at screen LM1
set rmargin at screen RM1

set label 1 "{/Symbol a = 0}" at -7, 0.3

savg = -0.1604039416865561   # -1/3 ln (phi)
sigma = 1

f(x) = 0.5 /(cosh(x - savg)/sigma)**2

#set xtics 2; set grid; set key b c

plot [-8:8][1e-6:1] \
    f(x) w l ls 10 notitle, \
    '../runsQ/KM/q7_v0dt1_i08_Kmult.txt' u ($1):($2)  ev 6    w p ls 1 t 'M - i = 42',  \
    '../runsQ/KM/q7_v0dt1_i20_Kmult.txt' u ($1):($2)  ev 6::2 w p ls 4 t '30', \
    '../runsQ/KM/q7_v0dt1_i44_Kmult.txt' u ($1):($2)  ev 6::4 w p ls 6 t '6' #, \
    exp(2*x) w l ls 100


unset grid
    
#==========================================================================

set lmargin at screen LM2
set rmargin at screen RM2

set label 1 "{/Symbol a = 1}"

savg = -0.3208078833731123   # -2/3 ln (phi)
sigma = 1

f(x) = 0.5 / (cosh(x - savg)/sigma)**2

plot [-8:8][1e-6:1] \
    f(x) w l ls 10 notitle, \
    '../runsQ/KM/q7b_g3500dt5_i16_Kmult.txt' u ($1):($2) ev 6    w p ls 1 t 'i - M = 6', \
    '../runsQ/KM/q7b_g3500dt5_i40_Kmult.txt' u ($1):($2) ev 6::2 w p ls 4 t '30', \
    '../runsQ/KM/q7b_g3500dt5_i52_Kmult.txt' u ($1):($2) ev 6::4 w p ls 6 t '42'
    
#==========================================================================



#==========================================================================



unset multiplot

set output








