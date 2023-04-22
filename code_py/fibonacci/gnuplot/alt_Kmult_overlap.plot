set term pdf fsize 7 enhanced color dashed lw 1 size 4, 8


set datafile commentschars '%#'


set style line 1 lt 1 pt 8 ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 4 ps 0.9    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 10 ps 1   lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4 ps 1    lc rgb "#800080"         lw 2


set style line 11 dt 3 pt 8 ps 1    lc rgb "#D00040"         lw 2
set style line 22 dt 3 pt 4 ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 3 pt 10 ps 1   lc rgb "orange"          lw 2
set style line 44 dt 3 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 45 dt 3 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 66 dt 3 pt 12 ps 1   lc rgb  "blue"           lw 2
set style line 77 dt 3 pt 4 ps 1    lc rgb "#800080"         lw 2

set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 4
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'alt_Kmult_overlap.pdf'

set multiplot

set xlabel font "Times-Italic"
set key font "Times-Italic"
#set label font "Times-Italic"

LM1=0.08;  RM1=0.98
LM2=0.58;  RM2=0.98

BM2=0.04;  TM2=0.48
BM1=0.54;  TM1=0.98

#==========================================================================

set tmargin at screen TM1
set bmargin at screen BM1

set logscale y

set xlabel '{/Symbol s}_i' offset 0,0.6
set ylabel 'probability'   offset 0, 0

set ytics 10
set format y "10^{%L}"

set key r spacing 1.2 samplen 2

set xtics 2

#==========================================================================

set lmargin at screen LM1
set rmargin at screen RM1


set grid

plot [-5:5][1e-4:1] \
    '../runsQ/KM/q7_v0dt1_i08_Kmult.txt' u ($1):($2)  ev 1  w l ls 1  t '{/Symbol a} = 0',  \
    '../runsQ/KM/q7_v0dt1_i20_Kmult.txt' u ($1):($2)  ev 1  w l ls 2  t '{/Symbol a} = 0', \
    '../runsQ/KM/q7_v0dt1_i44_Kmult.txt' u ($1):($2)  ev 1  w l ls 3  t '{/Symbol a} = 0', \
    '../runsQ/KM/q7b_g3500dt5_i16_Kmult.txt' u ($1):($2) ev 1 w l ls 6 t '{/Symbol a} = 1', \
    '../runsQ/KM/q7b_g3500dt5_i40_Kmult.txt' u ($1):($2) ev 1 w l ls 4 t '{/Symbol a} = 1', \
    '../runsQ/KM/q7b_g3500dt5_i52_Kmult.txt' u ($1):($2) ev 1 w l ls 5 t '{/Symbol a} = 1'

#==========================================================================

unset logscale
set xtics auto
set ytics auto

set format y "%g"


set tmargin at screen TM2
set bmargin at screen BM2

set xlabel 'X / k' offset 0, 0
set ylabel 'ln (P(X)) / k' offset 1,0
set ylabel font "Times-Italic"


p(k) = k**(-1)
q(k) = k**(-1)

set ytics 0.1
set xtics 0.05

set key 


#------------------------
set lmargin at screen LM1
set rmargin at screen RM1

k1=42; k2=34; k3=26
q1=30; q2=38; q3=46


plot [0:0.18][-0.5:0] \
       '../runsQ/StichedPDF/q7_v0dt1_i08_prob.txt' \
    i 0   u (log($1)*p(k1)):(log($1*$2)*q(k1)) w lp ls 1 t '{/Symbol a} = 0', \
 '' i 3   u (log($1)*p(k1)):(log($1*$2)*q(k1)) w lp ls 1 notitle, \
 '' i 6   u (log($1)*p(k1)):(log($1*$2)*q(k1)) w lp ls 1 notitle, \
 '' i 9   u (log($1)*p(k1)):(log($1*$2)*q(k1)) w lp ls 1 notitle, \
       '../runsQ/StichedPDF/q7_v0dt1_i16_prob.txt' \
    i 0   u (log($1)*p(k2)):(log($1*$2)*q(k2)) w lp ls 3  t '{/Symbol a} = 0', \
 '' i 3   u (log($1)*p(k2)):(log($1*$2)*q(k2)) w lp ls 3 notitle, \
 '' i 6   u (log($1)*p(k2)):(log($1*$2)*q(k2)) w lp ls 3 notitle, \
 '' i 9   u (log($1)*p(k2)):(log($1*$2)*q(k2)) w lp ls 3 notitle, \
       '../runsQ/StichedPDF/q7_v0dt1_i24_prob.txt' \
    i 0   u (log($1)*p(k3)):(log($1*$2)*q(k3)) w lp ls 2  t '{/Symbol a} = 0', \
 '' i 3   u (log($1)*p(k3)):(log($1*$2)*q(k3)) w lp ls 2 notitle, \
 '' i 6   u (log($1)*p(k3)):(log($1*$2)*q(k3)) w lp ls 2 notitle, \
 '' i 9   u (log($1)*p(k3)):(log($1*$2)*q(k3)) w lp ls 2 notitle, \
       '../runsQ/StichedPDF/q7b_g3500dt5_i56_prob.txt' \
    i 0   u (log($1)*p(q3)):(log($1*$2)*q(q3)) w lp ls 6 t '{/Symbol a} = 1', \
 '' i 3   u (log($1)*p(q3)):(log($1*$2)*q(q3)) w lp ls 6 notitle, \
 '' i 6   u (log($1)*p(q3)):(log($1*$2)*q(q3)) w lp ls 6 notitle, \
 '' i 9   u (log($1)*p(q3)):(log($1*$2)*q(q3)) w lp ls 6 notitle, \
       '../runsQ/StichedPDF/q7b_g3500dt5_i48_prob.txt' \
    i 0   u (log($1)*p(q2)):(log($1*$2)*q(q2)) w lp ls 4  t '{/Symbol a} = 1', \
 '' i 3   u (log($1)*p(q2)):(log($1*$2)*q(q2)) w lp ls 4 notitle, \
 '' i 6   u (log($1)*p(q2)):(log($1*$2)*q(q2)) w lp ls 4 notitle, \
 '' i 9   u (log($1)*p(q2)):(log($1*$2)*q(q2)) w lp ls 4 notitle, \
       '../runsQ/StichedPDF/q7b_g3500dt5_i40_prob.txt' \
    i 0   u (log($1)*p(q1)):(log($1*$2)*q(q1)) w lp ls 5  t '{/Symbol a} = 1', \
 '' i 3   u (log($1)*p(q1)):(log($1*$2)*q(q1)) w lp ls 5 notitle, \
 '' i 6   u (log($1)*p(q1)):(log($1*$2)*q(q1)) w lp ls 5 notitle, \
 '' i 9   u (log($1)*p(q1)):(log($1*$2)*q(q1)) w lp ls 5 notitle



#==========================================================================



unset multiplot

set output








