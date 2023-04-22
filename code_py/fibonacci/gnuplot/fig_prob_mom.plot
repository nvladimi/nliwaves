set term pdf fsize 7 enhanced color dashed lw 1 size 4, 4

set datafile commentschars '%#'

PS = 0.7

set style line 1 lt 1 pt  8 ps PS    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt  6 ps PS    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 10 ps PS    lc rgb "orange"          lw 2
set style line 4 lt 1 pt  4 ps PS    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt  6 ps PS    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps PS    lc rgb  "blue"           lw 2
set style line 7 lt 1 pt  4 ps PS    lc rgb "#800080"         lw 2

set style line 11 dt 3 pt  8 ps PS    lc rgb "#D00040"         lw 2
set style line 22 dt 3 pt  6 ps PS    lc rgb "#FF4000"         lw 2
set style line 33 dt 3 pt 10 ps PS    lc rgb "orange"          lw 2
set style line 44 dt 3 pt  4 ps PS    lc rgb "#60B000"         lw 2
set style line 45 dt 3 pt  6 ps PS    lc rgb "#008080"         lw 2
set style line 66 dt 3 pt 12 ps PS    lc rgb  "blue"           lw 2
set style line 77 dt 3 pt  4 ps PS    lc rgb "#800080"         lw 2

set style line 10   dt 2 pt 7 ps 0.5  lc rgb "black"       lw 3
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1

#-------------------------

set output 'fig_prob_mom.pdf'

set multiplot


set xlabel font "Times-Italic"
set key font "Times-Italic"


unset key

set logscale y
set key r spacing 1 samplen  2

LM1=0.10;  RM1=0.50
LM2=0.58;  RM2=0.98

BM1=0.08;  TM1=0.48
BM2=0.58;  TM2=0.98


#==========================================================================

set tmargin at screen TM2
set bmargin at screen BM2

set xlabel '{/Symbol |} a_i {/Symbol |}^2 / n_i' offset 0, 0.5
set ylabel 'probability' offset -1,0

set format y "10^{%L}"

set ytics 100
set xtics 10

#------------------------
set lmargin at screen LM1
set rmargin at screen RM1

set label 1 "{/Symbol a = 0}" at 5,0.4
set label 2 "p - i" at 1.4, 1e-3  font "Times-Italic"
set key b l 

plot [0:30][1e-6:1] \
       '../runsQ/MIx32x1/q7_v0dt1_i08_prob.txt' u ($1):($2) w lp ls 6 t '44', \
       '../runsQ/MIx32x1/q7_v0dt1_i16_prob.txt' u ($1):($2) w lp ls 1 t '36', \
       '../runsQ/MIx32x1/q7_v0dt1_i24_prob.txt' u ($1):($2) w lp ls 2 t '28', \
       '../runsQ/MIx32x1/q7_v0dt1_i32_prob.txt' u ($1):($2) w lp ls 3 t '20', \
       '../runsQ/MIx32x1/q7_v0dt1_i40_prob.txt' u ($1):($2) w lp ls 4 t '12', \
       '../runsQ/MIx32x1/q7_v0dt1_i48_prob.txt' u ($1):($2) w lp ls 5 t '4', \
exp(-x) w l ls 10 notitle

#------------------------

set lmargin at screen LM2
set rmargin at screen RM2

set label 1 "{/Symbol a = 1}"
unset ylabel
set key b l

set label 2 "i - p"

plot [0:30][1e-6:1] \
       '../runsQ/MIx32x1/q7b_g3500dt5_i16_prob.txt' u ($1):($2) w lp ls 5 t '4', \
       '../runsQ/MIx32x1/q7b_g3500dt5_i24_prob.txt' u ($1):($2) w lp ls 4 t '12', \
       '../runsQ/MIx32x1/q7b_g3500dt5_i32_prob.txt' u ($1):($2) w lp ls 3 t '20', \
       '../runsQ/MIx32x1/q7b_g3500dt5_i40_prob.txt' u ($1):($2) w lp ls 2 t '28', \
       '../runsQ/MIx32x1/q7b_g3500dt5_i48_prob.txt' u ($1):($2) w lp ls 1 t '32', \
       '../runsQ/MIx32x1/q7b_g3500dt5_i56_prob.txt' u ($1):($2) w lp ls 6 t '44', \
exp(-x) w l ls 10 notitle

#==========================================================================

set tmargin at screen TM1
set bmargin at screen BM1

set ylabel font "Times-Italic"

set ylabel '{/Symbol |} a_i {/Symbol |}^4 / 2 n_i^2,   {/Symbol |} a_i {/Symbol |} ^6 / 6 n_i^3' offset -2.5,0

set ytics (0.5, 1, 2, 4, 8)
set format y "%g"
set key b r

#------------------------

set lmargin at screen LM1
set rmargin at screen RM1

set label 1 "{/Symbol a = 0}" at -2,8

set xlabel '(p + 1) - i'

plot [-10:50][0.75:10] \
 '../runsQ/POST5/q5_v0m30_mom.txt'  u (-$1+31):($4)     w lp ls 6  t 'p = 30 of 40', \
 '../runsQ/POST5/q5_v0m30_mom.txt'  u (-$1+31):($5)     w lp ls 6  notitle, \
 '../runsQ/POSTmts/q7_v0dt1_mts.txt' u (-$1+51):($6/2/($4)**2) w lp ls 1 t 'p = 50 of 60', \
 '../runsQ/POSTmts/q7_v0dt1_mts.txt' u (-$1+51):($8/6/($4)**3) w lp ls 1 notitle, \
 exp( 0.017*x ) w l ls 10 notitle, \
 exp( 0.050*x ) w l ls 10 notitle


#-------------------------

set lmargin at screen LM2
set rmargin at screen RM2

set label 1 "{/Symbol a = 1}"

set xlabel 'i - (p + 1)'
unset ylabel

plot [-10:50][0.75:10] \
 '../runsQ/POST6/q6_m10gr140_mom.txt'    u ($1-11):($4)   w lp ls 6  t 'p = 10 of 40', \
 '../runsQ/POST6/q6_m10gr140_mom.txt'    u ($1-11):($5)   w lp ls 6  notitle, \
 '../runsQ/POSTmts/q7b_g3500dt5_mts.txt' u ($1-11):($6/2/($4)**2) w lp ls 1 t 'p = 10 of 60', \
 '../runsQ/POSTmts/q7b_g3500dt5_mts.txt' u ($1-11):($8/6/($4)**3) w lp ls 1 notitle, \
exp( 0.015*x ) w l ls 10 notitle, \
exp( 0.043*x ) w l ls 10 notitle

unset label

#==========================================================================

set tmargin at screen 0.96
set bmargin at screen 0.80

set xlabel '{/Symbol q}_i {/Symbol / p}' offset 3, 1
unset ylabel

set format y "%g"
unset logscale y

set ytics 0.1 offset 0.2,0.2 font 'Helvetica, 6'
set xtics 1 offset 0, 0.3 font 'Helvetica, 6'

#------------------------
set lmargin at screen 0.32
set rmargin at screen 0.48

unset key

#      '../runsQ/MIx32x1/q7_v0dt1_i08_prob.txt' u ($1):($2) w lp ls 6 t '44', \
       '../runsQ/MIx32x1/q7_v0dt1_i16_prob.txt' u ($1):($2) w lp ls 1 t '36', \
       '../runsQ/MIx32x1/q7_v0dt1_i24_prob.txt' u ($1):($2) w lp ls 2 t '28', \
       '../runsQ/MIx32x1/q7_v0dt1_i32_prob.txt' u ($1):($2) w lp ls 3 t '20', \
       '../runsQ/MIx32x1/q7_v0dt1_i40_prob.txt' u ($1):($2) w lp ls 4 t '12', \
       '../runsQ/MIx32x1/q7_v0dt1_i48_prob.txt' u ($1):($2) w lp ls 5 t '4', \


plot [-1:1][0.4:0.65] \
       '../runsQ/MIx32x1/q7_v0dt1_i08_probphi.txt' u ($1/pi):($2*pi) w l ls 6 t '44', \
0.5 w l ls 0 notitle

#------------------------

set lmargin at screen 0.80
set rmargin at screen 0.96

set ytics 0.1 offset 0.2,-0.2
set xtics 1 offset 0, 0.3

#plot [0:30][1e-6:1] \
       '../runsQ/MIx32x1/q7b_g3500dt5_i16_prob.txt' u ($1):($2) w lp ls 5 t '4', \
       '../runsQ/MIx32x1/q7b_g3500dt5_i24_prob.txt' u ($1):($2) w lp ls 4 t '12', \
       '../runsQ/MIx32x1/q7b_g3500dt5_i32_prob.txt' u ($1):($2) w lp ls 3 t '20', \
       '../runsQ/MIx32x1/q7b_g3500dt5_i40_prob.txt' u ($1):($2) w lp ls 2 t '28', \
       '../runsQ/MIx32x1/q7b_g3500dt5_i48_prob.txt' u ($1):($2) w lp ls 1 t '32', \
       '../runsQ/MIx32x1/q7b_g3500dt5_i56_prob.txt' u ($1):($2) w lp ls 6 t '44', \
exp(-x) w l ls 10 notitle


plot [-1:1][0.35:0.6] \
       '../runsQ/MIx32x1/q7b_g3500dt5_i56_probphi.txt' u ($1/pi):($2*pi) w l ls 6 t '44', \
0.5 w l

#==========================================================================


unset multiplot

set output








