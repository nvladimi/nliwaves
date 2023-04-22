set term pdf fsize 7 enhanced color dashed lw 1 size 6, 6

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

set output 'prob_alt.pdf'

set multiplot


set xlabel font "Times-Italic"
set key font "Times-Italic"


unset key

set logscale y
set key r spacing 1 samplen  2

#LM1=0.10;  RM1=0.50
#LM2=0.58;  RM2=0.98
#BM1=0.08;  TM1=0.48
#BM2=0.58;  TM2=0.98

LM1=0.06;  RM1=0.48
LM2=0.56;  RM2=0.98

BM1=0.06;  TM1=0.48
BM2=0.56;  TM2=0.98


#==========================================================================

set tmargin at screen TM2
set bmargin at screen BM2

set xlabel 'n_i / < n_i >' offset 0, 0.5
set ylabel 'probability' offset -1,0

set format y "10^{%L}"

set logscale x

#set ytics 100
set xtics 0.5,2

#------------------------
set lmargin at screen LM1
set rmargin at screen RM1

set label 1 "{/Symbol a = 0}" at 40,1.6e-2
set label 2 "i_p - i" at 4, 2e-5  font "Times-Italic"
set key b l 

set label 3 "lines  - 60 modes" at 20, 8e-3
set label 4 "points - 80 modes" at 20, 4e-3

plot [3:64][1e-7:0.05] \
       '../runsQ/MIx64x05/q7_v0dt1_i08_prob.txt' u ($1):($2) w l ls 6 t '44', \
       '../runsQ/MIx64x05/q7_v0dt1_i16_prob.txt' u ($1):($2) w l ls 1 t '36', \
       '../runsQ/MIx64x05/q7_v0dt1_i24_prob.txt' u ($1):($2) w l ls 2 t '28', \
       '../runsQ/MIx64x05/q7_v0dt1_i32_prob.txt' u ($1):($2) w l ls 3 t '20', \
       '../runsQ/MIx64x05/q7_v0dt1_i40_prob.txt' u ($1):($2) w l ls 4 t '12', \
       '../runsQ/MIx64x05/q7_v0dt1_i48_prob.txt' u ($1):($2) w l ls 5 t '4', \
       '../runsQ/MIx64x05/q8_v0dt1_i28_prob.txt' u ($1):($2) w p ls 6 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i36_prob.txt' u ($1):($2) w p ls 1 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i44_prob.txt' u ($1):($2) w p ls 2 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i52_prob.txt' u ($1):($2) w p ls 3 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i60_prob.txt' u ($1):($2) w p ls 4 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i68_prob.txt' u ($1):($2) w p ls 5 notitle, \
exp(-x) w l ls 10 t 'exp (-x)', \
30*x**(-4) w l ls 10 t 'p = -4'


unset label 3
unset label 4


#------------------------

set lmargin at screen LM2
set rmargin at screen RM2

set label 1 "{/Symbol a = 1}"
unset ylabel
set key b l

#set label 2 "i - i_p"


plot [3:64][1e-7:0.05] \
       '../runsQ/MIx64x05/q7b_g3500dt5_i16_prob.txt' u ($1):($2) w l ls 5 t '4', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i24_prob.txt' u ($1):($2) w l ls 4 t '12', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i32_prob.txt' u ($1):($2) w l ls 3 t '20', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i40_prob.txt' u ($1):($2) w l ls 2 t '28', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i48_prob.txt' u ($1):($2) w l ls 1 t '32', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i56_prob.txt' u ($1):($2) w l ls 6 t '44', \
exp(-x) w l ls 10 t 'exp (-x)', \
30*x**(-4) w l ls 10 t 'p = -4'




#==========================================================================

set tmargin at screen TM1
set bmargin at screen BM1

set logscale
set key r
unset label
set format y "%g"

set ylabel '- log10 ( probability )' offset -1,0

set ytics 0.5,2
set ytics (1,2,3,4,5,6,7,8)


#------------------------

set lmargin at screen LM1
set rmargin at screen RM1

q=0.4
a=1
set label "P = A exp ( - c n^q )" at 3.5, 6.3

plot [3:64][1.5:7] \
       '../runsQ/MIx64x05/q7_v0dt1_i08_prob.txt' u ($1):(-log10($2/a)) w l ls 6 t '44', \
       '../runsQ/MIx64x05/q7_v0dt1_i16_prob.txt' u ($1):(-log10($2/a)) w l ls 1 t '36', \
       '../runsQ/MIx64x05/q7_v0dt1_i24_prob.txt' u ($1):(-log10($2/a)) w l ls 2 t '28', \
       '../runsQ/MIx64x05/q7_v0dt1_i32_prob.txt' u ($1):(-log10($2/a)) w l ls 3 t '20', \
       '../runsQ/MIx64x05/q7_v0dt1_i40_prob.txt' u ($1):(-log10($2/a)) w l ls 4 t '12', \
       '../runsQ/MIx64x05/q7_v0dt1_i48_prob.txt' u ($1):(-log10($2/a)) w l ls 5 t '4', \
       '../runsQ/MIx64x05/q8_v0dt1_i28_prob.txt' u ($1):(-log10($2/a)) w p ls 6 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i36_prob.txt' u ($1):(-log10($2/a)) w p ls 1 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i44_prob.txt' u ($1):(-log10($2/a)) w p ls 2 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i52_prob.txt' u ($1):(-log10($2/a)) w p ls 3 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i60_prob.txt' u ($1):(-log10($2/a)) w p ls 4 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i68_prob.txt' u ($1):(-log10($2/a)) w p ls 5 notitle, \
0.4*x w l ls 10 t 'q = 1', 1.1*x**q w l ls 10 t "q = 0.4"


#-------------------------

set lmargin at screen LM2
set rmargin at screen RM2

plot [3:64][1.5:7] \
       '../runsQ/MIx64x05/q7b_g3500dt5_i16_prob.txt' u ($1):(-log10($2/a)) w l ls 5 t '4', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i24_prob.txt' u ($1):(-log10($2/a)) w l ls 4 t '12', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i32_prob.txt' u ($1):(-log10($2/a)) w l ls 3 t '20', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i40_prob.txt' u ($1):(-log10($2/a)) w l ls 2 t '28', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i48_prob.txt' u ($1):(-log10($2/a)) w l ls 1 t '32', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i56_prob.txt' u ($1):(-log10($2/a)) w l ls 6 t '44', \
0.4*x w l ls 10 t "q = 1", 1.1*x**q w l ls 10 t "q=0.4"

#


#==========================================================================

#==========================================================================


unset multiplot

set output








