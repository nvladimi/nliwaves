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

set output 'prob_comp.pdf'

set multiplot


set xlabel font "Times-Italic"
set ylabel font "Times-Italic"
set key font "Times-Italic"
set title font "Times-Italic"


set key r spacing 1 samplen  2

LM1=0.10;  RM1=0.50
LM2=0.58;  RM2=0.98

BM1=0.06;  TM1=0.46
BM2=0.57;  TM2=0.97

set title offset 0, -0.5
#==========================================================================

set tmargin at screen TM2
set bmargin at screen BM2

set xlabel 'u = ln (x) / k,    where   x = n_i / < n_i >' offset 0, 0
set ylabel 'v = ln [x P (x)] / k' offset -1,0

#set format y "10^{%L}"

#set logscale y


p(k) = k**(-1)
q(k) = k**(-1)


#------------------------
set lmargin at screen LM1
set rmargin at screen RM1

set title "{/Symbol a} = 0,  p = -1,  q = 1"

plot [0:0.2][-0.6:0] \
       '../runsQ/MIx64x05/q7_v0dt1_i08_prob.txt' u (log($1)*p(44)):(log($1*$2)*q(44)) w l ls 6 t '44', \
       '../runsQ/MIx64x05/q7_v0dt1_i16_prob.txt' u (log($1)*p(36)):(log($1*$2)*q(36)) w l ls 1 t '36', \
       '../runsQ/MIx64x05/q7_v0dt1_i24_prob.txt' u (log($1)*p(28)):(log($1*$2)*q(28)) w l ls 2 t '28', \
       '../runsQ/MIx64x05/q7_v0dt1_i32_prob.txt' u (log($1)*p(20)):(log($1*$2)*q(20)) w l ls 3 t '20', \
       '../runsQ/MIx64x05/q8_v0dt1_i28_prob.txt' u (log($1)*p(44)):(log($1*$2)*q(44)) w p ls 6 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i36_prob.txt' u (log($1)*p(36)):(log($1*$2)*q(36)) w p ls 1 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i44_prob.txt' u (log($1)*p(28)):(log($1*$2)*q(28)) w p ls 2 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i52_prob.txt' u (log($1)*p(20)):(log($1*$2)*q(20)) w p ls 3 notitle

#------------------------

set lmargin at screen LM2
set rmargin at screen RM2

set title "{/Symbol a} = 1,  p = -1,  q = 1"

plot [0:0.2][-0.6:0] \
       '../runsQ/MIx64x05/q7b_g3500dt5_i32_prob.txt' u (log($1)*p(20)):(log($1*$2)*q(20)) w l ls 3 t '20', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i40_prob.txt' u (log($1)*p(28)):(log($1*$2)*q(28)) w l ls 2 t '28', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i48_prob.txt' u (log($1)*p(32)):(log($1*$2)*q(32)) w l ls 1 t '32', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i56_prob.txt' u (log($1)*p(44)):(log($1*$2)*q(44)) w l ls 6 t '44'



#==========================================================================

set tmargin at screen TM1
set bmargin at screen BM1

#------------------------

set lmargin at screen LM1
set rmargin at screen RM1

p(k) = k**(-0.8)
q(k) = k**(-1)

set title "{/Symbol a} = 0,  p = -0.8,  q = -1"

set xlabel 'u = ln (x) k^{ p},    where   x = n_i / < n_i >' offset 0, 0
set ylabel 'v = ln [x P (x)] k^{ q}' offset -1,0


plot [][] \
       '../runsQ/MIx64x05/q7_v0dt1_i08_prob.txt' u (log($1)*p(44)):(log($1*$2)*q(44)) w l ls 6 t '44', \
       '../runsQ/MIx64x05/q7_v0dt1_i16_prob.txt' u (log($1)*p(36)):(log($1*$2)*q(36)) w l ls 1 t '36', \
       '../runsQ/MIx64x05/q7_v0dt1_i24_prob.txt' u (log($1)*p(28)):(log($1*$2)*q(28)) w l ls 2 t '28', \
       '../runsQ/MIx64x05/q7_v0dt1_i32_prob.txt' u (log($1)*p(20)):(log($1*$2)*q(20)) w l ls 3 t '20', \
       '../runsQ/MIx64x05/q8_v0dt1_i28_prob.txt' u (log($1)*p(44)):(log($1*$2)*q(44)) w p ls 6 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i36_prob.txt' u (log($1)*p(36)):(log($1*$2)*q(36)) w p ls 1 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i44_prob.txt' u (log($1)*p(28)):(log($1*$2)*q(28)) w p ls 2 notitle, \
       '../runsQ/MIx64x05/q8_v0dt1_i52_prob.txt' u (log($1)*p(20)):(log($1*$2)*q(20)) w p ls 3 notitle


#-------------------------

set lmargin at screen LM2
set rmargin at screen RM2


#set label 2 "i - i_p"

p(k) = k**(-0.9)
q(k) = k**(-1)

set title "{/Symbol a} = 1,  p = -0.9,  q = -1"


plot [][] \
       '../runsQ/MIx64x05/q7b_g3500dt5_i32_prob.txt' u (log($1)*p(20)):(log($1*$2)*q(20)) w l ls 3 t '20', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i40_prob.txt' u (log($1)*p(28)):(log($1*$2)*q(28)) w l ls 2 t '28', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i48_prob.txt' u (log($1)*p(32)):(log($1*$2)*q(32)) w l ls 1 t '32', \
       '../runsQ/MIx64x05/q7b_g3500dt5_i56_prob.txt' u (log($1)*p(44)):(log($1*$2)*q(44)) w l ls 6 t '44'


#==========================================================================

#==========================================================================


unset multiplot

set output








