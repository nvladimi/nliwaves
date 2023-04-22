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

#set grid
unset key

set logscale y

set xlabel 'n_k / < n_k >' offset 0, 0.3
set ylabel 'probability' offset 0.5,0
set xrange [0:35]
set ytics 10
set format y "10^{%L}"

#f(x) = exp(x)

f(x) = 1

I=2

set key r spacing 1 samplen  2

set yrange [1e-6:1]

a=100; q=-4
#set xtics 2


LM1=0.08;  RM1=0.48
LM2=0.58;  RM2=0.98

BM1=0.08;  TM1=0.48
BM2=0.58;  TM2=0.98


#==========================================================================

set tmargin at screen TM2
set bmargin at screen BM2

set lmargin at screen LM1
set rmargin at screen RM1

set ytics 100

#set title 'V = 1,   {/Symbol g}_L = 1.5,   i_p = 20   or    i_p = 30'

plot [][]  '../runsQ/PROB/q5_v0m20_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 11  notitle,  \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 22  notitle,  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 33  notitle,  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 44  notitle,  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 44  notitle,  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 44  notitle,  \
  '../runsQ/PROB/q5_v0m30_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 6  t "i_p - i = 25", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 7  t "i_p - i = 20",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 1  t "i_p - i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 2  t "i_p - i = 10",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 3  t "i_p - i =   5",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 4  t "i_p - i <   0",  \
exp(-x) w l ls 10

#  f20(x) w l ls 2 notitle

#==========================================================================

set lmargin at screen LM2
set rmargin at screen RM2

#set title 'V = F_i,   {/Symbol g}_R = 140,   i_p = 20   or i_p = 10'

plot [][]  '../runsQ/PROB/q6_m20gr140_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 44  notitle,  \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 44  notitle,  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 44  notitle,  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 33  notitle,  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 22  notitle,  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 11  notitle,  \
            '../runsQ/PROB/q6_m10gr140_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 4  t "i - i_p <   0", \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "i - i_p =   5",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 2  t "i - i_p = 10",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 1  t "i - i_p = 15",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 7  t "i - i_p = 20",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 6  t "i - i_p = 25", \
exp(-x) w l ls 10


#==========================================================================


set tmargin at screen TM1
set bmargin at screen BM1

set lmargin at screen LM1
set rmargin at screen RM1


set ytics (0.5, 1, 2, 4, 8)
set format y "%g"
set key b

#set title 'V = 1,   {/Symbol g}_L = 1.5,   dt = 0.1'
set ylabel '|a_i|^4 / 2 n_i^2,   |a_i|^6 / 6 n_i^3' offset 1,0

set xlabel '(i_p + 1) - i'

plot [-20:50][0.5:10.0] \
 '../runsQ/POST5/q5_v0m20_mom.txt'  u (-$1+21):($4)     w lp ls 6  t 'i_p = 20 of 40', \
 '../runsQ/POST5/q5_v0m30_mom.txt'  u (-$1+31):($4)     w lp ls 1  t 'i_p = 30 of 40', \
 '../runsQ/POST7/q7_v0dt1_mom.txt'  u (-$1+51):($4)     w lp ls 5  t 'i_p = 50 of 60', \
 '../runsQ/POST5/q5_v0m20_mom.txt'  u (-$1+21):($5)     w lp ls 6  notitle, \
 '../runsQ/POST5/q5_v0m30_mom.txt'  u (-$1+31):($5)     w lp ls 1  notitle, \
 '../runsQ/POST7/q7_v0dt1_mom.txt'  u (-$1+51):($5)     w lp ls 5  notitle, \
exp( 0.017*x ) w l ls 10 notitle, \
exp( 0.050*x ) w l ls 10 notitle


#==========================================================================

set lmargin at screen LM2
set rmargin at screen RM2


#set title 'V = Fi,  {/Symbol g}_R = 140,  dt = 10^{-4}  or  {/Symbol g}_R = 3500,  dt = 10^{-5}' 

set xlabel 'i - (i_p + 1)'

plot [-20:50][0.5:10] \
 '../runsQ/POST6/q6_m20gr140_mom.txt'    u ($1-21):($4)   w lp ls 6  t 'i_p = 20 of 40', \
 '../runsQ/POST6/q6_m10gr140_mom.txt'    u ($1-11):($4)   w lp ls 1  t 'i_p = 10 of 40', \
 '../runsQ/POST7b/q7b_g3500dt5_mom.txt'  u ($1-11):($4)   w lp ls 5  t 'i_p = 10 of 60', \
 '../runsQ/POST6/q6_m20gr140_mom.txt'    u ($1-21):($5)   w lp ls 6  notitle, \
 '../runsQ/POST6/q6_m10gr140_mom.txt'    u ($1-11):($5)   w lp ls 1  notitle, \
 '../runsQ/POST7b/q7b_g3500dt5_mom.txt'  u ($1-11):($5)   w lp ls 5  notitle, \
exp( 0.015*x ) w l ls 10 notitle, \
exp( 0.043*x ) w l ls 10 notitle

#==========================================================================


unset multiplot

set output








