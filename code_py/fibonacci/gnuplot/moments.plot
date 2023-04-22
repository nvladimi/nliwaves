set term pdf fsize 8 enhanced color dashed lw 1 size 8, 10

set datafile commentschars '%#'


set style line 1 lt 1 pt 1 ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 6 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 2 ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4 ps 1    lc rgb "#800080"         lw 2

set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'moments.pdf'

set multiplot
set size 1/2., 1/3.

set logscale y
set grid
set key spacing 1.2

# 1.i  2.Fi   3.nk  4.m2k  5.m3k   6.gk  7.fk  8.nnk  9.2re(ss1k)  10.2re(ss2k) 

set key l

#==========================================================================

set ytics (0.5, 1, 2, 4, 8)

set origin 0, 2/3.

set title 'V = 1,   {/Symbol g}_L = 1.5,   dt = 0.1'
set ylabel '|a_i|^4 / 2 n_i^2'

set xlabel '(i_p + 1) - i'

plot [-10:40][0.5:4.0] \
 '../runsQ/POST5/q5_v0m20_mom.txt'  u (-$1+21):($4)     w lp ls 6  t 'i_p = 20 out of 40', \
 '../runsQ/POST5/q5_v0m30_mom.txt'  u (-$1+31):($4)     w lp ls 1  t 'i_p = 30 out or 40', \
 '../runsQ/POST7/q7_v0dt1_mom.txt'  u (-$1+51):($4)     w lp ls 5  t 'i_p = 50 out of 60', \
 exp( 0.017*x ) w l ls 10


set origin 0.5, 2/3.

set title 'V = Fi,  {/Symbol g}_R = 140,  dt = 10^{-4}  or  {/Symbol g}_R = 3500,  dt = 10^{-5}' 

set xlabel 'i - (i_p + 1)'

plot [][0.5:4] \
 '../runsQ/POST6/q6_m20gr140_mom.txt'    u ($1-21):($4)   w lp ls 6  t 'i_p = 20 out of 40', \
 '../runsQ/POST6/q6_m10gr140_mom.txt'    u ($1-11):($4)   w lp ls 1  t 'i_p = 10 out of 40', \
 '../runsQ/POST7b/q7b_g3500dt5_mom.txt'  u ($1-11):($4)   w lp ls 5  t 'i_p = 10 out of 60', \
exp( 0.015*x ) w l ls 10


#==========================================================================


set origin 0, 1/3.

set title 'V = 1,   {/Symbol g}_L = 1.5,  dt = 0.1'
set ylabel '|a_i|^6 / 6 n_i^3'


set xlabel '(i_p + 1) - i'


plot [-10:40][0.5:10] \
  '../runsQ/POST5/q5_v0m20_mom.txt'  u (-$1+21):($5)     w lp ls 6  t 'i_p = 20 out of 40', \
  '../runsQ/POST5/q5_v0m30_mom.txt'  u (-$1+31):($5)     w lp ls 1  t 'i_p = 30 out of 40', \
  '../runsQ/POST7/q7_v0dt1_mom.txt'  u (-$1+51):($5)     w lp ls 5  t 'i_p = 50 out of 60', \
exp( 0.050*x ) w l ls 10


set origin 0.5, 1/3.

set title 'V = Fi,  {/Symbol g}_R = 140,  dt = 10^{-4}  or  {/Symbol g}_R = 3500,  dt = 10^{-5}' 

set xlabel 'i - (i_p + 1)'

plot [][0.5:8] \
 '../runsQ/POST6/q6_m20gr140_mom.txt'    u ($1-21):($5)   w lp ls 6  t 'i_p = 20 out of 40', \
 '../runsQ/POST6/q6_m10gr140_mom.txt'    u ($1-11):($5)   w lp ls 1  t 'i_p = 10 out of 40', \
 '../runsQ/POST7b/q7b_g3500dt5_mom.txt'  u ($1-11):($5)   w lp ls 5  t 'i_p = 10 out of 60', \
exp( 0.043*x ) w l ls 10


#==========================================================================

set ytics auto

set origin 0, 0.

set title 'V = 1,   {/Symbol g}_L = 1.5,  dt = 0.1'
set ylabel 'g_i'

set xlabel '(i_p + 1) - i'


plot [-10:40][0.01:10] \
  '../runsQ/POST5/q5_v0m20_mom.txt'  u (-$1+21):($6/$8)     w lp ls 6  t 'i_p = 20 out of 40', \
  '../runsQ/POST5/q5_v0m30_mom.txt'  u (-$1+31):($6/$8)     w lp ls 1  t 'i_p = 30 out of 40', \
  '../runsQ/POST7/q7_v0dt1_mom.txt'  u (-$1+51):($6/$8)     w lp ls 5  t 'i_p = 50 out of 60'


set origin 0.5, 0.

set title 'V = Fi,  {/Symbol g}_R = 140,  dt = 10^{-4}  or  {/Symbol g}_R = 3500,  dt = 10^{-5}' 

set xlabel 'i - (i_p + 1)'

plot [][0.01:8] \
 '../runsQ/POST6/q6_m20gr140_mom.txt'   u ($1-21):($6/$8)   w lp ls 6  t 'i_p = 20 out of 40', \
 '../runsQ/POST6/q6_m10gr140_mom.txt'   u ($1-11):($6/$8)   w lp ls 1  t 'i_p = 30 out of 40', \
 '../runsQ/POST7b/q7b_g3500dt5_mom.txt' u ($1-11):($6/$8)   w lp ls 5  t 'i_p = 10 out of 60'


#==========================================================================
unset multiplot

set output



