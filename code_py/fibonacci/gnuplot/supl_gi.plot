set term pdf fsize 8 enhanced color dashed lw 1 size 8, 2.8

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



set output 'supl_gi.pdf'

set multiplot
set size 1/3., 1

set logscale y
#set grid
set key spacing 1.2

# 1.i  2.Fi   3.nk  4.m2k  5.m3k   6.gk  7.fk  8.nnk  9.2re(ss1k)  10.2re(ss2k) 

set key b samplen 2

#==========================================================================

set lmargin 7
set bmargin 4


#==========================================================================

set ytics auto

set origin 0, 0.

#set title 'V = 1,   {/Symbol g}_L = 1.5,  dt = 0.1'
set ylabel 'g_i' offset 1,0

set xlabel '(i_p + 1) - i'

set label 1 'V = 1' at -15, 6

plot [-20:50][0.01:10] \
  '../runsQ/POST5/q5_v0m20_mom.txt'  u (-$1+21):($6/$8)     w lp ls 6  t 'i_p = 20 of 40', \
  '../runsQ/POST5/q5_v0m30_mom.txt'  u (-$1+31):($6/$8)     w lp ls 1  t 'i_p = 30 of 40', \
  '../runsQ/POST7/q7_v0dt1_mom.txt'  u (-$1+51):($6/$8)     w lp ls 5  t 'i_p = 50 of 60'


set origin 1/3., 0.

#set title 'V = Fi,  {/Symbol g}_R = 140,  dt = 10^{-4}  or  {/Symbol g}_R = 3500,  dt = 10^{-5}' 

set label 1 'V = Fi' at -15, 5

set xlabel 'i - (i_p + 1)'

plot [][0.01:8] \
 '../runsQ/POST6/q6_m20gr140_mom.txt'   u ($1-21):($6/$8)   w lp ls 6  t 'i_p = 20 of 40', \
 '../runsQ/POST6/q6_m10gr140_mom.txt'   u ($1-11):($6/$8)   w lp ls 1  t 'i_p = 30 of 40', \
 '../runsQ/POST7b/q7b_g3500dt5_mom.txt' u ($1-11):($6/$8)   w lp ls 5  t 'i_p = 10 of 60'



set label 1 'V = Fi^{1/2}' at -10, 0.47

ia=20
ib=30
unset logscale 

set origin 2/3., 0.
set ylabel 'g_i'

set xlabel 'i - i_p'

set key t c at -5, 0.43

plot [][0:0.5] \
  '../runsK/POST1/k_p1gBdt2_mom.txt'  u ($1-ia):($6/$8)   w lp ls 1 t  "{/Symbol g} =   3,  P = 0.1  ", \
  '../runsK/POST3/k_p2gAdt2_mom.txt'  u ($1-ia):($6/$8)   w lp ls 3 t  "{/Symbol g} =   1,  P = 0.01", \
  '../runsK/POST2/m_p0gZdt3_mom.txt'  u ($1-ib):($6/$8)   w lp ls 4 t  "{/Symbol g} = 30,  P = 1     ", \
  '../runsK/POST3/m_p1gAdt3_mom.txt'  u ($1-ib):($6/$8)   w lp ls 6 t  "{/Symbol g} = 10,  P = 0.1  "

#==========================================================================


#==========================================================================
unset multiplot

set output



