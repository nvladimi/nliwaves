set term pdf fsize 7 enhanced color dashed lw 1 size 4, 2

set datafile commentschars '%#'


set style line 1 lt 1 pt 4 ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 6 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 2 ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1.3   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4 ps 1    lc rgb "#800080"         lw 2


set style line 11 lt 3 dt 3 pt 1 ps 1    lc rgb "#D00040"         lw 3
set style line 22 lt 3 dt 3 pt 6 ps 1    lc rgb "#FF4000"         lw 3
set style line 33 lt 3 dt 3 pt 2 ps 1    lc rgb "orange"          lw 3
set style line 44 lt 3 dt 3 pt 4 ps 1    lc rgb "#60B000"         lw 3
set style line 55 lt 3 dt 3 pt 6 ps 1    lc rgb "#008080"         lw 3
set style line 66 lt 3 dt 3 pt 12 ps 1   lc rgb  "blue"           lw 3
set style line 77 lt 3 dt 3 pt 4 ps 1    lc rgb "#800080"         lw 3

set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'fig_sqrt_moments.pdf'

set multiplot
set size 1/2., 1



set xlabel font "Times-Italic"
set ylabel font "Times-Italic"
set key font "Times-Italic"

set xlabel 'i - p'

# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 
# 1.i  2.Fi   3.nk  4.m2k  5.m3k   5.gk  6.fk  8.nnk  9.2re(ss1k)  10.2re(ss2k) 


ia=20
ib=30

LM1=0.09;  RM1=0.48
LM2=0.59;  RM2=0.98

BM1=0.16;  TM1=0.96

set tmargin at screen TM1
set bmargin at screen BM1

#==========================================================================

set ylabel '{/Symbol |} a_i {/Symbol |}^4 / 2n_i^2' offset 0.5,0

set lmargin at screen LM1
set rmargin at screen RM1

set key l samplen 2 # at -35,1.485
set xlabel offset 0,0.5

set ytics 0.1
plot [][1:1.5] \
  '../runsK/POST1/k_p1gBdt2_mom.txt'  u ($1-ia):($4)   w lp ls 1 t  "40 modes", \
  '../runsK/POST2/m_p0gZdt3_mom.txt'  u ($1-ib):($4)   w lp ls 6 t  "60 modes"

#  '../runsK/POST1/k_p1gBdt2_mom.txt'  u ($1-ia):($4)   w lp ls 1 t  "{/Symbol g} =   3,  P = 0.1  ", \
#  '../runsK/POST2/m_p0gZdt3_mom.txt'  u ($1-ib):($4)   w lp ls 6 t  "{/Symbol g} = 30,  P = 1     "
# '../runsK/POST3/m_p1gAdt3_mom.txt'  u ($1-ib):($4)   w lp ls 4 t  "{/Symbol g} = 10,  P = 0.1  ", \
#  '../runsK/POST3/k_p2gAdt2_mom.txt'  u ($1-ia):($4)   w lp ls 3 t  "{/Symbol g} =   1,  P = 0.01"


#==========================================================================

set ylabel '{/Symbol |} a_i {/Symbol |}^6 / 6n_i^4' offset 0.5,0

set lmargin at screen LM2
set rmargin at screen RM2

unset key

plot [][1:1.5] \
  '../runsK/POST1/k_p1gBdt2_mom.txt'  u ($1-ia):($5)   w lp ls 1 t  "{/Symbol g} = 3,  P = 0.1  ", \
  '../runsK/POST2/m_p0gZdt3_mom.txt'  u ($1-ib):($5)   w lp ls 6 t  "{/Symbol g} = 30,  P = 1     "

#  '../runsK/POST3/m_p1gAdt3_mom.txt'  u ($1-ib):($5)   w lp ls 4 t  "{/Symbol g} = 10,  P = 0.1  "
# '../runsK/POST3/k_p2gAdt2_mom.txt'  u ($1-ia):($5)   w lp ls 3 t  "{/Symbol g} = 1,  P = 0.01", \

#==========================================================================


unset multiplot

set output



