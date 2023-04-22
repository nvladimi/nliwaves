set term pdf fsize 8 enhanced color dashed lw 1 size 8, 9

set datafile commentschars '%#'


set style line 1 lt 1 pt 1 ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 6 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 2 ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1   lc rgb  "blue"           lw 2
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



set output 'v_eq_sqrtFi_moments.pdf'

set multiplot
set size 1/2., 1/3.

#set logscale y
set grid

set key spacing 1.2

set xlabel 'i'
set ylabel 'n_i F_i'


# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 
# 1.i  2.Fi   3.nk  4.m2k  5.m3k   5.gk  6.fk  8.nnk  9.2re(ss1k)  10.2re(ss2k) 



#==========================================================================

set origin 0, 2/3.
set ylabel '|a_i|^4 / 2 n_i^2'

set key

plot [][1:1.2] \
  '../runsK/POST1/k_p1gBdt2_mom.txt'  u ($1):($4)   w lp ls 6 t  "{/Symbol g} = 3,  P = 0.1  ", \
  '../runsK/POST2/k_p2gBdt1_mom.txt'  u ($1):($4)   w lp ls 4 t  "{/Symbol g} = 3,  P = 0.01",  \
  '../runsK/POST3/k_p2gAdt2_mom.txt'  u ($1):($4)   w lp ls 1 t  "{/Symbol g} = 1,  P = 0.01"

set origin 0.5, 2/3.

plot [][1:1.2] \
  '../runsK/POST2/m_p0gZdt3_mom.txt'  u ($1):($4)   w lp ls 6 t  "{/Symbol g} = 30,  P = 1     ", \
  '../runsK/POST3/m_p1gAdt3_mom.txt'  u ($1):($4)   w lp ls 1 t  "{/Symbol g} = 10,  P = 0.1  ", \
  '../runsK/POST3/m_p2gBdt3_mom.txt'  u ($1):($4)   w lp ls 4 t  "{/Symbol g} =   3,  P = 0.01"


#==========================================================================


set origin 0, 1/3.
set ylabel '|a_i|^6 / 6 n_i^4'

set key

plot [][1:1.5] \
  '../runsK/POST1/k_p1gBdt2_mom.txt'  u ($1):($5)   w lp ls 6 t  "{/Symbol g} = 3,  P = 0.1  ", \
  '../runsK/POST2/k_p2gBdt1_mom.txt'  u ($1):($5)   w lp ls 4 t  "{/Symbol g} = 3,  P = 0.01",  \
  '../runsK/POST3/k_p2gAdt2_mom.txt'  u ($1):($5)   w lp ls 1 t  "{/Symbol g} = 1,  P = 0.01"

set origin 0.5, 1/3.

plot [][1:1.5] \
  '../runsK/POST2/m_p0gZdt3_mom.txt'  u ($1):($5)   w lp ls 6 t  "{/Symbol g} = 30,  P = 1     ", \
  '../runsK/POST3/m_p1gAdt3_mom.txt'  u ($1):($5)   w lp ls 1 t  "{/Symbol g} = 10,  P = 0.1  ", \
  '../runsK/POST3/m_p2gBdt3_mom.txt'  u ($1):($5)   w lp ls 4 t  "{/Symbol g} =   3,  P = 0.01"


#==========================================================================


set origin 0, 0.
set ylabel 'g_i'

set key

plot [][0:0.5] \
  '../runsK/POST1/k_p1gBdt2_mom.txt'  u ($1):($6/$8)   w lp ls 6 t  "{/Symbol g} = 3,  P = 0.1  ", \
  '../runsK/POST2/k_p2gBdt1_mom.txt'  u ($1):($6/$8)   w lp ls 4 t  "{/Symbol g} = 3,  P = 0.01",  \
  '../runsK/POST3/k_p2gAdt2_mom.txt'  u ($1):($6/$8)   w lp ls 1 t  "{/Symbol g} = 1,  P = 0.01"

set origin 0.5, 0.

plot [][0:0.5] \
  '../runsK/POST2/m_p0gZdt3_mom.txt'  u ($1):($6/$8)   w lp ls 6 t  "{/Symbol g} = 30,  P = 1     ", \
  '../runsK/POST3/m_p1gAdt3_mom.txt'  u ($1):($6/$8)   w lp ls 1 t  "{/Symbol g} = 10,  P = 0.1  ", \
  '../runsK/POST3/m_p2gBdt3_mom.txt'  u ($1):($6/$8)   w lp ls 4 t  "{/Symbol g} =   3,  P = 0.01"


#==========================================================================



unset multiplot

set output



