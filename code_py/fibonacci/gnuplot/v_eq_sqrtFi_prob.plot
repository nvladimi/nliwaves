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



set output 'v_eq_sqrtFi_prob.pdf'

set multiplot
set size 1/2., 1/3.

#set logscale y
set grid

set key spacing 1.2


set xlabel 'n_k / < n_k >' offset 0, 0.5
set ylabel '(probability) exp( n/<n_k> )'
set xrange [0:10]

f(x) = exp(x)

I=2

set tmargin 1.5
set title offset 0,-0.7
set bmargin 3


# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 

# '../runsK/POST1/k_p1gBdt2'
# '../runsK/POST2/k_p2gBdt1'
# '../runsK/POST3/k_p2gAdt2'


# '../runsK/POST2/m_p0gZdt3'
#  '../runsK/POST3/m_p1gAdt3'
#  '../runsK/POST3/m_p2gBdt3'


set key l

set xrange [0:9.7]

#==========================================================================

set origin 0, 2/3.

set title 'P = 0.1,   {/Symbol g} = 3,   dt = 0.01,  40 modes'

plot [][0.5:5]  '../runsK/PROB/k_p1gBdt2_probnk.txt' \
     u  ($1):(f($1)*$2) i I  w lp ls 1  t "i =   5", \
  '' u  ($1):(f($1)*$3) i I  w lp ls 2  t "i = 10", \
  '' u  ($1):(f($1)*$4) i I  w lp ls 3  t "i = 15", \
  '' u  ($1):(f($1)*$5) i I  w lp ls 4  t "i = 26", \
  '' u  ($1):(f($1)*$6) i I  w lp ls 5  t "i = 31", \
  '' u  ($1):(f($1)*$7) i I  w lp ls 6  t "i = 36"


set origin 0.5, 2/3.

set title 'P = 1,   {/Symbol g} = 30,   dt = 0.001,  60 modes'

plot [][0.8:2.4]  '../runsK/PROB/m_p0gZdt3_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "m = 10", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "m = 20",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "m = 25",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 4  t "m = 36",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 5  t "m = 41",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 6  t "m = 51"


#==========================================================================


set origin 0, 1/3.

set title 'P = 0.01,   {/Symbol g} = 3,   dt = 0.1,  40 modes'

plot [][0.5:5]  '../runsK/PROB/k_p2gBdt1_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "i = 10",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 4  t "i = 26",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 5  t "i = 31",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 6  t "i = 36"



set origin 0.5, 1/3.


set title 'P = 0.1,   {/Symbol g} = 10,   dt = 0.001,  60 modes'

set key l 

plot [][0.8:2.4]  '../runsK/PROB/m_p1gAdt3_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "m = 10", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "m = 20",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "m = 25",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 4  t "m = 36",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 5  t "m = 41",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 6  t "m = 51"


#==========================================================================


set origin 0, 0.

set title 'P = 0.01,   {/Symbol g} = 1,   dt = 0.01,  40 modes'

plot [][0.5:5]  '../runsK/PROB/k_p2gAdt2_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "m =   5", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "m = 10",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "m = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 4  t "m = 26",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 5  t "m = 31",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 6  t "m = 36"



set origin 0.5, 0.


set title 'P = 0.01,   {/Symbol g} = 3,    dt = 0.001,  60 modes'

plot [][0.8:2.4]  '../runsK/PROB/m_p2gBdt3_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "m = 10", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "m = 20",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "m = 25",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 4  t "m = 36",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 5  t "m = 41",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 6  t "m = 51"


#==========================================================================



unset multiplot

set output



