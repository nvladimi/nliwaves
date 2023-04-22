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



set output 'v_eq_sqrtFi_lever_prob.pdf'

set multiplot
set size 1/2., 1/3.

set grid
unset key

set tmargin 1.5
set title offset 0,-0.7
set bmargin 3

set xlabel 'n_k / < n_k >' offset 0, 0.5
set ylabel '(probability) exp( n/<n_k> )'
set xrange [0:9.7]
set ytics 1

f(x) = exp(x)

I=2

set key l

set yrange [0.5:4.5]

#==========================================================================

set origin 0, 2/3.

set title 'pumping at mode 3'

set label "in all runs:   {/Symbol g} = 1,   P = 0.01 Fi( 20 ) / Fi( pumping ),  dt = 10^{-3}" at 5,0.66 c

plot [][]  '../runsL/PROB/l1_typ303_probnk.txt' \
       i I    u  ($1):(f($1)*$3)  w lp ls 1  t "i = 10", \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 2  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 3  t "i = 20",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 4  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 5  t "i = 30",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 6  t "i = 35"

unset label

#--------------------------------


set origin 1/2., 2/3.

set title 'pumping at mode 38'

plot [][]  '../runsL/PROB/l1_typ338_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "i = 10",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 4  t "i = 20",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 5  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 6  t "i = 30"




#==========================================================================

set origin 0, 1/3.

set title 'pumping at mode 5'

plot [][]  '../runsL/PROB/l1_typ305_probnk.txt' \
       i I    u  ($1):(f($1)*$3)  w lp ls 1  t "i = 10", \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 2  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 3  t "i = 20",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 4  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 5  t "i = 30",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 6  t "i = 35"



#--------------------------------

set origin 1/2., 1/3.

set title 'pumping at mode 36'

plot [][]  '../runsL/PROB/l1_typ336_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "i = 10",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 4  t "i = 20",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 5  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 6  t "i = 30"

#==========================================================================

set origin 0, 0.

set title 'pumping at mode 10'

plot [][]  '../runsL/PROB/l1_typ310_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 2  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 3  t "i = 20",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 4  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 5  t "i = 30",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 6  t "i = 35"


#--------------------------------


set origin 1/2., 0.

set title 'pumping at mode 30'


plot [][]  '../runsL/PROB/l1_typ330_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "i = 10",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 4  t "i = 20",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 5  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 6  t "i = 35"


#==========================================================================
unset multiplot

set output








