set term pdf fsize 8 enhanced color dashed lw 1 size 8, 10

set datafile commentschars '%#'


set style line 1 lt 1 pt 8 ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 6 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 10 ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4 ps 1    lc rgb "#800080"         lw 2

set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 4
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'probnk.pdf'

set multiplot
set size 1/2., 1/3.

set grid
unset key

set tmargin 1.5
set title offset 0,-0.7
set bmargin 3

set logscale y

set xlabel 'n_k / < n_k >' offset 0, 0.5
set ylabel 'probability'
set xrange [0:35]
set ytics 10
set format y "10^{%L}"

#f(x) = exp(x)

f(x) = 1

I=2

set key r

set yrange [1e-6:1]

a=100; q=-4
#set xtics 2

#==========================================================================

set origin 0, 2/3.

set title 'V = 1,   {/Symbol g}_L = 1.5,    {/Symbol g}_R = 0,    i_p = 20'

plot [][]  '../runsQ/PROB/q5_v0m20_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "i = 10",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 4  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 5  t "i = 30",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 6  t "i = 35",  \
exp(-x) w l ls 10

unset label

#--------------------------------


set origin 1/2., 2/3.

set title 'V = 1,    {/Symbol g}_L = 1.5,   {/Symbol g}_R = 0,    i_p = 30'

f20(x) = exp(-x) + 0.04*exp(-0.45*x);

plot [][]  '../runsQ/PROB/q5_v0m30_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 7  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 7  t "i = 10",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 1  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 2  t "i = 20",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 3  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 6  t "i = 35",  \
exp(-x) w l ls 10

#  f20(x) w l ls 2 notitle

#==========================================================================

set origin 0, 1/3.

set title 'V = F_i,   {/Symbol g}_L = 0,    {/Symbol g}_R = 100,   i_p = 20'


#plot [][]  '../runsQ/PROB/q5_v1m20_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "i = 10",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 4  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 5  t "i = 30",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 6  t "i = 35", \
exp(-x) w l ls 10



#--------------------------------

set origin 1/2., 1/3.

set title 'V = F_i,   {/Symbol g}_L = 0,    {/Symbol g}_R = 100,   i_p = 10'

#plot [][]  '../runsQ/PROB/q5_v1m10_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 4  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 5  t "i = 20",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 6  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 7  t "i = 30",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 7  t "i = 35", \
exp(-x) w l ls 10

#==========================================================================

set origin 0, 1/3.


set title 'V = F_i,   {/Symbol g}_L = 0,    {/Symbol g}_R = 140,   i_p = 20'

plot [][]  '../runsQ/PROB/q6_m20gr140_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$3)  w lp ls 2  t "i = 10",  \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 3  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 4  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 5  t "i = 30",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 6  t "i = 35",  \
exp(-x) w l ls 10


#--------------------------------


set origin 1/2., 1/3.

set title 'V = F_i,   {/Symbol g}_L = 0,    {/Symbol g}_R = 140,   i_p = 10'

plot [][]  '../runsQ/PROB/q6_m10gr140_probnk.txt' \
       i I    u  ($1):(f($1)*$2)  w lp ls 1  t "i =   5", \
  ''   i I    u  ($1):(f($1)*$4)  w lp ls 4  t "i = 15",  \
  ''   i I    u  ($1):(f($1)*$5)  w lp ls 5  t "i = 20",  \
  ''   i I    u  ($1):(f($1)*$6)  w lp ls 6  t "i = 25",  \
  ''   i I    u  ($1):(f($1)*$7)  w lp ls 7  t "i = 30",  \
  ''   i I    u  ($1):(f($1)*$8)  w lp ls 7  t "i = 35", \
exp(-x) w l ls 10


#==========================================================================
unset multiplot

set output








