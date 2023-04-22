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



set output 'spectra.pdf'

set multiplot


#set logscale y
set grid
unset key
set format y "10^{%L}"

set rmargin 0
set tmargin 0
set bmargin 0

#==========================================================================

set lmargin at screen 0.06
set rmargin at screen 0.49
set tmargin at screen 0.96
set bmargin at screen 0.68


set title 'V = 1,   {/Symbol g}_L = 1.5,  {/Symbol g}_R = 0,    dt = 0.1'
set ylabel 'n_i F_i'

set logscale y
set xlabel 'i'
set key l spacing 1.75

set ytics 10

plot [][] \
 '../runsQ/POST5/q5_v0m20_spc.txt'  u ($1):($3)   w lp ls 6  t 'i_p = 20 out of 40', \
 '../runsQ/POST5/q5_v0m30_spc.txt'  u ($1):($3)   w lp ls 1  t 'i_p = 30 out of 40', \
 '../runsQ/POST7/q7_v0dt1_spc.txt'  u ($1):($3)   w lp ls 5  t 'i_p = 50 out of 60', \
 '' u ($1):(22*($2)**(1/3.)) w l ls 10 t '22 F_i^{1/3}' 



set lmargin at screen 0.56
set rmargin at screen 0.99

set key r

#set ytics ('1/8' 1/8.,  '1/2' 1/2., 2, 8, 32, 128)

set title 'V = Fi,  {/Symbol g}_L = 0,  {/Symbol g}_R = 140, dt = 10^{-4} or {/Symbol g}_R = 3500,  dt = 10^{-5}'

plot [][] \
 '../runsQ/POST6/q6_m20gr140_spc.txt'    u ($1):($3)   w lp ls 6  t 'i_p = 20 out of 40', \
 '../runsQ/POST6/q6_m10gr140_spc.txt'    u ($1):($3)   w lp ls 1  t 'i_p = 10 out of 40', \
 '../runsQ/POST7b/q7b_g3500dt5_spc.txt'  u ($1):($3)   w lp ls 5  t 'i_p = 10 out of 60', \
 '' u ($1):(50*($2)**(-1/3.)) w l ls 10   t '50 F_i^{-1/3}'         


set ytics auto

unset title

#==========================================================================

unset logscale
set format y "%g"


set lmargin at screen 0.06
set rmargin at screen 0.49
set tmargin at screen 0.64
set bmargin at screen 0.36

set key spacing 1.25

set ylabel '| {/Symbol x} |'

p=-1.5
q=-1.0

#set xlabel 'i - (i_p + 1)'

plot [][] \
  '../runsQ/POST5/q5_v0m20_spc.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 6 t 'i_p = 20 out of 40', \
  '../runsQ/POST5/q5_v0m30_spc.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 1 t 'i_p = 30 out of 40', \
  '../runsQ/POST7/q7_v0dt1_spc.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 5  t 'i_p = 50 out of 60'


set origin 0.5, 1/3.
set lmargin at screen 0.56
set rmargin at screen 0.99



p=-1.5
q=-2.0

plot [][0:0.3] \
 '../runsQ/POST6/q6_m20gr140_spc.txt'    u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 6  t 'i_p = 20 out of 40', \
 '../runsQ/POST6/q6_m10gr140_spc.txt'    u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 1  t 'i_p = 10 out of 40', \
 '../runsQ/POST7b/q7b_g3500dt5_spc.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 5  t 'i_p = 10 out of 60'


#==========================================================================

set origin 0, 0.



#==========================================================================
unset multiplot

set output



