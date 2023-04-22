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



set output 'v_eq_Fi.pdf'

set multiplot
set size 1/2., 1/3.

#set logscale y
set grid
#unset key

set xlabel 'i'
set ylabel 'n_i F_i'


# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 

set key spacing 1.2
#==========================================================================

set origin 0, 2/3.
set ylabel 'n_i F_i'

plot [][] \
  '../runsQ/POST1/q1_g01o30dt4_spc.txt'   u ($1):($3)   w lp ls 1 t '{/Symbol g}_L = 1    ', \
  '../runsQ/POST4/q4_gr30dt4_spc.txt'     u ($1):($3)   w lp ls 2 t '{/Symbol g}_L = 0.1 ', \
  '../runsQ/POST5/q5_v1L2r030_spc_s1.txt' u ($1):($3)   w lp ls 4 t '{/Symbol g}_L = 10^{-2}', \
  '../runsQ/POST5/q5_v1L3r030_spc_s1.txt' u ($1):($3)   w lp ls 6 t '{/Symbol g}_L = 10^{-3}'


set origin 0.5, 2/3.
set logscale y

set ytics ('1/16' 1/16., '1/8' 1/8., '1/4' 1/4., '1/2' 1/2., 1, 2, 4, 8)

set key spacing 1.5

plot [][1/16.:8] \
  '../runsQ/POST5/q5_v1L0r030_spc_s1.txt'  u ($1):($3)  w lp ls 1  t '{/Symbol g}_R =   30', \
  '../runsQ/POST5/q5_v1L0r060_spc_s1.txt'  u ($1):($3)  w lp ls 2  t '{/Symbol g}_R =   60', \
  '../runsQ/POST6/q6_m20gr140_spc_s1.txt'  u ($1):($3)  w lp ls 3  t '{/Symbol g}_R = 140' , \
  '../runsQ/POST5/q5_v1L0r200_spc_s1.txt'  u ($1):($3)  w lp ls 4  t '{/Symbol g}_R = 200', \
  '../runsQ/POST5/q5_v1L0r300_spc_s1.txt'  u ($1):($3)  w lp ls 6  t '{/Symbol g}_R = 300', \
 '' u ($1):(50*($2)**(-1/3.)) w l ls 10   t '50 F_i^{-1/3}'         
  

unset logscale
set ytics auto

#  '../runsQ/POST6/q6_m20gr140_spc_s1.txt'  u ($1):($3)  w lp ls 7  t '{/Symbol g}_R = 150', \


#==========================================================================

set origin 0, 1/3.
set ylabel '{/Symbol P}'

set key b r spacing 1.2


plot [][-1.5:0.5] \
 '../runsQ/POST1/q1_g01o30dt4_spc.txt'  u ($1):($8)    w lp ls 1 t '{/Symbol g}_L = 1,   ', \
 '../runsQ/POST4/q4_gr30dt4_spc.txt'    u ($1):($8)    w lp ls 2 t '{/Symbol g}_L = 0.1  ', \
 '../runsQ/POST5/q5_v1L2r030_spc_s1.txt' u ($1):($8)   w lp ls 4 t '{/Symbol g}_L = 10^{-2}', \
 '../runsQ/POST5/q5_v1L3r030_spc_s1.txt' u ($1):($8)   w lp ls 6 t '{/Symbol g}_L = 10^{-3}'
 


set origin 0.5, 1/3.

plot [][-0.1:0.1] \
  '../runsQ/POST5/q5_v1L0r030_spc_s1.txt'  u ($1):($8)  w lp ls 1  t '{/Symbol g}_R =   30', \
  '../runsQ/POST5/q5_v1L0r060_spc_s1.txt'  u ($1):($8)  w lp ls 2  t '{/Symbol g}_R =   60', \
  '../runsQ/POST6/q6_m20gr140_spc_s1.txt'  u ($1):($8)  w lp ls 3  t '{/Symbol g}_R = 140' , \
  '../runsQ/POST5/q5_v1L0r200_spc_s1.txt'  u ($1):($8)  w lp ls 4  t '{/Symbol g}_R = 200', \
  '../runsQ/POST5/q5_v1L0r300_spc_s1.txt'  u ($1):($8)  w lp ls 6  t '{/Symbol g}_R = 300'

#  '../runsQ/POST6/q6_m20gr140_spc_s1.txt'  u ($1):($8)  w lp ls 7  t '{/Symbol g}_R = 150', \

 
#==========================================================================


set origin 0, 0
set ylabel '| {/Symbol x} |'

# ksi = Pi (V F)^{-1} n^{-3/2}

set ytics 0.1

p=-1.5
q=-2.0

set key t r

plot [][0:0.4] \
 '../runsQ/POST1/q1_g01o30dt4_spc.txt'   u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 1 t '{/Symbol g}_L = 1    ', \
 '../runsQ/POST4/q4_gr30dt4_spc.txt'     u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 2 t '{/Symbol g}_L = 0.1  ', \
 '../runsQ/POST5/q5_v1L2r030_spc_s1.txt' u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 4 t '{/Symbol g}_L = 10^{-2}', \
 '../runsQ/POST5/q5_v1L3r030_spc_s1.txt' u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 6 t '{/Symbol g}_L = 10^{-3}'



set origin 0.5, 0
set key l 

plot [][0:0.4] \
  '../runsQ/POST5/q5_v1L0r030_spc_s1.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 1  t '{/Symbol g}_R =   30', \
  '../runsQ/POST5/q5_v1L0r060_spc_s1.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 2  t '{/Symbol g}_R =   60', \
  '../runsQ/POST6/q6_m20gr140_spc_s1.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 3  t '{/Symbol g}_R = 140' , \
  '../runsQ/POST5/q5_v1L0r200_spc_s1.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 4  t '{/Symbol g}_R = 200', \
  '../runsQ/POST5/q5_v1L0r300_spc_s1.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 6  t '{/Symbol g}_R = 300'

#  '../runsQ/POST6/q6_m20gr140_spc_s1.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 7  t '{/Symbol g}_R = 150', \





#==========================================================================


unset multiplot

set output



