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



set output 'v_eq_sqrtFi_lever.pdf'

set multiplot
set size 1/2., 1/3.

#set logscale y
set grid

set xlabel 'i'

unset key

# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 



#==========================================================================

set origin 0, 2/3.
set ylabel 'n_i F_i  {/Symbol P}_@p^{-2/3}'

#set key

p2=(0.01*6765)**(-2/3.)

# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 


plot [][0:4] \
  '../runsL/POST1/l1_typ305_spc.txt'   u ($1):($3*p2)   w lp ls 1 t  "i_p =   5", \
  '../runsL/POST1/l1_typ310_spc.txt'   u ($1):($3*p2)   w lp ls 2 t  "i_p = 10", \
  '../runsK/POST3/k_p2gAdt2_spc.txt'   u ($1):($3*p2)   w lp ls 3 t  "i_p = 20", \
  '../runsL/POST1/l1_typ330_spc.txt'   u ($1):($3*p2)   w lp ls 4 t  "i_p = 30", \
  '../runsL/POST1/l1_typ336_spc.txt'   u ($1):($3*p2)   w lp ls 6 t  "i_p = 36"

#  '../runsL/POST1/l1_typ303_spc.txt'   u ($1):($3*p2)   w lp ls 1 t  "i_p =   3"
#  '../runsL/POST1/l1_typ338_spc.txt'   u ($1):($3*p2)   w lp ls 7 t  "i_p = 38"
#  '../runsK/POST3/k_p2gAdt2_spc.txt'   u ($1):($3*p2)   w lp ls 7 t  "i_p = 20"

#==========================================================================

set origin 0, 1/3.
set ylabel '{/Symbol P} / {/Symbol P}_p '
#set yrange [-2:0]


p2=(0.01*6765)

plot [][-1:1] \
  '../runsL/POST1/l1_typ305_spc.txt'   u ($1):($8/p2)   w lp ls 1 t  "i_p =   5", \
  '../runsL/POST1/l1_typ310_spc.txt'   u ($1):($8/p2)   w lp ls 2 t  "i_p = 10", \
  '../runsK/POST3/k_p2gAdt2_spc.txt'   u ($1):($8/p2)   w lp ls 3 t  "i_p = 20", \
  '../runsL/POST1/l1_typ330_spc.txt'   u ($1):($8/p2)   w lp ls 4 t  "i_p = 30", \
  '../runsL/POST1/l1_typ336_spc.txt'   u ($1):($8/p2)   w lp ls 6 t  "i_p = 36"
 
#==========================================================================

set origin 0, 0.

set ylabel '| {/Symbol x} |'

# ksi = Pi (V F)^{-1} n^{-3/2}

p=-1.5
q=-1.5

plot [][0:1] \
  '../runsL/POST1/l1_typ305_spc.txt'   u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 1 t  "i_p =   5", \
  '../runsL/POST1/l1_typ310_spc.txt'   u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 2 t  "i_p = 10", \
  '../runsK/POST3/k_p2gAdt2_spc.txt'   u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 3 t  "i_p = 20", \
  '../runsL/POST1/l1_typ330_spc.txt'   u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 4 t  "i_p = 30", \
  '../runsL/POST1/l1_typ336_spc.txt'   u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 6 t  "i_p = 36"
  
#==========================================================================

# 1.i  2.Fi   3.nk  4.m2k  5.m3k   5.gk  6.fk  8.nnk  9.2re(ss1k)  10.2re(ss2k) 


set origin 0.5, 2/3.

set ylabel '|a_i|^4 / 2 n_i^2'
set ytics 0.05

plot [][1:1.2] \
  '../runsL/POST1/l1_typ305_mom.txt'   u ($1):($4)   w lp ls 1 t  "i_p =   5", \
  '../runsL/POST1/l1_typ310_mom.txt'   u ($1):($4)   w lp ls 2 t  "i_p = 10", \
  '../runsK/POST3/k_p2gAdt2_mom.txt'   u ($1):($4)   w lp ls 3 t  "i_p = 20", \
  '../runsL/POST1/l1_typ330_mom.txt'   u ($1):($4)   w lp ls 4 t  "i_p = 30", \
  '../runsL/POST1/l1_typ336_mom.txt'   u ($1):($4)   w lp ls 6 t  "i_p = 36"
 
#==========================================================================



set origin 0.5, 1/3.

set ylabel '|a_i|^6 / 6 n_i^3'

plot [][1:1.2] \
  '../runsL/POST1/l1_typ305_mom.txt'   u ($1):($5)   w lp ls 1 t  "i_p =   5", \
  '../runsL/POST1/l1_typ310_mom.txt'   u ($1):($5)   w lp ls 2 t  "i_p = 10", \
  '../runsK/POST3/k_p2gAdt2_mom.txt'   u ($1):($5)   w lp ls 3 t  "i_p = 20", \
  '../runsL/POST1/l1_typ330_mom.txt'   u ($1):($5)   w lp ls 4 t  "i_p = 30", \
  '../runsL/POST1/l1_typ336_mom.txt'   u ($1):($5)   w lp ls 6 t  "i_p = 36"
 
#==========================================================================


set origin 0.5, 0

set ylabel 'g_i'

plot [][-0.05:0.25] \
  '../runsL/POST1/l1_typ305_mom.txt'   u ($1):($6/$8)   w lp ls 1 t  "i_p =   5", \
  '../runsL/POST1/l1_typ310_mom.txt'   u ($1):($6/$8)   w lp ls 2 t  "i_p = 10", \
  '../runsK/POST3/k_p2gAdt2_mom.txt'   u ($1):($6/$8)   w lp ls 3 t  "i_p = 20", \
  '../runsL/POST1/l1_typ330_mom.txt'   u ($1):($6/$8)   w lp ls 4 t  "i_p = 30", \
  '../runsL/POST1/l1_typ336_mom.txt'   u ($1):($6/$8)   w lp ls 6 t  "i_p = 36"
 
#==========================================================================



unset multiplot

set output



