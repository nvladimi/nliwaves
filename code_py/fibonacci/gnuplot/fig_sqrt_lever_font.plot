set term pdf fsize 8 enhanced color dashed lw 1.5 size 8, 2.6

set datafile commentschars '%#'


set style line 1 lt 1 pt  8 ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt  6 ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt  4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt  6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 10 ps 1   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt  4 ps 1    lc rgb "#800080"         lw 2


set style line 11 lt 3 dt 3 pt  8 ps 1    lc rgb "#D00040"         lw 3
set style line 22 lt 3 dt 3 pt 12 ps 1    lc rgb "#FF4000"         lw 3
set style line 33 lt 3 dt 3 pt  6 ps 1    lc rgb "orange"          lw 3
set style line 44 lt 3 dt 3 pt  4 ps 1    lc rgb "#60B000"         lw 3
set style line 55 lt 3 dt 3 pt  6 ps 1    lc rgb "#008080"         lw 3
set style line 66 lt 3 dt 3 pt 10 ps 1    lc rgb  "blue"           lw 3
set style line 77 lt 3 dt 3 pt  4 ps 1    lc rgb "#800080"         lw 3

set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'fig_sqrt_lever_font.pdf'

set multiplot
set size 1/3., 1


set xlabel font "Times-Italic"
set ylabel font "Times-Italic"
set key font "Times-Italic"


set xlabel 'i'

unset key

# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 



LM1=0.04;  RM1=0.32
LM2=0.38;  RM2=0.66
LM3=0.71;  RM3=0.99

BM1=0.14;  TM1=0.96

set tmargin at screen TM1
set bmargin at screen BM1

set xtics 10

#==========================================================================

set lmargin at screen LM1
set rmargin at screen RM1

set ylabel 'n_i F_i  {/Symbol P}_@p^{-2/3}'


set ytics 1
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

set ylabel '{/Symbol P} / {/Symbol P}_p' offset 1,0

set lmargin at screen LM2
set rmargin at screen RM2

set ytics auto

p2=(0.01*6765)

set key left spacing 1.0

plot [][-1:1] \
  '../runsL/POST1/l1_typ305_spc.txt'   u ($1):($8/p2)   w lp ls 1 t  "p =   5", \
  '../runsL/POST1/l1_typ310_spc.txt'   u ($1):($8/p2)   w lp ls 2 t  "p = 10", \
  '../runsK/POST3/k_p2gAdt2_spc.txt'   u ($1):($8/p2)   w lp ls 3 t  "p = 20", \
  '../runsL/POST1/l1_typ330_spc.txt'   u ($1):($8/p2)   w lp ls 4 t  "p = 30", \
  '../runsL/POST1/l1_typ336_spc.txt'   u ($1):($8/p2)   w lp ls 6 t  "p = 36"

unset key

#==========================================================================


set ylabel '{/Symbol | x |}' offset 1.5,0

set lmargin at screen LM3
set rmargin at screen RM3

unset key

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

set tmargin at screen 0.93  
set bmargin at screen 0.44
set lmargin at screen 0.78
set rmargin at screen 0.94

unset ylabel
unset key

p=-1.5
q=-1.5


set logscale

set xtics (1, 2, 4, 10, 20, 40)
set ytics (0.1, 0.2, 0.5,  1, 2)
set xlabel 'i-1,   m-i' offset 0,0.5

m1=1; m2a=40; m2b=60;

set key t r 

plot [1:40][0.03:2] \
  '../runsK/POST3/k_p2gAdt2_spc.txt'   u ($1-m1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls  3 notitle, \
  '../runsK/POST3/k_p2gAdt2_spc.txt'   u (m2a-$1):(abs($8)*($3/$2)**p*($2)**q) w lp ls 33  notitle, \
  '../runsL/POST1/l1_typ330_spc.txt'   u ($1-m1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 4 notitle, \
  '../runsL/POST1/l1_typ336_spc.txt'   u ($1-m1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 6 notitle, \
  '../runsL/POST1/l1_typ305_spc.txt'   u (m2a-$1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 11 notitle, \
  '../runsL/POST1/l1_typ310_spc.txt'   u (m2a-$1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 22 notitle, \
 2/x w l ls 10 t '\~ 1 / i'


#==========================================================================

unset multiplot

set output



