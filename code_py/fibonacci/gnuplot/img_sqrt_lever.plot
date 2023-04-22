set term pdf fsize 7 enhanced color dashed lw 1.5 size 2.0, 2.0

set datafile commentschars '%#'

PS =0.4

set style line 1 lt 1 pt  7 ps PS    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt  7 ps PS    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt  7 ps PS   lc rgb "orange"          lw 2
set style line 4 lt 1 pt  7 ps PS    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt  7 ps PS    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt  7 ps PS   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt  7 ps PS    lc rgb "#800080"         lw 2


set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'img_sqrt_lever.pdf'

set multiplot
set size 1, 1



unset key

set xtics 10

#==========================================================================


set ylabel 'occupation numbers'
set xlabel 'mode number' offset 0, 0.5
set tmargin 0.7
set rmargin 1.5
set lmargin 5.5

set ytics 1
p2=(0.01*6765)**(-2/3.)

# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 


plot [][0:4] \
  '../runsL/POST1/l1_typ305_spc.txt'   u ($1):($3*p2)   w lp ls 1 t  "i_p =   5", \
  '../runsL/POST1/l1_typ310_spc.txt'   u ($1):($3*p2)   w lp ls 2 t  "i_p = 10", \
  '../runsK/POST3/k_p2gAdt2_spc.txt'   u ($1):($3*p2)   w lp ls 3 t  "i_p = 20", \
  '../runsL/POST1/l1_typ330_spc.txt'   u ($1):($3*p2)   w lp ls 4 t  "i_p = 30", \
  '../runsL/POST1/l1_typ336_spc.txt'   u ($1):($3*p2)   w lp ls 6 t  "i_p = 36"

#==========================================================================



#==========================================================================

unset multiplot

set output



