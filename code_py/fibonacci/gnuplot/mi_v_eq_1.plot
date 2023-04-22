set term pdf fsize 8 enhanced color dashed lw 1 size 8, 3

set datafile commentschars '%#'


set style line 1 lt 1 pt 6 ps 1.5  lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 12 ps 1  lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 2 ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 4 ps 1.3   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4 ps 1    lc rgb "#800080"         lw 2

set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'mi_v_eq_1.pdf'

set multiplot


#set logscale y
set grid
unset key

set xlabel 'i_p - i'
set ylabel 'MI,  {/Symbol D}S' 

#1.i0   2.MI      3.S4D       4.S1     5.S2     6.S3     7.S4

s0 = 0.058648  # entropy shift due to finite bin size (dn = 1)
#s0=0

#==========================================================================



set lmargin at screen 0.08
set rmargin at screen 0.40
set tmargin at screen 0.96
set bmargin at screen 0.14


set key outside r # at -5, 0.15

  plot 'mi_data.txt'  \
       u (30-$1):($2)    i 0  w lp ls 1  t '3-mode MI',\
    '' u (30-$1):($2)    i 1  w lp ls 6  t '2-mode MI',\
    '' u (30-$1):($3-3*s0) i 0  w lp ls 5  t 'S(n_{n-2}, n_{i-1}, n_{i}, {/Symbol q}_3)',\
    '' u (30-$1):($3-2*s0) i 1  w lp ls 4  t 'S(n_{i-1}, n_{i}, {/Symbol q}_2)',\
    '' u (30-$1):($4-s0) i 0  w lp ls 3  t 'S(n_{i-2})',\
    '' u (30-$1):($5-s0) i 0  w lp ls 3  t 'S(n_{i-1})',\
    '' u (30-$1):($6-s0) i 0  w lp ls 3  t 'S(n_i)' , \
    '' u (30-$1):($7)    i 0  w lp ls 2  t 'S({/Symbol q}_3)'

unset key

#==========================================================================


set lmargin at screen 0.48
set rmargin at screen 0.80

set ytics (0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
set logscale y
set label 'dashed lines:'     at 16, 0.009
set label '0.04 exp (0.07 x)'  at 16, 0.009*0.8
set label '0.01 exp (0.08 x)'  at 16, 0.009*0.64
set ylabel 'MI'

  plot[0:25][0.005:0.3] 'mi_data.txt'  \
       u (30-$1):($2) i 0  w lp ls 1 t '3-mode MI',\
    '' u (30-$1):($2) i 1  w lp ls 6 t '2-mode MI', \
    0.04*exp(0.07*x)  w l ls 10, \
    0.01*exp(0.08*x)  w l ls 10
    

#==========================================================================




#==========================================================================
unset multiplot

set output



