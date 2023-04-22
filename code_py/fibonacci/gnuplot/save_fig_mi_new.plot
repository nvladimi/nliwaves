set term pdf fsize 8 enhanced color dashed lw 1 size 7, 7

set datafile commentschars '%#'


set style line 1 lt 1 pt 6  ps 1.5  lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 4  ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 10  ps 1   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4  ps 1    lc rgb "#800080"         lw 2


set style line 11 dt 2 pt 6  ps 1.0  lc rgb "#D00040"         lw 2
set style line 22 dt 2 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 2 pt 4 ps 1     lc rgb "orange"          lw 2
set style line 44 dt 2 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 55 dt 2 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 66 dt 2 pt 10  ps 1.3 lc rgb  "blue"           lw 2
set style line 77 dt 2 pt 2  ps 1    lc rgb "#800080"         lw 2


set style line 111 dt 4 pt 6  ps 0.5  lc rgb "#D00040"         lw 2
set style line 222 dt 4 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 333 dt 4 pt 4 ps 1     lc rgb "orange"          lw 2
set style line 444 dt 4 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 555 dt 4 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 666 dt 4 pt 10  ps 1.3 lc rgb  "blue"           lw 2
set style line 777 dt 4 pt 2  ps 1    lc rgb "#800080"         lw 2




set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'fig_mi_new.pdf'

set multiplot

set grid
unset key

TM1=0.95;  BM1=0.55
TM2=0.45;  BM2=0.05

LM1=0.08;  RM1=0.48
LM2=0.56;  RM2=0.98


#==========================================================================


set tmargin at screen TM1
set bmargin at screen BM1

set lmargin at screen LM1
set rmargin at screen RM1

set key l spacing 1.2

set title "V = 1, effect of ensemble size" 

set ylabel 'MI'
set xlabel 'i_p - i'

set key l spacing 1.2 at 4, 0.15

set xtics 2,8
set ytics 0.04

plot[0:64][-0.04:0.16]  \
 '../runsQ/MIx32x1/mi_mi_Q8.txt' \
      u (50-$1):3 i 0 w lp ls 4 t 'I_{12}', \
   '' u (50-$1):4 i 0 w lp ls 6 t 'I_{23}', \
   '' u (50-$1):5 i 0 w lp ls 3 t 'I_{31}', \
   '' u (50-$1):6 i 0 w lp ls 1 t 'II', \
   '' u (70-$1):3 i 2 w lp ls 4 notitle, \
   '' u (70-$1):4 i 2 w lp ls 6 notitle, \
   '' u (70-$1):5 i 2 w lp ls 3 notitle, \
   '' u (70-$1):6 i 2 w lp ls 1 notitle, \
 '../runsQ/MIx48/mi_mi_Q8.txt' \
      u (50-$1):3 i 0 w lp ls 44 notitle, \
   '' u (50-$1):4 i 0 w lp ls 66 notitle, \
   '' u (50-$1):5 i 0 w lp ls 33 notitle, \
   '' u (50-$1):6 i 0 w lp ls 11 notitle, \
   '' u (70-$1):3 i 2 w lp ls 44 notitle, \
   '' u (70-$1):4 i 2 w lp ls 66 notitle, \
   '' u (70-$1):5 i 2 w lp ls 33 notitle, \
   '' u (70-$1):6 i 2 w lp ls 11 notitle, \
 '../runsQ/MIx48reduced/mi_mi_Q8.txt' \
      u (50-$1):3 i 0 w lp ls 444 notitle, \
   '' u (50-$1):4 i 0 w lp ls 666 notitle, \
   '' u (50-$1):5 i 0 w lp ls 333 notitle, \
   '' u (50-$1):6 i 0 w lp ls 111 notitle, \
   '' u (70-$1):3 i 2 w lp ls 444 notitle, \
   '' u (70-$1):4 i 2 w lp ls 666 notitle, \
   '' u (70-$1):5 i 2 w lp ls 333 notitle, \
   '' u (70-$1):6 i 2 w lp ls 111 notitle
  


#==========================================================================

set title "V = Fi,  effect of ensemble size" 

set lmargin at screen LM2
set rmargin at screen RM2

set ylabel 'MI'
set xlabel 'i - i_p'

set xtics 2,4
set ytics 0.02
set key spacing 1.5 at 8, 0.114


plot[4:48][-0.02:0.12]  \
'../runsQ/MIx32x1/mi_mi_Q8.txt' \
      u ($1-10):3 i 1 w lp ls 6 t 'I_{12}', \
   '' u ($1-10):4 i 1 w lp ls 4 t 'I_{23}', \
   '' u ($1-10):5 i 1 w lp ls 3 t 'I_{31}', \
   '' u ($1-10):6 i 1 w lp ls 1 t 'II', \
'../runsQ/MIx48/mi_mi_Q8.txt' \
      u ($1-10):3 i 1 w lp ls 66 notitle, \
   '' u ($1-10):4 i 1 w lp ls 44 notitle, \
   '' u ($1-10):5 i 1 w lp ls 33 notitle, \
   '' u ($1-10):6 i 1 w lp ls 11 notitle, \
  'mi_data.txt'  \
      u ($1-10):($2) i 7  w lp ls 4 notitle, \
'../runsQ/MIx48reduced/mi_mi_Q8.txt' \
      u ($1-10):3 i 1 w lp ls 665 notitle, \
   '' u ($1-10):4 i 1 w lp ls 444 notitle, \
   '' u ($1-10):5 i 1 w lp ls 333 notitle, \
   '' u ($1-10):6 i 1 w lp ls 111 notitle
   
#==========================================================================
#==========================================================================


set tmargin at screen TM2
set bmargin at screen BM2

set lmargin at screen LM1
set rmargin at screen RM1

set key l spacing 1.2

set title "V = 1, effect of bin size" 

set ylabel 'MI'
set xlabel 'i_p - i'

set key l spacing 1.2 at 4, 0.23

set xtics 2,8
set ytics 0.08

plot[0:64][-0.20:0.24]  \
 '../runsQ/MIx32x1/mi_mi_Q8.txt' \
      u (50-$1):3 i 0 w lp ls 4 t 'I_{12}', \
   '' u (50-$1):4 i 0 w lp ls 6 t 'I_{23}', \
   '' u (50-$1):5 i 0 w lp ls 3 t 'I_{31}', \
   '' u (50-$1):6 i 0 w lp ls 1 t 'II', \
   '' u (70-$1):3 i 2 w lp ls 4 notitle, \
   '' u (70-$1):4 i 2 w lp ls 6 notitle, \
   '' u (70-$1):5 i 2 w lp ls 3 notitle, \
   '' u (70-$1):6 i 2 w lp ls 1 notitle, \
 '../runsQ/MIx32x4/mi_mi_Q8.txt' \
      u (50-$1):3 i 0 w lp ls 44 notitle, \
   '' u (50-$1):4 i 0 w lp ls 66 notitle, \
   '' u (50-$1):5 i 0 w lp ls 33 notitle, \
   '' u (50-$1):6 i 0 w lp ls 11 notitle, \
   '' u (70-$1):3 i 2 w lp ls 44 notitle, \
   '' u (70-$1):4 i 2 w lp ls 66 notitle, \
   '' u (70-$1):5 i 2 w lp ls 33 notitle, \
   '' u (70-$1):6 i 2 w lp ls 11 notitle, \
 '../runsQ/MIx32x8/mi_mi_Q8.txt' \
      u (50-$1):3 i 0 w lp ls 444 notitle, \
   '' u (50-$1):4 i 0 w lp ls 666 notitle, \
   '' u (50-$1):5 i 0 w lp ls 333 notitle, \
   '' u (50-$1):6 i 0 w lp ls 111 notitle #, \
   '' u (70-$1):3 i 2 w lp ls 444 notitle, \
   '' u (70-$1):4 i 2 w lp ls 666 notitle, \
   '' u (70-$1):5 i 2 w lp ls 333 notitle, \
   '' u (70-$1):6 i 2 w lp ls 111 notitle
  


#==========================================================================

set title "V = Fi,  effect of bin size" 

set lmargin at screen LM2
set rmargin at screen RM2

set ylabel 'MI'
set xlabel 'i - i_p'

set xtics 2,4
set ytics 0.04
set key spacing 1.5 at 8, 0.19


plot[4:48][-0.12:0.20]  \
'../runsQ/MIx32x1/mi_mi_Q8.txt' \
      u ($1-10):3 i 1 w lp ls 6 t 'I_{12}', \
   '' u ($1-10):4 i 1 w lp ls 4 t 'I_{23}', \
   '' u ($1-10):5 i 1 w lp ls 3 t 'I_{31}', \
   '' u ($1-10):6 i 1 w lp ls 1 t 'II', \
'../runsQ/MIx32x4/mi_mi_Q8.txt' \
      u ($1-10):3 i 1 w lp ls 66 notitle, \
   '' u ($1-10):4 i 1 w lp ls 44 notitle, \
   '' u ($1-10):5 i 1 w lp ls 33 notitle, \
   '' u ($1-10):6 i 1 w lp ls 11 notitle, \
'../runsQ/MIx32x8/mi_mi_Q8.txt' \
      u ($1-10):3 i 1 w lp ls 665 notitle, \
   '' u ($1-10):4 i 1 w lp ls 444 notitle, \
   '' u ($1-10):5 i 1 w lp ls 333 notitle, \
   '' u ($1-10):6 i 1 w lp ls 111 notitle
   


#==========================================================================
unset multiplot

set output



