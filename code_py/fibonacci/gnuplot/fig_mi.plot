set term pdf fsize 7 enhanced color dashed lw 1 size 4, 4

set datafile commentschars '%#'



set style line 1 dt 1 pt 6  ps 1.5  lc rgb "#D00040"         lw 2
set style line 2 dt 1 pt 1  ps 1    lc rgb "#FF4000"         lw 2
set style line 3 dt 1 pt 1  ps 1    lc rgb "orange"          lw 2
set style line 4 dt 1 pt 1  ps 1    lc rgb "#60B000"         lw 2
set style line 5 dt 1 pt 1  ps 1    lc rgb "#008080"         lw 2
set style line 6 dt 1 pt 4  ps 1.5  lc rgb  "blue"           lw 2
set style line 7 dt 1 pt 1  ps 1    lc rgb "#800080"         lw 2


set style line 11 dt 2 pt 6  ps 0.8  lc rgb "#D00040"         lw 2
set style line 22 dt 2 pt 1  ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 2 pt 1  ps 1    lc rgb "orange"          lw 2
set style line 44 dt 2 pt 1  ps 1    lc rgb "#60B000"         lw 2
set style line 55 dt 2 pt 1  ps 1    lc rgb "#008080"         lw 2
set style line 66 dt 2 pt 4  ps 0.8  lc rgb  "blue"           lw 2
set style line 77 dt 2 pt 1  ps 1    lc rgb "#800080"         lw 2

set style line 111 dt 4 pt 1  ps 1.5  lc rgb "#D00040"         lw 2
set style line 222 dt 4 pt 1  ps 1    lc rgb "#FF4000"         lw 2
set style line 333 dt 4 pt 1  ps 1    lc rgb "orange"          lw 2
set style line 444 dt 4 pt 1  ps 1    lc rgb "#60B000"         lw 2
set style line 555 dt 4 pt 1  ps 1    lc rgb "#008080"         lw 2
set style line 666 dt 4 pt 2  ps 1.5  lc rgb  "blue"           lw 2
set style line 777 dt 4 pt 1  ps 1    lc rgb "#800080"         lw 2


set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1


set output 'fig_mi.pdf'

set multiplot


#==========================================================================


LM1=0.10;  RM1=0.48
LM2=0.60;  RM2=0.98
TM1=0.98;  BM1=0.58
TM2=0.48;  BM2=0.08

set key samplen 3

set xlabel font "Times-Italic"
set ylabel font "Times-Italic"
set key font "Times-Italic"


set xlabel "{/Symbol |} i - p {/Symbol |}" offset 0, 0.3

ip1=50
ip2=70
ip3=10
set xtics 20
set xrange [0:64]

#==========================================================================


set tmargin at screen TM1
set bmargin at screen BM1

set lmargin at screen LM1
set rmargin at screen RM1

set ytics 0.1

set ylabel 'S' offset 0,0
set key b l spacing 1.2

plot [][]  \
  '../runsQ/MIx32x1/mi_S_Q8.txt' i 2  u (ip2-$1):($8-$19) w lp ls 1 t '{/Symbol a} = 0,  80 modes', \
  '../runsQ/MIx32x1/mi_S_Q8.txt' i 1  u ($1-ip3):($8-$19) w lp ls 6 t '{/Symbol a} = 1,  60 modes'
  
unset label      
    
unset title

#--------------------------------------------------------------------------


set lmargin at screen LM2
set rmargin at screen RM2

set key t c  spacing 1.2

set ylabel 'I' offset 1,0

set ytics 0.05

set key l

plot [][]  '../runsQ/MIx32x1/mi_mi_Q8.txt' \
       i 2  u (ip2-$1):($4) w lp ls 1     t   '{/Symbol a} = 0,  I_{23}', \
    '' i 2  u (ip2-$1):($5) w lp ls 111   t   '{/Symbol a} = 0,  I_{31}', \
    '' i 1  u ($1-ip3):($4) w lp ls 6     t   '{/Symbol a} = 1,  I_{23}', \
    '' i 1  u ($1-ip3):($5) w lp ls 666   t   '{/Symbol a} = 1,  I_{31}'


#--------------------------------------------------------------------------


set tmargin at screen TM2
set bmargin at screen BM2

set lmargin at screen LM1
set rmargin at screen RM1

set ylabel 'I_{123}' offset 0,0
set ytics 0.1

set key b r

plot [][]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 2  u (ip2-$1):($2) w lp ls 1 t '{/Symbol a} = 0', \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 1  u ($1-ip3):($2) w lp ls 6 t '{/Symbol a} = 1'

#--------------------------------------------------------------------------


set lmargin at screen LM2
set rmargin at screen RM2

set ylabel 'II / I_{123}' offset 0,0

set ytics 0.2
#set ytics auto

plot [][]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 2  u (ip2-$1):($6/$2) w lp ls 1  t '{/Symbol a} = 0', \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 1  u ($1-ip3):($6/$2) w lp ls 6  t '{/Symbol a} = 1', \
0 w l ls 0 notitle

#==========================================================================
unset multiplot

set output



