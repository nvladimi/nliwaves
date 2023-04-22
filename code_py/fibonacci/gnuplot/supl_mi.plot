set term pdf fsize 7 enhanced color dashed lw 1 size 8, 6

set datafile commentschars '%#'



set style line 1 dt 4 pt 4  ps 1.5    lc rgb "#D00040"         lw 2
set style line 2 dt 1 pt 1  ps 1    lc rgb "#FF4000"         lw 2
set style line 3 dt 1 pt 1  ps 1    lc rgb "orange"          lw 2
set style line 4 dt 1 pt 1  ps 1    lc rgb "#60B000"         lw 2
set style line 5 dt 1 pt 1  ps 1   lc rgb "#008080"         lw 2
set style line 6 dt 4 pt 12  ps 1.5  lc rgb  "blue"           lw 2
set style line 7 dt 1 pt 1  ps 1    lc rgb "#800080"         lw 2

set style line 11 dt 2 pt 6  ps 1  lc rgb "#D00040"         lw 2
set style line 22 dt 2 pt 1  ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 2 pt 4  ps 1    lc rgb "orange"          lw 2
set style line 44 dt 2 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 55 dt 2 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 66 dt 2 pt 1  ps 1  lc rgb  "blue"           lw 2
set style line 77 dt 2 pt 2  ps 1    lc rgb "#800080"         lw 2

set style line 111 dt 4 pt 7  ps 0.5  lc rgb "#D00040"         lw 2
set style line 222 dt 4 pt 1  ps 1    lc rgb "#FF4000"         lw 2
set style line 333 dt 4 pt 7  ps 1    lc rgb "orange"          lw 2
set style line 444 dt 4 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 555 dt 4 pt 7  ps 1    lc rgb "#008080"         lw 2
set style line 666 dt 4 pt 7 ps 0.5  lc rgb  "blue"           lw 2
set style line 777 dt 4 pt 2  ps 1    lc rgb "#800080"         lw 2


set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'supl_mi.pdf'

set multiplot

LM1=0.04;  RM1=0.24
LM2=0.29;  RM2=0.49
LM3=0.54;  RM3=0.74
LM4=0.79;  RM4=0.99

TM1=0.99;  BM1=0.72
TM2=0.66;  BM2=0.39
TM3=0.33;  BM3=0.06

set key samplen 3

set xlabel font "Times-Italic"
set ylabel font "Times-Italic"
set key font "Times-Italic"

LS1 = 1
LS11 = 111
LS2 = 6
LS22 = 666

#==========================================================================

set xlabel "i - i_p" offset 0, 0.3

ip1=30

set xtics 10

set tmargin at screen TM1
set bmargin at screen BM1

#==========================================================================

set lmargin at screen LM1
set rmargin at screen RM1
 
set ytics 0.005

set ylabel 'S' offset 3,0
set key c c
set label "{/Symbol a = 1/2}" at 17,-0.0008 


plot [-30:30][:0]  \
  '../runsK/MI_16x02b/mi_S_Km.txt' i 0  u ($1-ip1):($8-$19) w lp ls LS1   t '52 M, bin = 1/2', \
  '../runsK/MI_16x04b/mi_S_Km.txt' i 0  u ($1-ip1):($8-$19) w lp ls LS11  t '52 M, bin = 1/4', \
  '../runsK/MI_16x02a/mi_S_Km.txt' i 0  u ($1-ip1):($8-$19) w lp ls LS2   t '28 M, bin = 1/2', \
  '../runsK/MI_16x04a/mi_S_Km.txt' i 0  u ($1-ip1):($8-$19) w lp ls LS22  t '28 M, bin = 1/4'

unset label
unset key

#--------------------------------------------------------------------------


set lmargin at screen LM2
set rmargin at screen RM2

set ylabel 'I_{23}' offset 4,0

set ytics 0.001

  plot [-30:30][0:0.003]  \
  '../runsK/MI_16x02b/mi_mi_Km.txt' i 0  u ($1-ip1):($4) w lp ls LS1   t '52 M, bin = 1/2', \
  '../runsK/MI_16x04b/mi_mi_Km.txt' i 0  u ($1-ip1):($4) w lp ls LS11  t '52 M, bin = 1/4', \
  '../runsK/MI_16x02a/mi_mi_Km.txt' i 0  u ($1-ip1):($4) w lp ls LS2   t '28 M, bin = 1/2', \
  '../runsK/MI_16x04a/mi_mi_Km.txt' i 0  u ($1-ip1):($4) w lp ls LS22  t '28 M, bin = 1/4'
  
#--------------------------------------------------------------------------


set lmargin at screen LM3
set rmargin at screen RM3

set ylabel 'I_{123}' offset 2,0

set ytics 0.01

plot[-30:30][0:0.03]  \
   '../runsK/MI_16x02b/mi_mi_Km.txt' u ($1-ip1):2 i 0 w lp ls LS1   t '52 M, bin = 1/2', \
   '../runsK/MI_16x04b/mi_mi_Km.txt' u ($1-ip1):2 i 0 w lp ls LS11  t '52 M, bin = 1/4', \
   '../runsK/MI_16x02a/mi_mi_Km.txt' u ($1-ip1):2 i 0 w lp ls LS2   t '28 M, bin = 1/2', \
   '../runsK/MI_16x04a/mi_mi_Km.txt' u ($1-ip1):2 i 0 w lp ls LS22  t '28 M, bin = 1/4'

#--------------------------------------------------------------------------


set lmargin at screen LM4
set rmargin at screen RM4

set ylabel 'II' offset 2,0

set ytics 0.01

plot[-30:30][-0.03:0]  \
   '../runsK/MI_16x02b/mi_mi_Km.txt' u ($1-ip1):6 i 0 w lp ls LS1   t '52M, bin = 1/2', \
   '../runsK/MI_16x04b/mi_mi_Km.txt' u ($1-ip1):6 i 0 w lp ls LS11  t '52M, bin = 1/4', \
   '../runsK/MI_16x02a/mi_mi_Km.txt' u ($1-ip1):6 i 0 w lp ls LS2   t '28M, bin = 1/2', \
   '../runsK/MI_16x04a/mi_mi_Km.txt' u ($1-ip1):6 i 0 w lp ls LS22  t '28M, bin = 1/4'


#==========================================================================



set xlabel "i_p - i" offset 0, 0.3

ip1=50
set xtics 2,8
set xrange [0:44]

set tmargin at screen TM2
set bmargin at screen BM2


#==========================================================================

set lmargin at screen LM1
set rmargin at screen RM1

set ytics 0.1

set ylabel 'S' offset 2,0
set key b l  #spacing 1.5 # at -5, 0.15

set label "{/Symbol a = 0}" at 36, -0.03

plot [][]  \
  '../runsQ/MIx32x1/mi_S_Q8.txt' i 0    u (ip1-$1):($8-$19) w lp ls LS1   t '60 M, bin = 1   ', \
  '../runsQ/MIx32x4/mi_S_Q8.txt' i 0    u (ip1-$1):($8-$19) w lp ls LS11  t '60 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_S_Q8.txt' i 0   u (ip1-$1):($8-$19) w lp ls LS2   t '30 M, bin = 1   ', \
  '../runsQ/MIx32x4a/mi_S_Q8.txt' i 0   u (ip1-$1):($8-$19) w lp ls LS22  t '30 M, bin = 1/4'

unset label
unset key    


#--------------------------------------------------------------------------


set lmargin at screen LM2
set rmargin at screen RM2

set ylabel 'I_{23}' offset 1,0

set ytics 0.05

plot [][]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 0    u (ip1-$1):($4) w lp ls LS1 t '60 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 0    u (ip1-$1):($4) w lp ls LS11  t '60 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 0   u (ip1-$1):($4) w lp ls LS2 t '30 M, bin =  1 ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 0   u (ip1-$1):($4) w lp ls LS22  t '30 M, bin = 1/4'

#--------------------------------------------------------------------------

set lmargin at screen LM3
set rmargin at screen RM3

set ylabel 'I_{123}' offset 1,0

set ytics 0.2

plot [][]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 0    u (ip1-$1):($2) w lp ls LS1 t '60 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 0    u (ip1-$1):($2) w lp ls LS11  t '60 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 0   u (ip1-$1):($2) w lp ls LS2 t '30 M, bin =  1 ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 0   u (ip1-$1):($2) w lp ls LS22  t '30 M, bin = 1/4'

#--------------------------------------------------------------------------

set lmargin at screen LM4
set rmargin at screen RM4

set ylabel 'II' offset 2,0

set ytics 0.05

plot [][]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 0    u (ip1-$1):($6) w lp ls LS1 t '60 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 0    u (ip1-$1):($6) w lp ls LS11  t '60 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 0   u (ip1-$1):($6) w lp ls LS2 t '30 M, bin =  1 ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 0   u (ip1-$1):($6) w lp ls LS22  t '30 M, bin = 1/4'

#=========================================================================


set xlabel "i - i_p" offset 0, 0.3

ip1=10
set xtics 6,8
set xrange [4:48]

set tmargin at screen TM3
set bmargin at screen BM3

#==========================================================================


set lmargin at screen LM1
set rmargin at screen RM1

set ytics 0.1
set label "{/Symbol a = 1}" at 40, -0.03

set ylabel 'S' offset 2,0
set key b l  #spacing 1.5 # at -5, 0.15

plot [][]  \
  '../runsQ/MIx32x1/mi_S_Q8.txt' i 1    u ($1-ip1):($8-$19) w lp ls LS1   t '100 M, bin = 1   ', \
  '../runsQ/MIx32x4/mi_S_Q8.txt' i 1    u ($1-ip1):($8-$19) w lp ls LS11  t '100 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_S_Q8.txt' i 1   u ($1-ip1):($8-$19) w lp ls LS2   t '  30 M, bin = 1   ', \
  '../runsQ/MIx32x4a/mi_S_Q8.txt' i 1   u ($1-ip1):($8-$19) w lp ls LS22  t '  30 M, bin = 1/4'

unset label      
unset key

unset title

#--------------------------------------------------------------------------


set lmargin at screen LM2
set rmargin at screen RM2

set ylabel 'I_{23}' offset 1,0

set ytics 0.05

plot [][]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 1    u ($1-ip1):($4) w lp ls LS1   t '100 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 1    u ($1-ip1):($4) w lp ls LS11  t '100 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 1   u ($1-ip1):($4) w lp ls LS2   t '30 M, bin =  1    ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 1   u ($1-ip1):($4) w lp ls LS22  t '30 M, bin = 1/4   '


#--------------------------------------------------------------------------


set lmargin at screen LM3
set rmargin at screen RM3

set ylabel 'I_{123}' offset 1,0

set ytics 0.2

plot [][]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 1    u ($1-ip1):($2) w lp ls LS1 t '100 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 1    u ($1-ip1):($2) w lp ls LS11  t '100 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 1   u ($1-ip1):($2) w lp ls LS2 t '30 M, bin =  1 ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 1   u ($1-ip1):($2) w lp ls LS22  t '30 M, bin = 1/4'

#--------------------------------------------------------------------------


set lmargin at screen LM4
set rmargin at screen RM4

set ylabel 'II' offset 2,0

set ytics 0.05

plot [][]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 1    u ($1-ip1):($6) w lp ls LS1 t '100 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 1    u ($1-ip1):($6) w lp ls LS11  t '100 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 1   u ($1-ip1):($6) w lp ls LS2 t '30 M, bin =  1 ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 1   u ($1-ip1):($6) w lp ls LS22  t '30 M, bin = 1/4'


#=========================================================================
unset multiplot

set output



