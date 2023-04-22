set term pdf fsize 8 enhanced color dashed lw 2 size 12, 3

set datafile commentschars '%#'



set style line 1 dt 4 pt 6  ps 1    lc rgb "#D00040"         lw 2
set style line 2 dt 4 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 dt 4 pt 4  ps 1.0    lc rgb "orange"          lw 2
set style line 4 dt 4 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 5 dt 4 pt 6  ps 1.0   lc rgb "#008080"         lw 2
set style line 6 dt 4 pt 10  ps 1  lc rgb  "blue"           lw 2
set style line 7 dt 4 pt 4  ps 1    lc rgb "#800080"         lw 2

set style line 11 dt 2 pt 6  ps 1  lc rgb "#D00040"         lw 2
set style line 22 dt 2 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 2 pt 4 ps 1.5    lc rgb "orange"          lw 2
set style line 44 dt 2 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 55 dt 2 pt 6  ps 1.5    lc rgb "#008080"         lw 2
set style line 66 dt 2 pt 10  ps 1.3  lc rgb  "blue"           lw 2
set style line 77 dt 2 pt 2  ps 1    lc rgb "#800080"         lw 2

set style line 111 dt 1 pt 6  ps 1  lc rgb "#D00040"         lw 2
set style line 222 dt 1 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 333 dt 1 pt 4 ps 2.0    lc rgb "orange"          lw 2
set style line 444 dt 1 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 555 dt 1 pt 6  ps 2.0    lc rgb "#008080"         lw 2
set style line 666 dt 1 pt 10  ps 1.3  lc rgb  "blue"           lw 2
set style line 777 dt 1 pt 2  ps 1    lc rgb "#800080"         lw 2


set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'mi_sqrt.pdf'

set multiplot


#==========================================================================


#LM1=0.08;  RM1=0.48
#LM2=0.58;  RM2=0.98
#TM1=0.98;  BM1=0.58
#TM2=0.48;  BM2=0.08

#LM1=0.07;  RM1=0.33
#LM2=0.40;  RM2=0.67
#LM3=0.74;  RM3=0.99
#BM1=0.16;  TM1=0.96


LM1=0.04;  RM1=0.24
LM2=0.29;  RM2=0.49
LM3=0.54;  RM3=0.74
LM4=0.79;  RM4=0.99
TM1=0.96;  BM1=0.16


set key samplen 3

set xlabel font "Times-Italic"
set ylabel font "Times-Italic"
set key font "Times-Italic"

set xlabel "i - i_p" offset 0, 0.3

ip1=30

set grid

#==========================================================================


set tmargin at screen TM1
set bmargin at screen BM1

set lmargin at screen LM1
set rmargin at screen RM1
 

set ytics 0.002

set ylabel 'S' offset 2,0
set key b c  #spacing 1.5 # at -5, 0.15

set xtics 10

  plot [-30:30][]  \
  '../runsK/MI_32x4/mi_S_Km.txt' i 0     u ($1-ip1):($8-$19) w l ls 555 t '7 M, bin = 1/4  ', \
  '../runsK/MI_32x8/mi_S_Km.txt' i 0     u ($1-ip1):($8-$19) w l ls 55  t '7 M, bin = 1/8  ', \
  '../runsK/MI_16x16/mi_S_Km.txt' i 0    u ($1-ip1):($8-$19) w l ls 5   t '7 M, bin = 1/16', \
  '../runsK/MI_16x04a/mi_S_Km.txt' i 0   u ($1-ip1):($8-$19) w p ls 333 t '4 M, bin = 1/4  ', \
  '../runsK/MI_16x08a/mi_S_Km.txt' i 0   u ($1-ip1):($8-$19) w p ls 33  t '4 M, bin = 1/8  ', \
  '../runsK/MI_16x16a/mi_S_Km.txt' i 0   u ($1-ip1):($8-$19) w p ls 3   t '4 M, bin = 1/16'#, \
  '../runsK/MI_16x04b/mi_S_Km.txt' i 0   u ($1-ip1):($8-$19) w p ls 111 t 'bin =   1/4', \
  '../runsK/MI_16x08b/mi_S_Km.txt' i 0   u ($1-ip1):($8-$19) w p ls 11  t 'bin =   1/8', \
  '../runsK/MI_16x16b/mi_S_Km.txt' i 0   u ($1-ip1):($8-$19) w p ls 1   t 'bin = 1/16'


      
    
unset title

#--------------------------------------------------------------------------


set lmargin at screen LM2
set rmargin at screen RM2

set key t c  spacing 1.2

set ylabel 'I_{23}' offset 0.5,0

set ytics 0.02

unset key

  plot [-30:30][]  \
  '../runsK/MI_32x4/mi_mi_Km.txt' i 0    u ($1-ip1):($4) w lp ls 555 t '7M, bin =   1/4', \
  '../runsK/MI_32x8/mi_mi_Km.txt' i 0    u ($1-ip1):($4) w lp ls 55  t '7M, bin =   1/8', \
  '../runsK/MI_16x16/mi_mi_Km.txt' i 0   u ($1-ip1):($4) w lp ls 5   t '7M, bin = 1/16', \
  '../runsK/MI_16x04a/mi_mi_Km.txt' i 0   u ($1-ip1):($4) w lp ls 333 t '4M, bin =   1/4', \
  '../runsK/MI_16x08a/mi_mi_Km.txt' i 0   u ($1-ip1):($4) w lp ls 33  t '4M, bin =   1/8', \
  '../runsK/MI_16x16a/mi_mi_Km.txt' i 0   u ($1-ip1):($4) w lp ls 3   t '4M, bin = 1/16'#, \
  '../runsK/MI_16x04b/mi_mi_Km.txt' i 0   u ($1-ip1):($4) w lp ls 111 t 'bin =   1/4', \
  '../runsK/MI_16x08b/mi_mi_Km.txt' i 0   u ($1-ip1):($4) w lp ls 11  t 'bin =   1/8', \
  '../runsK/MI_16x16b/mi_mi_Km.txt' i 0   u ($1-ip1):($4) w lp ls 1   t 'bin = 1/16'


#--------------------------------------------------------------------------


set lmargin at screen LM3
set rmargin at screen RM3

set key b c  spacing 1.2

set ylabel 'I_{123}' offset 0,0

set ytics 0.2

unset key

plot[-30:30][]  \
    '../runsK/MI_32x4/mi_mi_Km.txt'   u ($1-ip1):2 i 0 w lp ls 555 t '7M, bin =   1/4', \
    '../runsK/MI_32x8/mi_mi_Km.txt'   u ($1-ip1):2 i 0 w lp ls 55  t '7M, bin =   1/8', \
    '../runsK/MI_16x16/mi_mi_Km.txt'  u ($1-ip1):2 i 0 w lp ls 5   t '7M, bin = 1/16',  \
    '../runsK/MI_16x04a/mi_mi_Km.txt'  u ($1-ip1):2 i 0 w lp ls 333 t '7M, bin =   1/4', \
    '../runsK/MI_16x08a/mi_mi_Km.txt'  u ($1-ip1):2 i 0 w lp ls 33  t '7M, bin =   1/8', \
    '../runsK/MI_16x16a/mi_mi_Km.txt'  u ($1-ip1):2 i 0 w lp ls 3   t '7M, bin = 1/16' 


unset logscale
unset label
set ytics auto



#--------------------------------------------------------------------------


set lmargin at screen LM4
set rmargin at screen RM4

set key b c  spacing 1.2

set ylabel 'II' offset 1,0

set ytics 0.2

unset key

plot[-30:30][]  \
    '../runsK/MI_32x4/mi_mi_Km.txt'  u ($1-ip1):6 i 0 w lp ls 555 t '7M, bin =   1/4', \
    '../runsK/MI_32x8/mi_mi_Km.txt'  u ($1-ip1):6 i 0 w lp ls 55  t '7M, bin =   1/8', \
    '../runsK/MI_16x16/mi_mi_Km.txt' u ($1-ip1):6 i 0 w lp ls 5   t '7M, bin = 1/16',  \
    '../runsK/MI_16x04a/mi_mi_Km.txt'  u ($1-ip1):6 i 0 w lp ls 333 t '7M, bin =   1/4', \
    '../runsK/MI_16x08a/mi_mi_Km.txt'  u ($1-ip1):6 i 0 w lp ls 33  t '7M, bin =   1/8', \
    '../runsK/MI_16x16a/mi_mi_Km.txt'  u ($1-ip1):6 i 0 w lp ls 3   t '7M, bin = 1/16' 


unset logscale
unset label
set ytics auto




#==========================================================================
unset multiplot

set output



