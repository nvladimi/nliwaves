set term pdf fsize 7 enhanced color dashed lw 1 size 4, 4

set datafile commentschars '%#'



set style line 1 dt 1 pt 4  ps 1.2  lc rgb "#D00040"         lw 2
set style line 2 dt 1 pt 1  ps 1    lc rgb "#FF4000"         lw 2
set style line 3 dt 1 pt 1  ps 1    lc rgb "orange"          lw 2
set style line 4 dt 1 pt 1  ps 1    lc rgb "#60B000"         lw 2
set style line 5 dt 1 pt 8  ps 1.5  lc rgb "#008080"         lw 2
set style line 6 dt 1 pt 6  ps 1.2  lc rgb  "blue"           lw 2
set style line 7 dt 1 pt 1  ps 1    lc rgb "#800080"         lw 2


set style line 11 dt 2 pt 4  ps 0.8  lc rgb "#D00040"         lw 2
set style line 22 dt 2 pt 1  ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 2 pt 1  ps 1    lc rgb "orange"          lw 2
set style line 44 dt 2 pt 1  ps 1    lc rgb "#60B000"         lw 2
set style line 55 dt 2 pt 8  ps 1.0  lc rgb "#008080"         lw 2
set style line 66 dt 2 pt 6  ps 0.8  lc rgb  "blue"           lw 2
set style line 77 dt 2 pt 1  ps 1    lc rgb "#800080"         lw 2

set style line 111 dt 4 pt 4  ps 0.8  lc rgb "#D00040"         lw 2
set style line 222 dt 4 pt 1  ps 1    lc rgb "#FF4000"         lw 2
set style line 333 dt 4 pt 1  ps 1.0  lc rgb "orange"          lw 2
set style line 444 dt 4 pt 1  ps 1    lc rgb "#60B000"         lw 2
set style line 555 dt 4 pt 8  ps 1    lc rgb "#008080"         lw 2
set style line 666 dt 4 pt 6  ps 0.8  lc rgb  "blue"           lw 2
set style line 777 dt 4 pt 1  ps 1    lc rgb "#800080"         lw 2


set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'fig_sqrt_mi.pdf'

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

set xlabel "i - i_p" offset 0, 0.3

ip1=30


#==========================================================================


set tmargin at screen TM1
set bmargin at screen BM1

set lmargin at screen LM1
set rmargin at screen RM1
 

set ytics 0.002

set ylabel 'S' offset 2,0
set key t c  #spacing 1.5 # at -5, 0.15

set xtics 10

  plot [-30:30][-0.01:0]  \
  '../runsK/MI_32x4/mi_S_Km.txt' i 0 u ($1-ip1):($8-$19) w lp ls 6 t '{/Symbol a = 1/2},  60 modes'
      
    
unset title

#--------------------------------------------------------------------------


set lmargin at screen LM2
set rmargin at screen RM2

set key t c  spacing 1.2

set ylabel 'I_{23}' offset 0.5,0

set ytics 0.005


  plot [-30:30][0:0.02] '../runsK/MI_32x4/mi_mi_Km.txt' \
        i 0  u ($1-ip1):($4) w lp ls 6   t 'I_{23}', \
    ''  i 0  u ($1-ip1):($3) w l ls 66  t 'I_{12}', \
    ''  i 0  u ($1-ip1):($5) w l ls 666 t 'I_{31}'

unset key

#--------------------------------------------------------------------------

set tmargin at screen TM2
set bmargin at screen BM2

set lmargin at screen LM1
set rmargin at screen RM1

set key b c  spacing 1.2

set ylabel 'I_{123}' offset 0,0

set ytics 0.05

unset key

plot[-30:30][]  \
    '../runsK/MI_32x4/mi_mi_Km.txt'   u ($1-ip1):2 i 0 w lp ls 6 t '7M, bin =   1/4'


unset logscale
unset label
set ytics auto



#--------------------------------------------------------------------------


set lmargin at screen LM2
set rmargin at screen RM2

set key b c  spacing 1.2

set ylabel 'II' offset 1,0

set ytics 0.03

unset key

plot[-30:30][-0.09:0]  \
    '../runsK/MI_32x4/mi_mi_Km.txt'  u ($1-ip1):6 i 0 w lp ls 6 t '7M, bin =   1/4'
unset logscale
unset label
set ytics auto




#==========================================================================
unset multiplot

set output



