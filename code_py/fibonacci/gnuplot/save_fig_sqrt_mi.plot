set term pdf fsize 7 enhanced color dashed lw 2 size 6, 2

set datafile commentschars '%#'



set style line 1 lt 1 pt 6  ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 4  ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 10  ps 1  lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4  ps 1    lc rgb "#800080"         lw 2


set style line 11 dt 2 pt 6  ps 1  lc rgb "#D00040"         lw 2
set style line 22 dt 2 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 2 pt 4 ps 1    lc rgb "orange"          lw 2
set style line 44 dt 2 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 55 dt 2 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 66 dt 2 pt 10  ps 1.3  lc rgb  "blue"           lw 2
set style line 77 dt 2 pt 2  ps 1    lc rgb "#800080"         lw 2

set style line 111 dt 4 pt 6  ps 1  lc rgb "#D00040"         lw 2
set style line 222 dt 4 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 333 dt 4 pt 4 ps 1    lc rgb "orange"          lw 2
set style line 444 dt 4 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 555 dt 4 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 666 dt 4 pt 10  ps 1.3  lc rgb  "blue"           lw 2
set style line 777 dt 4 pt 2  ps 1    lc rgb "#800080"         lw 2





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

LM1=0.07;  RM1=0.33
LM2=0.40;  RM2=0.67
LM3=0.74;  RM3=0.99



BM1=0.16;  TM1=0.96

set tmargin at screen TM1
set bmargin at screen BM1

set key samplen 3

set xlabel font "Times-Italic"
set ylabel font "Times-Italic"
set key font "Times-Italic"


#==========================================================================

set lmargin at screen LM1
set rmargin at screen RM1

ip1=30 
ip2=30

I1=12
I2=13
I3=14
I4=15

set ytics 0.002

set xlabel "i - i_p" offset 0, 0.3
set ylabel 'S' offset 2,0
set key b c  #spacing 1.5 # at -5, 0.15

set xtics 10



  plot [-30:30][]  \
  '../runsK/MI_32x4/mi_S_Km.txt' i 0    u ($1-30):($8-$19) w l ls 555 t 'bin =   1/4', \
  '../runsK/MI_32x8/mi_S_Km.txt' i 0    u ($1-30):($8-$19) w l ls 55  t 'bin =   1/8', \
  '../runsK/MI_16x16/mi_S_Km.txt' i 0   u ($1-30):($8-$19) w l ls 5   t 'bin = 1/16'


      
    
unset title

#--------------------------------------------------------------------------

set lmargin at screen LM2
set rmargin at screen RM2

set key t c  spacing 1.2

set xlabel 'i - i_p' offset 0,0.3

set ylabel 'MI' offset 1,0

set ytics 0.02



plot[-30:30][]  \
'../runsK/MI_32x4/mi_mi_Km.txt' \
      u ($1-30):3 i 0 w l ls 666 notitle, \
   '' u ($1-30):4 i 0 w l ls 445 notitle, \
   '' u ($1-30):5 i 0 w l ls 333 notitle, \
'../runsK/MI_32x8/mi_mi_Km.txt' \
      u ($1-30):3 i 0 w l ls 66 notitle, \
   '' u ($1-30):4 i 0 w l ls 44 notitle, \
   '' u ($1-30):5 i 0 w l ls 33 notitle, \
'../runsK/MI_16x16/mi_mi_Km.txt' \
      u ($1-30):3 i 0 w l ls 6 t 'I_{12}', \
   '' u ($1-30):4 i 0 w l ls 4 t 'I_{23}', \
   '' u ($1-30):5 i 0 w l ls 3 t 'I_{31}'


#--------------------------------------------------------------------------

set lmargin at screen LM3
set rmargin at screen RM3

set key b c  spacing 1.2

set xlabel 'i - i_p' offset 0,0.3

set ylabel 'II' offset 1,0

set ytics 0.2

unset key

plot[-30:30][]  \
    '../runsK/MI_32x4/mi_mi_Km.txt'  u ($1-30):6 i 0 w l ls 111 t 'bin =   1/4', \
    '../runsK/MI_32x8/mi_mi_Km.txt'  u ($1-30):6 i 0 w l ls 11  t 'bin =   1/8', \
    '../runsK/MI_16x16/mi_mi_Km.txt' u ($1-30):6 i 0 w l ls 1   t 'bin = 1/16' 


unset logscale
unset label
set ytics auto




#==========================================================================
unset multiplot

set output



