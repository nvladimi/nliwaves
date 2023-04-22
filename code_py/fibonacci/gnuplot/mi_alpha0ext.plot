set term pdf fsize 8 enhanced color dashed lw 2 size 12, 3

set datafile commentschars '%#'



set style line 1 dt 4 pt 2  ps 1.0    lc rgb "#D00040"         lw 2
set style line 2 dt 4 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 dt 4 pt 4  ps 1.0    lc rgb "orange"          lw 2
set style line 4 dt 4 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 5 dt 4 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 6 dt 4 pt 10  ps 1  lc rgb  "blue"           lw 2
set style line 7 dt 4 pt 4  ps 1    lc rgb "#800080"         lw 2


set style line 11 dt 2 pt 2  ps 1.5  lc rgb "#D00040"         lw 2
set style line 22 dt 2 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 2 pt 4 ps 1.5    lc rgb "orange"          lw 2
set style line 44 dt 2 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 55 dt 2 pt 6  ps 1.5    lc rgb "#008080"         lw 2
set style line 66 dt 2 pt 1  ps 1.3  lc rgb  "blue"           lw 2
set style line 77 dt 2 pt 2  ps 1    lc rgb "#800080"         lw 2

set style line 111 dt 1 pt 2  ps 2.0  lc rgb "#D00040"         lw 2
set style line 222 dt 1 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 333 dt 1 pt 4 ps 2.0    lc rgb "orange"          lw 2
set style line 444 dt 1 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 555 dt 1 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 666 dt 1 pt 10  ps 1.3  lc rgb  "blue"           lw 2
set style line 777 dt 1 pt 2  ps 1    lc rgb "#800080"         lw 2





set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'mi_alpha0ext.pdf'

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


set xlabel "i_p - i" offset 0, 0.3

ip1=50
ip2=70
set xtics 2,8
set xrange [0:64]
set grid

#==========================================================================


set tmargin at screen TM1
set bmargin at screen BM1

set lmargin at screen LM1
set rmargin at screen RM1

set ytics 0.1

set label "{/Symbol a = 0}" at 54, -0.05

set ylabel 'S' offset 1,0
set key b l  #spacing 1.5 # at -5, 0.15

plot [][]  \
  '../runsQ/MIx32x1/mi_S_Q8.txt' i 0    u (ip1-$1):($8-$19) w l ls 555 t '60 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_S_Q8.txt' i 0    u (ip1-$1):($8-$19) w l ls 55  t '60 M, bin = 1/4', \
  '../runsQ/MIx32x8/mi_S_Q8.txt' i 0    u (ip1-$1):($8-$19) w l ls 5   t '60 M, bin = 1/8', \
  '../runsQ/MIx32x1a/mi_S_Q8.txt' i 0   u (ip1-$1):($8-$19) w p ls 333 t '30 M, bin =  1  ', \
  '../runsQ/MIx32x4a/mi_S_Q8.txt' i 0   u (ip1-$1):($8-$19) w p ls 33  t '30 M, bin = 1/4', \
  '../runsQ/MIx32x8a/mi_S_Q8.txt' i 0   u (ip1-$1):($8-$19) w p ls 3   t '30 M, bin = 1/8', \
  '../runsQ/MIx32x1/mi_S_Q8.txt' i 2    u (ip2-$1):($8-$19) w l ls 666 t '240 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_S_Q8.txt' i 2    u (ip2-$1):($8-$19) w l ls 66  t '240 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_S_Q8.txt' i 2   u (ip2-$1):($8-$19) w p ls 111 t '120 M, bin =  1  ', \
  '../runsQ/MIx32x4a/mi_S_Q8.txt' i 2   u (ip2-$1):($8-$19) w p ls 11  t '120 M, bin = 1/4', \
  '../runsQ/MIx32x8a/mi_S_Q8.txt' i 2   u (ip2-$1):($8-$19) w lp ls 1  t '120 M, bin = 1/8'

unset label      
    
unset title

#--------------------------------------------------------------------------


set lmargin at screen LM2
set rmargin at screen RM2

set key t c  spacing 1.2


set ylabel 'I_{23}' offset 0.5,0

set ytics 0.1

unset key 

plot [][]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 0    u (ip1-$1):($4) w l ls 555 t '60 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 0    u (ip1-$1):($4) w lp ls 55  t '60 M, bin = 1/4', \
  '../runsQ/MIx32x8/mi_mi_Q8.txt' i 0    u (ip1-$1):($4) w lp ls 5   t '60 M, bin = 1/8', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 0    u (ip1-$1):($4) w p ls 333 t '30 M, bin =  1 ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 0    u (ip1-$1):($4) w lp ls 33  t '30 M, bin = 1/4', \
  '../runsQ/MIx32x8a/mi_mi_Q8.txt' i 0    u (ip1-$1):($4) w lp ls 3   t '30 M, bin = 1/8', \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 2    u (ip2-$1):($4) w l ls 666 t '240 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 2    u (ip2-$1):($4) w lp ls 66  t '240 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 2   u (ip2-$1):($4) w p ls 111 t '120 M, bin =  1  ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 2   u (ip2-$1):($4) w lp ls 11  t '120 M, bin = 1/4', \
  '../runsQ/MIx32x8a/mi_mi_Q8.txt' i 2   u (ip2-$1):($4) w lp ls 1  t '120 M, bin = 1/8'



#--------------------------------------------------------------------------


set lmargin at screen LM3
set rmargin at screen RM3

set key b c  spacing 1.2


set ylabel 'I_{123}' offset 0,0

set ytics 0.2

unset key

plot [][]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 0    u (ip1-$1):($2) w l ls 555 t '60 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 0    u (ip1-$1):($2) w lp ls 55  t '60 M, bin = 1/4', \
  '../runsQ/MIx32x8/mi_mi_Q8.txt' i 0    u (ip1-$1):($2) w lp ls 5   t '60 M, bin = 1/8', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 0    u (ip1-$1):($2) w p ls 333 t '30 M, bin =  1 ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 0    u (ip1-$1):($2) w lp ls 33  t '30 M, bin = 1/4', \
  '../runsQ/MIx32x8a/mi_mi_Q8.txt' i 0    u (ip1-$1):($2) w lp ls 3   t '30 M, bin = 1/8', \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 2    u (ip2-$1):($2) w l ls 666 t '240 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 2    u (ip2-$1):($2) w lp ls 66  t '240 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 2   u (ip2-$1):($2) w p ls 111 t '120 M, bin =  1  ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 2   u (ip2-$1):($2) w lp ls 11  t '120 M, bin = 1/4', \
  '../runsQ/MIx32x8a/mi_mi_Q8.txt' i 2   u (ip2-$1):($2) w lp ls 1  t '120 M, bin = 1/8'

unset logscale
unset label
set ytics auto



#--------------------------------------------------------------------------


set lmargin at screen LM4
set rmargin at screen RM4

set key b c  spacing 1.2

set ylabel 'II' offset 1,0

set ytics 0.1

unset key


plot [][-0.3:]  \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 0    u (ip1-$1):($6) w l ls 555 t '60 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 0    u (ip1-$1):($6) w lp ls 55  t '60 M, bin = 1/4', \
  '../runsQ/MIx32x8/mi_mi_Q8.txt' i 0    u (ip1-$1):($6) w lp ls 5   t '60 M, bin = 1/8', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 0    u (ip1-$1):($6) w p ls 333 t '30 M, bin =  1 ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 0    u (ip1-$1):($6) w lp ls 33  t '30 M, bin = 1/4', \
  '../runsQ/MIx32x8a/mi_mi_Q8.txt' i 0    u (ip1-$1):($6) w lp ls 3   t '30 M, bin = 1/8', \
  '../runsQ/MIx32x1/mi_mi_Q8.txt' i 2    u (ip2-$1):($6) w l ls 666 t '240 M, bin =  1  ', \
  '../runsQ/MIx32x4/mi_mi_Q8.txt' i 2    u (ip2-$1):($6) w lp ls 66  t '240 M, bin = 1/4', \
  '../runsQ/MIx32x1a/mi_mi_Q8.txt' i 2   u (ip2-$1):($6) w p ls 111 t '120 M, bin =  1  ', \
  '../runsQ/MIx32x4a/mi_mi_Q8.txt' i 2   u (ip2-$1):($6) w lp ls 11  t '120 M, bin = 1/4', \
  '../runsQ/MIx32x8a/mi_mi_Q8.txt' i 2   u (ip2-$1):($6) w lp ls 1  t '120 M, bin = 1/8'


unset logscale
unset label
set ytics auto




#==========================================================================
unset multiplot

set output



