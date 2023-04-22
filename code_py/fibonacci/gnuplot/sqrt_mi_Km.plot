set term pdf fsize 7 enhanced color dashed lw 1 size 8, 5.5

set datafile commentschars '%#'



set style line 1 lt 1 pt 6  ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 4  ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 10  ps 1  lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4  ps 1    lc rgb "#800080"         lw 2


set style line 11 dt 3 pt 6  ps 1  lc rgb "#D00040"         lw 2
set style line 22 dt 3 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 3 pt 4 ps 1    lc rgb "orange"          lw 2
set style line 44 dt 3 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 55 dt 3 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 66 dt 3 pt 10  ps 1.3  lc rgb  "blue"           lw 2
set style line 77 dt 3 pt 2  ps 1    lc rgb "#800080"         lw 2





set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'sqrt_mi_Km.pdf'

set multiplot




# 1.i  2.I123   3.I12  4.I23  5.I31  6.II
#
# 1.i  2.S123   3.S12  4.S23  5.S31  6.S1  7.S2  8.S3
# 9.S1o  10.S2o  11.S3o  12.Sao  13.Sa123  14.Sa12  15.Sa23  16.Sa31
# 17.S1s  18.S2s  19.S3s


#==========================================================================


LM1=0.07;  RM1=0.33
LM2=0.40;  RM2=0.67
LM3=0.74;  RM3=0.99

BM1=0.08;  TM1=0.48
BM2=0.58;  TM2=0.98



set tmargin at screen TM2
set bmargin at screen BM2

set key samplen 2


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
set ytics auto

set xlabel "i - i_p" offset 0, 0.3
set ylabel 'S_1' offset 0,0
#set key t c  spacing 1.5 # at -5, 0.15

set xtics -22,8
set xtics -24,8
set xtics 10


unset key 

  plot [-30:30][]  \
  '../runsK/MI_48x1/mi_S_Km.txt' i 0    u ($1-30):($8-$19) w l ls 55, \
  '../runsK/MI_32x4/mi_S_Km.txt' i 0    u ($1-30):($8-$19) w l ls 5, \
  '../runsK/MI_32x8/mi_S_Km.txt' i 0    u ($1-30):($8-$19) w l ls 5, \
  '../runsK/MI_16x16/mi_S_Km.txt' i 0    u ($1-30):($8-$19) w l ls 5


      
    
unset title

#--------------------------------------------------------------------------

set lmargin at screen LM2
set rmargin at screen RM2

set key b c  spacing 1.2

set xlabel 'i - i_p' offset 0,0

set ylabel 'S_{ij} + 1.39 i' offset -1.5,0

set ytics auto

a=1.39*0

plot[-30:30][]  \
'../runsK/MI_48x1/mi_S_Km.txt' \
      u ($1-30):($3+a*$1) i 0 w l ls 66 t 'S_{12}', \
   '' u ($1-30):($4+a*$1) i 0 w l ls 44 t 'S_{23}', \
   '' u ($1-30):($5+a*$1) i 0 w l ls 33 t 'S_{31}', \
'../runsK/MI_32x4/mi_S_Km.txt' \
      u ($1-30):($3+a*$1) i 0 w l ls 6 notitle, \
   '' u ($1-30):($4+a*$1) i 0 w l ls 4 notitle, \
   '' u ($1-30):($5+a*$1) i 0 w l ls 3 notitle, \
'../runsK/MI_32x8/mi_S_Km.txt' \
      u ($1-30):($3+a*$1) i 0 w l ls 6 notitle, \
   '' u ($1-30):($4+a*$1) i 0 w l ls 4 notitle, \
   '' u ($1-30):($5+a*$1) i 0 w l ls 3 notitle, \
'../runsK/MI_16x16/mi_S_Km.txt' \
      u ($1-30):($3+a*$1) i 0 w l ls 6 t 'S_{12}', \
   '' u ($1-30):($4+a*$1) i 0 w l ls 4 t 'S_{23}', \
   '' u ($1-30):($5+a*$1) i 0 w l ls 3 t 'S_{31}'


#--------------------------------------------------------------------------

set lmargin at screen LM3
set rmargin at screen RM3

set key b c  spacing 1.2

set xlabel 'i - i_p' offset 0,0

set ylabel 'S_{123} + 2.08 i ' offset -1.5,0

set ytics auto

a=2.08*0

plot[-30:30][]  \
'../runsK/MI_48x1/mi_S_Km.txt'   u ($1-30):($2+a*$1) i 0 w l ls 11 t 'S_{123}',  \
'../runsK/MI_32x4/mi_S_Km.txt'   u ($1-30):($2+a*$1) i 0 w l ls 1 notitle, \
'../runsK/MI_32x8/mi_S_Km.txt'   u ($1-30):($2+a*$1) i 0 w l ls 1 notitle, \
'../runsK/MI_16x16/mi_S_Km.txt'  u ($1-30):($2+a*$1) i 0 w l ls 1 notitle



#==========================================================================

set tmargin at screen TM1
set bmargin at screen BM1

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
set ytics auto

set xlabel "i - i_p" offset 0, 0.3
set ylabel 'I_{123}' offset 0,0
#set key t c  spacing 1.5 # at -5, 0.15

set xtics -22,8
set xtics -24,8
set xtics 10


unset key 

  plot [-30:30][]  \
  '../runsK/MI_48x1/mi_mi_Km.txt' i 0    u ($1-30):($2) w l ls 11, \
  '../runsK/MI_32x4/mi_mi_Km.txt' i 0    u ($1-30):($2) w l ls 1, \
  '../runsK/MI_32x8/mi_mi_Km.txt' i 0    u ($1-30):($2) w l ls 1, \
  '../runsK/MI_16x16/mi_mi_Km.txt' i 0   u ($1-30):($2) w l ls 1


#--------------------------------------------------------------------------

set lmargin at screen LM2
set rmargin at screen RM2

set key t c  spacing 1.2

set xlabel 'i - i_p' offset 0,0

set ylabel 'MI' offset 0,0

set ytics 0.01



plot[-30:30][]  \
'../runsK/MI_48x1/mi_mi_Km.txt' \
      u ($1-30):3 i 0 w l ls 66 t 'I_{12}', \
   '' u ($1-30):4 i 0 w l ls 44 t 'I_{23}', \
   '' u ($1-30):5 i 0 w l ls 33 t 'I_{31}', \
'../runsK/MI_32x4/mi_mi_Km.txt' \
      u ($1-30):3 i 0 w l ls 6 notitle, \
   '' u ($1-30):4 i 0 w l ls 4 notitle, \
   '' u ($1-30):5 i 0 w l ls 3 notitle, \
'../runsK/MI_32x8/mi_mi_Km.txt' \
      u ($1-30):3 i 0 w l ls 6 notitle, \
   '' u ($1-30):4 i 0 w l ls 4 notitle, \
   '' u ($1-30):5 i 0 w l ls 3 notitle, \
'../runsK/MI_16x16/mi_mi_Km.txt' \
      u ($1-30):3 i 0 w l ls 6 notitle, \
   '' u ($1-30):4 i 0 w l ls 4 notitle, \
   '' u ($1-30):5 i 0 w l ls 3 notitle


#--------------------------------------------------------------------------

set lmargin at screen LM3
set rmargin at screen RM3

set key b c  spacing 1.2

set xlabel 'i - i_p' offset 0,0

set ylabel 'MI' offset 0,0

set ytics 0.1


plot[-30:30][]  \
'../runsK/MI_48x1/mi_mi_Km.txt' \
   u ($1-30):6 i 0 w l ls 11 t 'II',  \
'../runsK/MI_32x4/mi_mi_Km.txt' \
    u ($1-30):6 i 0 w l ls 1 notitle, \
'../runsK/MI_32x8/mi_mi_Km.txt' \
    u ($1-30):6 i 0 w l ls 1 notitle, \
'../runsK/MI_16x16/mi_mi_Km.txt' \
    u ($1-30):6 i 0 w l ls 1 notitle




#==========================================================================
unset multiplot

set output



