set term pdf fsize 8 enhanced color dashed lw 1 size 8, 8

set datafile commentschars '%#'


set style line 1 lt 1 pt 1 ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 6 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 2 ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4 ps 1    lc rgb "#800080"         lw 2

set style line 10   dt 2 pt 7 ps 0.5  lc rgb "black"        lw 2
set style line 100  lt 3 pt 6 ps 2    lc rgb "#D00040"      lw 2
set style line 200  lt 3 pt 6 ps 2    lc rgb "blue"         lw 2
set style line 101  dt 3 pt 6 ps 2    lc rgb "#D00040"      lw 2
set style line 202  dt 3 pt 6 ps 2    lc rgb "blue"         lw 2
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"        lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1

set style line 101  dt 2 pt 6 ps 2    lc rgb "#D00040"      lw 2
set style line 202  dt 2 pt 6 ps 2    lc rgb "blue"         lw 2


set output 'runsQ8_mts.pdf'

set multiplot
set size 1/2., 1/2.

set grid
unset key


set lmargin 12

#==========================================================================

set origin 0, 1/2.

set xlabel 'p'
set ylabel '{/Symbol D z}_{p}'
set key t l

d3a=1.01563649285418
d3b=1.98765

plot [0:12][] \
   'dzeta_data.txt'  i 0 w lp ls 101 t 'old', \
   'dzeta_data.txt'  i 1 w lp ls 202 t 'old', \
   'dzeta_data.txt' u ($1):($1*d3a/3. + $3)  i 2 w lp ls 100 t 'V = 1,  60 modes', \
   'dzeta_data.txt' u ($1):($1*d3b/3. + $3)  i 3 w lp ls 200 t 'V = Fi, 60 modes', \
   -0.001*x w l ls 9

#--------------------

set origin 0.5, 1/2.

set xlabel 'i'
set ylabel 'S_i' offset -2,0
set key b c spacing 1.5
s0 = 0.058648


plot [][-0.10:0.02] 'mi_data.txt' \
     u ($1):($6-s0)    i  2  w  p ls 100 t 'V = 1,   60 modes', \
  '' u ($1-1):($5-s0)  i  2  w  p ls 100 notitle, \
  '' u ($1-2):($4-s0)  i  2  w  p ls 100 notitle, \
  '' u ($1):($6-s0)    i 10  w  p ls 200 t 'V = Fi,  60 modes', \
  '' u ($1-1):($5-s0)  i 10  w  p ls 200 notitle, \
  '' u ($1-2):($4-s0)  i 10  w  p ls 200 notitle, \
  -0.001*(50-x) w l ls 101, \
  -0.001*(x-10) w l ls 202


     
#--------------------

set origin 0, 0

set xlabel 'ln F_i'
set ylabel 'ln < |a_i|^p >'
set key l b spacing 1.2

set label 1 'V = 1,   60 modes' at 16,30

ip = 50

plot [-2:30][-120:40]  '../runsQ/POSTmts/q7_v0dt1_mts.txt' \
        u (log($2)):(log($3))    w  p ls 1 t 'p_1 ', \
   ''   u (log($2)):(log($4))    w  p ls 3 t 'p_2 ', \
   ''   u (log($2)):(log($5))    w  p ls 4 t 'p_3 ', \
   ''   u (log($2)):(log($6))    w  p ls 5 t 'p_4 ', \
   ''   u (log($2)):(log($7))    w  p ls 6 t 'p_5 ', \
   ''   u (log($2)):(log($8))    w  p ls 7 t 'p_6 ', \
   ''   u (log($2)):(log($9))    w  p ls 1 t 'p_7 ', \
   ''   u (log($2)):(log($10))   w  p ls 3 t 'p_8 ', \
   ''   u (log($2)):(log($11))   w  p ls 4 t 'p_9 ', \
   ''   u (log($2)):(log($12))   w  p ls 5 t 'p_{10} ', \
   ''   u (log($2)):(log($13))   w  p ls 6 t 'p_{11} ', \
   ''   u (log($2)):(log($14))   w  p ls 7 t 'p_{12} ', \
 -7.96136625519147 -0.32608976608663  *x  w l ls 1 notitle, \
 -6.2083320493598  -0.66442813076154  *x  w l ls 3 notitle, \
 -4.02166682126255 -1.01563649285418  *x  w l ls 4 notitle, \
 -1.41460541034588 -1.38159388396599  *x  w l ls 5 notitle, \
  1.66331888426289 -1.7672513542457   *x  w l ls 6 notitle, \
  5.30039039591479 -2.17919497707706  *x  w l ls 7 notitle, \
  9.5047450050387  -2.61794603796502  *x  w l ls 1 notitle, \
 14.1538188513629  -3.07543260051455  *x  w l ls 3 notitle, \
 19.0869153309202  -3.54188354945607  *x  w l ls 4 notitle, \
 24.1812260850932  -4.01063447684968  *x  w l ls 5 notitle, \
 29.3618957291013  -4.47824389844392  *x  w l ls 6 notitle, \
 34.5873828294463  -4.94325399177845  *x  w l ls 7 notitle


#--------------------

set origin 1/2., 0.

set xlabel 'ln F_i'
set ylabel 'ln < |a_i|^p >'

set label 1 'V = Fi,  60 modes' at 16,20

set key spacing 1.2

ip = 10 

plot [-2:30][-200:25]  '../runsQ/POSTmts/q7b_g3500dt5_mts.txt' \
        u (log($2)):(log($3))    w  p ls 1 t 'p_1 ', \
   ''   u (log($2)):(log($4))    w  p ls 3 t 'p_2 ', \
   ''   u (log($2)):(log($5))    w  p ls 4 t 'p_3 ', \
   ''   u (log($2)):(log($6))    w  p ls 5 t 'p_4 ', \
   ''   u (log($2)):(log($7))    w  p ls 6 t 'p_5 ', \
   ''   u (log($2)):(log($8))    w  p ls 7 t 'p_6 ', \
   ''   u (log($2)):(log($9))    w  p ls 1 t 'p_7 ', \
   ''   u (log($2)):(log($10))   w  p ls 3 t 'p_8 ', \
   ''   u (log($2)):(log($11))   w  p ls 4 t 'p_9 ', \
   ''   u (log($2)):(log($12))   w  p ls 5 t 'p_{10} ', \
   ''   u (log($2)):(log($13))   w  p ls 6 t 'p_{11} ', \
   ''   u (log($2)):(log($14))   w  p ls 7 t 'p_{12} ', \
   -7.34488 - 0.673823*x w l ls 1 notitle, \
   -5.3157  - 1.33604*x  w l ls 3 notitle, \
   -3.19166 - 1.98765*x  w l ls 4 notitle, \
   -1.00679 - 2.6294 *x  w l ls 5 notitle, \
    1.22076 - 3.26202*x  w l ls 6 notitle, \
    3.48181 - 3.88629*x  w l ls 7 notitle, \
    5.7723  - 4.50303*x  w l ls 1 notitle, \
    8.09084 - 5.113  *x  w l ls 3 notitle, \
    10.4381 - 5.717  *x  w l ls 4 notitle, \
    12.8182 - 6.31603*x  w l ls 5 notitle, \
    15.2393 - 6.9114 *x  w l ls 6 notitle, \
    17.7095 - 7.50448*x  w l ls 7 notitle



#==========================================================================
set xlabel 'i_p - i'
set ylabel '{/Symbol m}_{2n}'

set logscale y
set format y "10^{%L}"

#--------------------

set origin 0, 2/3.

set title 'V = 1,  80 modes,   {/Symbol g}_R = 1.5,   dt = 10^{-1}'

ip = 70

set key t l

#plot [][1e-1:1e8]  '../runsQ/POST8a/q8_v0dt1_mom.txt' \
        u (ip - $1):($4)   w lp ls 1 t '{/Symbol m}_4 ', \
   ''   u (ip - $1):($5)   w lp ls 3 t '{/Symbol m}_6 ', \
   ''   u (ip - $1):($6)   w lp ls 4 t '{/Symbol m}_8 ', \
   ''   u (ip - $1):($7)   w lp ls 6 t '{/Symbol m}_{10}', \
   ''   u (ip - $1):($8)   w lp ls 7 t '{/Symbol m}_{12}', \
   exp(x*0.017) w l ls 1, \
   exp(x*0.060) w l ls 3, \
   exp(x*0.200) w l ls 4, \
   exp(x*0.300) w l ls 6, \
   exp(x*0.400) w l ls 7


#==========================================================================
unset multiplot

set output



