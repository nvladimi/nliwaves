set term pdf fsize 8 enhanced color dashed lw 1 size 8, 10

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



set output 'runsQ8_mom.pdf'

set multiplot
set size 1/2., 1/3.

set grid
unset key


# 1.i  2.Fi  3.nk    4.m2k  5.m3k  6.m4k  7.m5k  8.m6k    9.gk  10.fk    11.nnk  12.2re(ss1k)  13.2re(ss2k) 

# POST8/q7b_g3500dt5_mom.txt
# POST8/q7_v0dt1_mom.txt
# POST8/q8_v0dt1_mom.txt

set lmargin 12

#==========================================================================

set origin 0.5, 2/3.

set xlabel '2n'
set ylabel '{/Symbol D z}_{2n}'
set key t c

h0=0.36
plot [0:12][-0.05:0.3] \
   'dzeta_data.txt'  i 0 w lp ls 100 t 'V = 1,  60 modes', \
   'dzeta_data.txt'  i 1 w lp ls 200 t 'V = Fi, 60 modes', \
   -0.001*x w l ls 9

#--------------------

set origin 0.5, 1/3.

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

plot [][1e-1:1e8]  '../runsQ/POST8a/q8_v0dt1_mom.txt' \
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
     
#--------------------

set origin 0, 1/3.

set title 'V = 1,   60 modes,   {/Symbol g}_R = 1.5,   dt = 10^{-1}'

ip = 50

plot [][1e-1:1e6]  '../runsQ/POST8a/q7_v0dt1_mom.txt' \
        u (ip - $1):($4)   w lp ls 1 t '{/Symbol m}_4 ', \
   ''   u (ip - $1):($5)   w lp ls 3 t '{/Symbol m}_6 ', \
   ''   u (ip - $1):($6)   w lp ls 4 t '{/Symbol m}_8 ', \
   ''   u (ip - $1):($7)   w lp ls 6 t '{/Symbol m}_{10}', \
   ''   u (ip - $1):($8)   w lp ls 7 t '{/Symbol m}_{12}', \
   exp(x*0.017) w l ls 1, \
   exp(x*0.054) w l ls 3, \
   exp(x*0.105) w l ls 4, \
   exp(x*0.180) w l ls 6, \
   exp(x*0.280) w l ls 7
      

#--------------------

set origin 0, 0.

set xlabel 'i - i_p'

set title 'V = Fi,  60 modes,   {/Symbol g}_R = 3500,   dt = 10^{-5}'

ip = 10 
plot [][1e-1:1e4]  '../runsQ/POST8a/q7b_g3500dt5_mom.txt' \
        u ($1 - ip):($4)   w  p ls 1 t '{/Symbol m}_4 ', \
   ''   u ($1 - ip):($5)   w  p ls 3 t '{/Symbol m}_6 ', \
   ''   u ($1 - ip):($6)   w  p ls 4 t '{/Symbol m}_8 ', \
   ''   u ($1 - ip):($7)   w  p ls 6 t '{/Symbol m}_{10}', \
   ''   u ($1 - ip):($8)   w  p ls 7 t '{/Symbol m}_{12}', \
   exp(x*0.014) w l ls 1, \
   exp(x*0.041) w l ls 3, \
   exp(x*0.080) w l ls 4, \
   exp(x*0.128) w l ls 6, \
   exp(x*0.185) w l ls 7
   

#==========================================================================
unset multiplot

set output



