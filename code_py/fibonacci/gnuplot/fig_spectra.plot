set term pdf fsize 7 enhanced color dashed lw 1 size 4, 4

set datafile commentschars '%#'

PS=0.7

set style line 1 lt 1 pt  4 ps PS   lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt  6 ps PS   lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt  2 ps PS   lc rgb "orange"          lw 2
set style line 4 lt 1 pt  4 ps PS   lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt  6 ps PS   lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps PS   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt  4 ps PS   lc rgb "#800080"         lw 2

set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1

set output 'fig_spectra.pdf'

set multiplot

LM1=0.07;  RM1=0.49
LM2=0.57;  RM2=0.99

BM1=0.06;  TM1=0.48
BM2=0.56;  TM2=0.98



set xlabel font "Times-Italic"
set ylabel font "Times-Italic"
set key font "Times-Italic"

set xlabel 'i'  offset 0,1

#==========================================================================

set tmargin at screen TM2
set bmargin at screen BM2

set lmargin at screen LM1
set rmargin at screen RM1

set logscale y
set format y "10^{%L}"

#set title 'V = 1,   {/Symbol g}_L = 1.5,  {/Symbol g}_R = 0,    dt = 0.1'
set ylabel 'n_i F_i' offset -0.5, 0

set key b r spacing 1.1 samplen 2

set ytics 10
set xtics 20

plot [][8:1e6] \
 '../runsQ/POSTspc/q7_spc.txt'      u ($1):(22*($2)**(1/3.)) w l ls 10 t '\~ F_i^{1/3}', \
 '../runsQ/POST5/q5_v0m20_spc.txt'  u ($1):($3)   w lp ls 6  t 'p = 20 of 40', \
 '../runsQ/POST5/q5_v0m30_spc.txt'  u ($1):($3)   w lp ls 1  t 'p = 30 of 40', \
 '../runsQ/POSTspc/q7_spc.txt'      u ($1):($3)   w lp ls 5  t 'p = 50 of 60' #, \
 '' u ($1):(22*($2)**(1/3.)) w l ls 10 t '22 F_i^{1/3}' , \

# '../runsQ/POST7/q7_v0dt1_spc.txt'  u ($1):($3)   w l  ls 3  t 'i_p = 50 of 60'

set format y "%g"


set lmargin at screen LM2
set rmargin at screen RM2

set key b l spacing 1.15

set ylabel offset 1.5,0.5

#set title 'V = Fi,  {/Symbol g}_L = 0,
#{/Symbol g}_R = 140, dt = 10^{-4} or {/Symbol g}_R = 3500,  dt = 10^{-5}'

plot [][0.003:100] \
 '../runsQ/POSTspc/q7b_spc.txt'  u ($1):(50*($2)**(-1/3.)) w l ls 10   t '\~ F_i^{-1/3}', \
 '../runsQ/POST6/q6_m20gr140_spc.txt'    u ($1):($3)   w lp ls 6  t 'p = 20 of 40', \
 '../runsQ/POST6/q6_m10gr140_spc.txt'    u ($1):($3)   w lp ls 1  t 'p = 10 of 40', \
 '../runsQ/POSTspc/q7b_spc.txt'          u ($1):($3)   w  p ls 5  t 'p = 10 of 60' #, \
'' u ($1):(50*($2)**(-1/3.)) w l ls 10   t '50 F_i^{-1/3}'         

 #'../runsQ/POST7b/q7b_g3500dt5_spc.txt'  u ($1):($3)   w l  ls 3  t 'i_p = 10 of 60'

set ytics auto

unset title
unset key
set format y "%g"


#==========================================================================


set tmargin at screen TM1
set bmargin at screen BM1


unset logscale
set format y "%g"


set lmargin at screen LM1
set rmargin at screen RM1

#set key spacing 1.15

set ylabel '{/Symbol | x |}' offset 1,0

p=-1.5
q=-1.0

set ytics 0.5

plot [][] \
  '../runsQ/POST5/q5_v0m20_spc.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  \
  w lp ls 6 t 'i_p = 20 of 40', \
  '../runsQ/POST5/q5_v0m30_spc.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)  \
  w lp ls 1 t 'i_p = 30 of 40', \
  '../runsQ/POSTspc/q7_spc.txt'      u ($1):(abs($8)*($3/$2)**p*($2)**q)  \
  w lp ls 5 t 'i_p = 50 of 60'

# '../runsQ/POST7/q7_v0dt1_spc.txt' u ($1):(abs($8)*($3/$2)**p*($2)**q) w l ls 3 t 'i_p = 50 of 60'



set lmargin at screen LM2
set rmargin at screen RM2
set ytics 0.1

p=-1.5
q=-2.0

plot [][0:0.3] \
 '../runsQ/POST6/q6_m20gr140_spc.txt'    u ($1):(abs($8)*($3/$2)**p*($2)**q) \
 w lp ls 6  t 'i_p = 20 of 40', \
 '../runsQ/POST6/q6_m10gr140_spc.txt'    u ($1):(abs($8)*($3/$2)**p*($2)**q)  \
 w lp ls 1  t 'i_p = 10 of 40', \
 '../runsQ/POSTspc/q7b_spc.txt'         u ($1):(abs($8)*($3/$2)**p*($2)**q)  \
 w lp ls 5  t 'i_p = 10 of 60'

# '../runsQ/POST7b/q7b_g3500dt5_spc.txt'       u ($1):(abs($8)*($3/$2)**p*($2)**q)  \
 w l ls 3  t 'i_p = 10 of 60'

#==========================================================================

#==========================================================================

set tmargin at screen 0.97
set bmargin at screen 0.80

unset logscale

set lmargin at screen 0.08
set rmargin at screen 0.25

set y2label '{/Symbol P / P}_p' offset -1.5,0

pp = 67.65

set xtics font "Helvetica, 6"
set ytics font "Helvetica, 6"

unset xlabel
unset ylabel

unset ytics
set y2tics (-1, -0.5, 0) font "Helvetica, 6"

#set ytics (0, 0.25, '' 0.5, 0.75, 1)

plot [][-1.2:0.2] \
  '../runsQ/POST5/q5_v0m20_spc.txt'  u ($1):($8/pp)  w lp ls 6 t 'i_p = 20 of 40', \
  '../runsQ/POST5/q5_v0m30_spc.txt'  u ($1):($8/pp)  w lp ls 1 t 'i_p = 30 of 40', \
   '../runsQ/POSTspc/q7_spc.txt'     u ($1):($8/pp)  w lp ls 5 t 'i_p = 50 of 60'

unset y2label

unset y2tics
set ytics 0.5 font "Helvetica, 6"

set lmargin at screen 0.80
set rmargin at screen 0.97

set ylabel '{/Symbol P / P}_p' offset -1,0


plot [][-0.2:1.2] \
 '../runsQ/POST6/q6_m20gr140_spc.txt'  u ($1):($8/pp)  w lp ls 6  t 'i_p = 20 of 40', \
 '../runsQ/POST6/q6_m10gr140_spc.txt'  u ($1):($8/pp)  w lp ls 1  t 'i_p = 10 of 40', \
  '../runsQ/POSTspc/q7b_spc.txt'       u ($1):($8/pp)  w lp ls 5  t 'i_p = 50 of 60'


# '../runsQ/POST7b/q7b_g3500dt5_spc.txt'  u ($1):($8)  w l ls 3  t 'i_p = 10 of 60'
 
unset key





#==========================================================================








unset multiplot

set output



