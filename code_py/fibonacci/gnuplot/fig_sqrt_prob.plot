set term pdf fsize 7 enhanced color dashed lw 1 size 4, 2

set datafile commentschars '%#'

PS = 1

set style line 1 lt 1 pt 4 ps PS    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 2 ps PS    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 8 ps PS    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 10 ps PS    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps PS    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps PS   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4 ps PS    lc rgb "#800080"         lw 2


set style line 11 lt 3 dt 3 pt 1 ps PS    lc rgb "#D00040"         lw 3
set style line 22 lt 3 dt 3 pt 6 ps PS    lc rgb "#FF4000"         lw 3
set style line 33 lt 3 dt 3 pt 2 ps PS    lc rgb "orange"          lw 3
set style line 44 lt 3 dt 3 pt 4 ps PS    lc rgb "#60B000"         lw 3
set style line 55 lt 3 dt 3 pt 6 ps PS    lc rgb "#008080"         lw 3
set style line 66 lt 3 dt 3 pt 12 ps PS   lc rgb  "blue"           lw 3
set style line 77 lt 3 dt 3 pt 4 ps PS    lc rgb "#800080"         lw 3

set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1


set output 'fig_sqrt_prob.pdf'

set multiplot

set xlabel font "Times-Italic"
set key font "Times-Italic"
set label font "Times-Italic"

LM1=0.08;  RM1=0.48
LM2=0.58;  RM2=0.98

BM1=0.16;  TM1=0.96

set tmargin at screen TM1
set bmargin at screen BM1

set xlabel '{/Symbol |}a_i{/Symbol |}^2 / n_i' offset 0, 0.5
set xtics 4

EV=1

#==========================================================================
set lmargin at screen LM1
set rmargin at screen RM1

set ylabel 'probability' offset 0.3,0
set logscale y
set format y "10^{%L}"
set ytics offset 0.6, 0

set key r spacing 1 samplen 2 #at 2.5, 4.8

plot [0:12][1e-5:1] \
     '../runsK/MI_16x02b/m_p0gZdt3_i06_prob.txt' u ($1):($4) ev EV w lp ls 1 t 'i - p = -24', \
     '../runsK/MI_16x02b/m_p0gZdt3_i16_prob.txt' u ($1):($4) ev EV w lp ls 2 t '-14', \
     '../runsK/MI_16x02b/m_p0gZdt3_i26_prob.txt' u ($1):($4) ev EV w lp ls 3 t '-4', \
     '../runsK/MI_16x02b/m_p0gZdt3_i36_prob.txt' u ($1):($4) ev EV w lp ls 4 t  '6', \
     '../runsK/MI_16x02b/m_p0gZdt3_i46_prob.txt' u ($1):($4) ev EV w lp ls 5 t '16', \
     '../runsK/MI_16x02b/m_p0gZdt3_i56_prob.txt' u ($1):($4) ev EV w lp ls 6 t '26', \
exp(-x) w l ls 10 notitle

unset key

unset logscale

#------------------------

set lmargin at screen LM2
set rmargin at screen RM2

set ylabel 'compensated probability' offset -0.5,0
set ytics 1 offset 0,0
set format y "%g"


f(x) = exp(x)


plot [0:12][0.5:5] \
     '../runsK/MI_16x02b/m_p0gZdt3_i06_prob.txt' u ($1):(f($1)*($4)) ev EV w lp ls 1 t 'i - i_p = -24', \
     '../runsK/MI_16x02b/m_p0gZdt3_i16_prob.txt' u ($1):(f($1)*($4)) ev EV w lp ls 2 t '-14', \
     '../runsK/MI_16x02b/m_p0gZdt3_i26_prob.txt' u ($1):(f($1)*($4)) ev EV w lp ls 3 t '-4', \
     '../runsK/MI_16x02b/m_p0gZdt3_i36_prob.txt' u ($1):(f($1)*($4)) ev EV w lp ls 4 t  '6', \
     '../runsK/MI_16x02b/m_p0gZdt3_i46_prob.txt' u ($1):(f($1)*($4)) ev EV w lp ls 5 t '16', \
     '../runsK/MI_16x02b/m_p0gZdt3_i56_prob.txt' u ($1):(f($1)*($4)) ev EV w lp ls 6 t '26'


#==========================================================================
set lmargin at screen 0.63
set rmargin at screen 0.81

set tmargin at screen 0.92
set bmargin at screen 0.56

unset ylabel
set xlabel '{/Symbol q}_i {/Symbol / p}' offset 0,1 

unset ytics
set ytics 0.1 offset 0.2, 0 font 'Helvetica, 6'
set xtics 1 offset 0, 0.3 font 'Helvetica, 6'



plot [-1:1][0.4:0.66] \
     '../runsK/MI_16x02b/m_p0gZdt3_i06_probphi.txt' u ($1/pi):($2*pi) ev EV w  l ls 1 t 'i - i_p = -24', \
     '../runsK/MI_16x02b/m_p0gZdt3_i16_probphi.txt' u ($1/pi):($2*pi) ev EV w  l ls 2 t '-14', \
     '../runsK/MI_16x02b/m_p0gZdt3_i26_probphi.txt' u ($1/pi):($2*pi) ev EV w  l ls 3 t '-4', \
     '../runsK/MI_16x02b/m_p0gZdt3_i36_probphi.txt' u ($1/pi):($2*pi) ev EV w  l ls 4 t  '6', \
     '../runsK/MI_16x02b/m_p0gZdt3_i46_probphi.txt' u ($1/pi):($2*pi) ev EV w  l ls 5 t '16', \
     '../runsK/MI_16x02b/m_p0gZdt3_i56_probphi.txt' u ($1/pi):($2*pi) ev EV w  l ls 6 t '26', \
     0.5 w l ls 0


#==========================================================================



unset multiplot

set output



