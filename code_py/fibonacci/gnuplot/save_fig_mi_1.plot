set term pdf fsize 8 enhanced color dashed lw 1 size 8, 6.5

set datafile commentschars '%#'


set style line 1 lt 1 pt 6  ps 1.5  lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 2  ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4  ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 4  ps 1.3  lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4  ps 1    lc rgb "#800080"         lw 2


set style line 11 dt 3 pt 6  ps 1.5  lc rgb "#D00040"         lw 2
set style line 22 dt 3 pt 12 ps 1    lc rgb "#FF4000"         lw 2
set style line 33 dt 3 pt 2  ps 1    lc rgb "orange"          lw 2
set style line 44 dt 3 pt 4  ps 1    lc rgb "#60B000"         lw 2
set style line 55 dt 3 pt 6  ps 1    lc rgb "#008080"         lw 2
set style line 66 dt 3 pt 4  ps 1.3  lc rgb  "blue"           lw 2
set style line 77 dt 3 pt 2  ps 1    lc rgb "#800080"         lw 2




set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'fig_mi_save.pdf'

set multiplot


#set logscale y
set grid
unset key

#1.i0   2.MI      3.S4D       4.S1     5.S2     6.S3     7.S4

s0 = 0.058648  # entropy shift due to finite bin size (dn = 1)
#s0=0

set xtics -26,4

#==========================================================================

set lmargin at screen 0.08
set rmargin at screen 0.40
set tmargin at screen 0.94
set bmargin at screen 0.56

ip1=30    #ip1=50 for early time
ip2=50
ip3=70

I1=0      #I1=4 for early time
I2=1      #I1=5 for early time
I3=2
I4=3

set ylabel 'MI,  {/Symbol D}S' 
set xlabel 'i_p - i'

set key outside r spacing 1.2 # at -5, 0.15

set title "V = 1, {/Symbol P} = 67.65, {/Symbol g}_L = 1.5: pumping at 30 out of 40 (dashed) or at 50 out ot 60 (solid)" offset 25, 0

  plot [0:]'mi_data.txt'  \
       u (ip1-$1):($2)      i I1  w lp ls 11   notitle,\
    '' u (ip1-$1):($2)      i I2  w lp ls 66   notitle,\
    '' u (ip1-$1):($3-3*s0) i I1  w lp ls 55   notitle,\
    '' u (ip1-$1):($3-2*s0) i I2  w lp ls 44   notitle,\
    '' u (ip1-$1):($4-s0)   i I1  w lp ls 33   notitle,\
    '' u (ip1-$1):($5-s0)   i I1  w lp ls 33   notitle,\
    '' u (ip1-$1):($6-s0)   i I1  w lp ls 33   notitle, \
    '' u (ip1-$1):($7)      i I1  w lp ls 22   notitle, \
    '' u (ip2-$1):($2)      i I3  w lp ls 1  t '3-mode MI',\
    '' u (ip2-$1):($2)      i I4  w lp ls 6  t '2-mode MI',\
    '' u (ip2-$1):($3-3*s0) i I3  w lp ls 5  t 'S(n_{n-2}, n_{i-1}, n_{i}, {/Symbol q}_3)',\
    '' u (ip2-$1):($3-2*s0) i I4  w lp ls 4  t 'S(n_{i-1}, n_{i}, {/Symbol q}_2)',\
    '' u (ip2-$1):($4-s0)   i I3  w lp ls 3  t 'S(n_{i-2})',\
    '' u (ip2-$1):($5-s0)   i I3  w lp ls 3  t 'S(n_{i-1})',\
    '' u (ip2-$1):($6-s0)   i I3  w lp ls 3  t 'S(n_i)' , \
    '' u (ip2-$1):($7)      i I3  w lp ls 2  t 'S({/Symbol q}_3)'
    
unset title
unset key

#--------------------------------------------------------------------------

set lmargin at screen 0.48
set rmargin at screen 0.80

set ytics (0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
set logscale y

set ylabel 'MI'


  plot[0:64][0.005:0.5] 'mi_data.txt'  \
       u (ip1-$1):($2) i I1  w lp ls 11 t '3-mode MI',\
    '' u (ip1-$1):($2) i I2  w lp ls 66 t '2-mode MI', \
    '' u (ip2-$1):($2) i I3  w lp ls 1  t '3-mode MI',\
    '' u (ip2-$1):($2) i I4  w lp ls 6  t '2-mode MI', \
    '' u (ip3-$1):($2) i 16  w lp ls 7  t '3-mode MI',\
    '' u (ip3-$1):($2) i 17  w lp ls 7  t '2-mode MI', \
    '' u (ip3-$1):($2) i 18  w lp ls 77  t 'old 3-mode MI',\
    '' u (ip3-$1):($2) i 19  w lp ls 77  t 'old 2-mode MI'#, \
    0.04*exp(0.07*x)  w l ls 100, \
    0.01*exp(0.08*x)  w l ls 101
    
unset logscale
unset label
set ytics auto

#==========================================================================


set lmargin at screen 0.08
set rmargin at screen 0.40
set tmargin at screen 0.44
set bmargin at screen 0.07

ip1=10
ip2=10

I1=6
I2=7
I3=10
I4=11

set ylabel 'MI,  {/Symbol D}S'
set xlabel 'i - i_p'

set key outside r # at -5, 0.15
set title "V = Fi, {/Symbol P} = 67.65: pumping at 10 out of 40 with {/Symbol g}_R = 140 (dashed) or at 10 out ot 60 with {/Symbol g}_R = 3500 (solid)" offset 30, 0


  plot [4:48]'mi_data.txt'  \
       u ($1-ip1):($2)      i I1  w lp ls 11  notitle,\
    '' u ($1-ip1):($2)      i I2  w lp ls 66  notitle,\
    '' u ($1-ip1):($3-3*s0) i I1  w lp ls 55  notitle,\
    '' u ($1-ip1):($3-2*s0) i I2  w lp ls 44  notitle,\
    '' u ($1-ip1):($4-s0)   i I1  w lp ls 33  notitle,\
    '' u ($1-ip1):($5-s0)   i I1  w lp ls 33  notitle,\
    '' u ($1-ip1):($6-s0)   i I1  w lp ls 33  notitle, \
    '' u ($1-ip1):($7)      i I1  w lp ls 22  notitle, \
    '' u ($1-ip2):($2)      i I3  w lp ls 1  t '3-mode MI',\
    '' u ($1-ip2):($2)      i I4  w lp ls 6  t '2-mode MI',\
    '' u ($1-ip2):($3-3*s0) i I3  w lp ls 5  t 'S(n_{n-2}, n_{i-1}, n_{i}, {/Symbol q}_3)',\
    '' u ($1-ip2):($3-2*s0) i I4  w lp ls 4  t 'S(n_{i-1}, n_{i}, {/Symbol q}_2)',\
    '' u ($1-ip2):($4-s0)   i I3  w lp ls 3  t 'S(n_{i-2})',\
    '' u ($1-ip2):($5-s0)   i I3  w lp ls 3  t 'S(n_{i-1})',\
    '' u ($1-ip2):($6-s0)   i I3  w lp ls 3  t 'S(n_i)' , \
    '' u ($1-ip2):($7)      i I3  w lp ls 2  t 'S({/Symbol q}_3)'
    
unset title
unset key

#--------------------------------------------------------------------------

set lmargin at screen 0.48
set rmargin at screen 0.80

set ytics (0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
set logscale y

set ylabel 'MI'


  plot[4:48][0.005:0.5] 'mi_data.txt'  \
       u ($1-ip1):($2) i I1  w lp ls 11 t '3-mode MI',\
    '' u ($1-ip1):($2) i I2  w lp ls 66 t '2-mode MI', \
    '' u ($1-ip2):($2) i I3  w lp ls 1  t '3-mode MI',\
    '' u ($1-ip2):($2) i I4  w lp ls 6  t '2-mode MI'#, \
    0.04*exp(0.07*x)  w l ls 10, \
    0.01*exp(0.08*x)  w l ls 10
    
unset logscale
unset label
set ytics auto

#==========================================================================




#==========================================================================
unset multiplot

set output



