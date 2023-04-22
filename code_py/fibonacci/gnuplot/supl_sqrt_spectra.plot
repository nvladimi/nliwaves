set term pdf fsize 8 enhanced color dashed lw 1 size 8, 4

set datafile commentschars '%#'


set style line 1 lt 1 pt 1 ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 6 ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 2 ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4 ps 1    lc rgb "#800080"         lw 2


set style line 11 lt 3 dt 3 pt 1 ps 1    lc rgb "#D00040"         lw 3
set style line 22 lt 3 dt 3 pt 6 ps 1    lc rgb "#FF4000"         lw 3
set style line 33 lt 3 dt 3 pt 2 ps 1    lc rgb "orange"          lw 3
set style line 44 lt 3 dt 3 pt 4 ps 1    lc rgb "#60B000"         lw 3
set style line 55 lt 3 dt 3 pt 6 ps 1    lc rgb "#008080"         lw 3
set style line 66 lt 3 dt 3 pt 12 ps 1   lc rgb  "blue"           lw 3
set style line 77 lt 3 dt 3 pt 4 ps 1    lc rgb "#800080"         lw 3

set style line 10   dt 2 pt 7 ps 0.5    lc rgb "black"       lw 2
set style line 100  lt 3 pt 6 ps 1    lc rgb "black"       lw 1
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1



set output 'supl_sqrt_spectra.pdf'

set multiplot
set size 1/2., 1

#set logscale y
set grid

set key spacing 1.2



# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 


Fia=6765
Fib=832040

ia=20
ib=30

#==========================================================================

set origin 0, 0

set ylabel 'n_i F_i  {/Symbol P}_@p^{-2/3}'
set xlabel 'i - i_p'



set key

#set logscale y

Fi=6765
p0a=(1.00*Fia)**(-2/3.)
p1a=(0.10*Fia)**(-2/3.)
p2a=(0.01*Fia)**(-2/3.)

p0b=(1.00*Fib)**(-2/3.)
p1b=(0.10*Fib)**(-2/3.)
p2b=(0.01*Fib)**(-2/3.)


plot [][0:6] \
  '../runsK/POST1/k_p1gBdt2_favgT.txt'  u ($1-ia):($3*p1a)   w lp ls 1 t  "{/Symbol g} = 3,  P = 0.1  ", \
  '../runsK/POST3/k_p2gAdt2_spc.txt'    u ($1-ia):($3*p2a)   w lp ls 3 t  "{/Symbol g} = 1,  P = 0.01", \
  '../runsK/POST2/m_p0gZdt3_favgT.txt'  u ($1-ib):($3*p0b)   w lp ls 4 t  "{/Symbol g} = 30,  P = 1     ", \
  '../runsK/POST3/m_p1gAdt3_spc.txt'    u ($1-ib):($3*p1b)   w lp ls 6 t  "{/Symbol g} = 10,  P = 0.1  "


#==========================================================================

set origin 0.5, 0

set lmargin 8
set bmargin 4

set ylabel '| {/Symbol x} |' offset 1,0

# ksi = Pi (V F)^{-1} n^{-3/2}

#set ytics 0.1

p=-1.5
q=-1.5


set key r #at 35, 0.95

set logscale

set xtics (1, 2, 3,  5, 10, 15, 20, 25, 30)
set ytics (0.01, 0.02, 0.05, 0.1, 0.2, 0.5,  1, 2, 5)
set xlabel 'i-1,   m-i'

m1=1; m2a=40; m2b=60;

plot [1:30][0.03:2] \
  '../runsK/POST1/k_p1gBdt2_favgT.txt'  u ($1-m1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 1 t  "{/Symbol g} = 3,  P = 0.1  ", \
  '../runsK/POST3/k_p2gAdt2_spc.txt'    u ($1-m1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 3 t  "{/Symbol g} = 1,  P = 0.01", \
  '../runsK/POST1/k_p1gBdt2_favgT.txt'  u (m2a-$1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 11  notitle, \
  '../runsK/POST3/k_p2gAdt2_spc.txt'    u (m2a-$1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 33  notitle, \
  '../runsK/POST2/m_p0gZdt3_favgT.txt'  u ($1-m1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 4 t  "{/Symbol g} = 30,  P = 1     ", \
  '../runsK/POST3/m_p1gAdt3_spc.txt'    u ($1-m1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 6 t  "{/Symbol g} = 10,  P = 0.1  ", \
  '../runsK/POST2/m_p0gZdt3_favgT.txt'  u (m2b-$1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 44 notitle, \
  '../runsK/POST3/m_p1gAdt3_spc.txt'    u (m2b-$1):(abs($8)*($3/$2)**p*($2)**q)  w lp ls 66 notitle, \
 1.5/x w l ls 10 t '1.5 / i'



#==========================================================================



unset multiplot

set output



