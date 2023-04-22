set term pdf fsize 8 enhanced color dashed lw 1 size 8, 9

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



set output 'v_eq_sqrtFi.pdf'

set multiplot
set size 1/2., 1/3.

#set logscale y
set grid

set key spacing 1.2

set xlabel 'i'
set ylabel 'n_i F_i'


# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 



#==========================================================================

set origin 0, 2/3.
set ylabel 'n_i F_i  {/Symbol P}_@p^{-2/3}'

set key

#set logscale y

p0=(1.00*6765)**(-2/3.)
p1=(0.10*6765)**(-2/3.)
p2=(0.01*6765)**(-2/3.)

# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP 

set yrange [0:6]
#set yrange [0.1:10]


plot [][] \
  '../runsK/POST1/k_p0gBdt3_favgT.txt'  u ($1):($3*p0)   w lp ls 1 t  "{/Symbol g} = 3,  P = 1    ", \
  '../runsK/POST1/k_p1gBdt2_favgT.txt'  u ($1):($3*p1)   w lp ls 6 t  "{/Symbol g} = 3,  P = 0.1  ", \
  '../runsK/POST2/k_p2gBdt1_favgT.txt'  u ($1):($3*p2)   w lp ls 4 t  "{/Symbol g} = 3,  P = 0.01"

#  '../runsK/POST3/k_p1gBdt3_spc.txt'    u ($1):($3*p1)   w lp ls 5 t  "{/Symbol g} = 3,  P = 0.1  ", \

#  '' u ($1):(22*($2)**(1/3.)) w l ls 10 t '22 F_i^{1/3}' 


set origin 0.5, 2/3.

plot [][] \
  '../runsK/POST1/k_p2gDdt1_favgT.txt'  u ($1):($3*p2)   w lp ls 1 t  "{/Symbol g} = 0.3,  P = 0.01", \
  '../runsK/POST1/k_p2gCdt1_favgT.txt'  u ($1):($3*p2)   w lp ls 6 t  "{/Symbol g} = 1,    P = 0.01", \
  '../runsK/POST2/k_p2gBdt1_favgT.txt'  u ($1):($3*p2)   w lp ls 4 t  "{/Symbol g} = 3,     P = 0.01"
 
# '../runsK/POST3/k_p2gAdt2_spc.txt'    u ($1):($3*p2)   w lp ls 7 t  "{/Symbol g} = 1,    P = 0.01"


unset logscale

#==========================================================================


set origin 0, 1/3.
set ylabel '{/Symbol P} / {/Symbol P}_p '
#set yrange [-2:0]

set key l


p0=(1.00*6765)
p1=(0.10*6765)
p2=(0.01*6765)


plot [][-0.6:0.6] \
  '../runsK/POST1/k_p0gBdt3_favgT.txt'  u ($1):($8/p0)   w lp ls 1 t  "{/Symbol g} = 3,  P = 1    ", \
  '../runsK/POST1/k_p1gBdt2_favgT.txt'  u ($1):($8/p1)   w lp ls 6 t  "{/Symbol g} = 3,  P = 0.1  ", \
  '../runsK/POST2/k_p2gBdt1_favgT.txt'  u ($1):($8/p2)   w lp ls 4 t  "{/Symbol g} = 3,  P = 0.01"

# '../runsK/POST3/k_p1gBdt3_spc.txt'    u ($1):($8/p1)   w lp ls 5 t  "{/Symbol g} = 3,  P = 0.1  ", \


set origin 0.5, 1/3.

plot [][-0.6:0.6] \
  '../runsK/POST1/k_p2gDdt1_favgT.txt'  u ($1):($8/p2)   w lp ls 1 t  "{/Symbol g} = 0.3,  P = 0.01", \
  '../runsK/POST1/k_p2gCdt1_favgT.txt'  u ($1):($8/p2)   w lp ls 6 t  "{/Symbol g} = 1,    P = 0.01", \
  '../runsK/POST2/k_p2gBdt1_favgT.txt'  u ($1):($8/p2)   w lp ls 4 t  "{/Symbol g} = 3,    P = 0.01", \
  '../runsK/POST3/k_p2gAdt2_spc.txt'    u ($1):($8/p2)   w lp ls 7 t  "{/Symbol g} = 1,    P = 0.01"


#==========================================================================

set ylabel '| {/Symbol x} |'

# ksi = Pi (V F)^{-1} n^{-3/2}

#set ytics 0.1

p=-1.5
q=-1.5


set origin 0, 0.
set key r at 35, 0.95

plot [][0:1] \
  '../runsK/POST1/k_p0gBdt3_favgT.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 1 t  "{/Symbol g} = 3,  P = 1    ", \
  '../runsK/POST1/k_p1gBdt2_favgT.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 6 t  "{/Symbol g} = 3,  P = 0.1  ", \
  '../runsK/POST2/k_p2gBdt1_favgT.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 4 t  "{/Symbol g} = 3,  P = 0.01"

# '../runsK/POST3/k_p1gBdt3_spc.txt'    u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 5 t  "{/Symbol g} = 3,  P = 0.1  ", \


set origin 0.5, 0.


plot [][0:1] \
  '../runsK/POST1/k_p2gDdt1_favgT.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 1 t  "{/Symbol g} = 0.3,  P = 0.01", \
  '../runsK/POST1/k_p2gCdt1_favgT.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 6 t  "{/Symbol g} = 1,    P = 0.01", \
  '../runsK/POST2/k_p2gBdt1_favgT.txt'  u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 4 t  "{/Symbol g} = 3,     P = 0.01"

#  '../runsK/POST3/k_p2gAdt2_spc.txt'    u ($1):(abs($8)*($3/$2)**p*($2)**q)   w lp ls 7 t  "{/Symbol g} = 1,    P = 0.01"




#==========================================================================



unset multiplot

set output



