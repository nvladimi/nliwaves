set term pdf fsize 8 enhanced color dashed lw 1 size 8, 7

set datafile commentschars '%#'


set style line 1 lt 1 pt 1  ps 1    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 1  ps 1    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 2  ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 8  ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 10 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1    lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4  ps 1    lc rgb "#800080"         lw 2


set style line 10   dt 2 pt 7 ps 0.5   lc rgb "black"       lw 2
set style line 100  lt 1 pt 7 ps 0.5   lc rgb "black"       lw 3
set style line 1000 lt 4 pt 6 ps 1     lc rgb "black"       lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 1
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1


set output 'mi_vs_t.pdf'

set multiplot

set size 0.5, 0.5

#set logscale y
set grid
unset key

#1.i0   2.MI      3.S4D       4.S1     5.S2     6.S3     7.S4

s0 = 0.058648  # entropy shift due to finite bin size (dn = 1)
#s0=0


set xlabel '| i  - i_p |'
set ylabel 'MI'
set xtics -26,4

set key l

#set yrange [0:0.35]

#==========================================================================

set origin 0,0.5

ip=50

set title "V = 1,  {/Symbol g}_L = 1.5,  pumping at 50 out ot 60 

plot  [0:44]  'mi_vs_t_data.txt'  \
       u (ip-$1):($2) i 0  w lp ls 1  t 't = [0.0, 0.6] x 10^6', \
    '' u (ip-$1):($2) i 1  w lp ls 3  t 't = [0.6, 1.2] x 10^6', \
    '' u (ip-$1):($2) i 2  w lp ls 4  t 't = [1.2, 1.8] x 10^6', \
    '' u (ip-$1):($2) i 3  w lp ls 5  t 't = [1.8, 2.4] x 10^6', \
    '' u (ip-$1):($2) i 4  w lp ls 6  t 't = [2.4, 3.0] x 10^6', \
     'mi_data.txt'\
       u (ip-$1):($2)  i 2   w lp ls 100  t 'overall', \
    ''  u (70-$1):($2)  i 16   w lp ls 7  t '80 modes'

    
 #   '' u (ip-$1):($2) i 14 w lp ls 100  t 'overall'


unset logscale
unset label
set ytics auto

#==========================================================================


set origin 0.5, 0.5

set title "V = Fi,  {/Symbol g}_R = 3500,  pumping at 10 out ot 60 

set key spacing 1.2

  ip=10

plot  [4:] 'mi_vs_t_data.txt'  \
       u ($1-ip):($2) i 5  w lp ls 1  t  't =     [0, 100]', \
    '' u ($1-ip):($2) i 6  w lp ls 3  t  't = [100, 200]', \
    '' u ($1-ip):($2) i 7  w lp ls 4  t  't = [200, 300]', \
    '' u ($1-ip):($2) i 8  w lp ls 5  t  't = [300, 400]', \
    '' u ($1-ip):($2) i 9  w lp ls 6  t  't = [400, 500]', \
        'mi_data.txt' \
       u ($1-ip):($2)  i 10   w lp ls 100  t 'overall'

#  '' u ($1-ip):($2) i 10  w lp ls 100  t  'overall'
 
#==========================================================================




set origin 0,0

ip=50

set title "V = 1,   {/Symbol g}_L = 1.5,  pumping at 50 out ot 60 

plot [0:44] 'mi_vs_t_data.txt'  \
       u (ip-$1):($2) i 17  w lp ls 1  t 't = [2.8, 3.2] x 10^6', \
    '' u (ip-$1):($2) i 16  w lp ls 3  t 't = [2.4, 3.2] x 10^6', \
    '' u (ip-$1):($2) i 15  w lp ls 4  t 't = [1.6, 3.2] x 10^6', \
    '' u (ip-$1):($2) i 14  w lp ls 6  t 't = [   0, 3.2] x 10^6', \
        'mi_data.txt' \
       u (ip-$1):($2) i 2  w lp ls 100  t  'overall', \
    ''  u (70-$1):($2)  i 16   w lp ls 7  t '80 modes'

 
unset logscale
unset label
set ytics auto

#==========================================================================


set origin 0.5, 0

set title "V = Fi,  {/Symbol g}_R = 3500,  pumping at 10 out ot 60 

set key spacing 1.2

  ip=10

plot  [4:] 'mi_vs_t_data.txt'  \
       u ($1-ip):($2) i 13  w lp ls 1  t  't = [444, 512]', \
    '' u ($1-ip):($2) i 12  w lp ls 3  t  't = [384, 512]', \
    '' u ($1-ip):($2) i 11  w lp ls 4  t  't = [256, 512]', \
    '' u ($1-ip):($2) i 10  w lp ls 6  t  't = [    0, 512]', \
        'mi_data.txt'\
       u ($1-ip):($2)  i 10   w lp ls 100  t 'overall'
 




#==========================================================================
unset multiplot

set output



