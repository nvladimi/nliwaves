set term pdf fsize 8 enhanced color dashed dl 1.5 lw 1 size 5.0, 2.4

set datafile commentschars '%#'


PS = 1

set style line 1 lt 1 pt 1 ps PS    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 6 ps PS    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 2 ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4 ps 1    lc rgb "#800080"         lw 2

set style line 100  lt 3 pt 7 ps PS    lc rgb "#D00040"      lw 3
set style line 200  lt 3 pt 5 ps PS    lc rgb "blue"         lw 3
set style line 300  dt 1 pt 7 ps PS   lc rgb "orange"     lw 3
set style line 400  dt 1 pt 5 ps PS   lc rgb "#008080"     lw 3

set style line 101  dt 4 pt 6 ps PS    lc rgb "#D00040"     lw 4
set style line 202  dt 4 pt 4 ps PS   lc rgb "blue"         lw 4
set style line 303  dt 4 pt 6 ps PS   lc rgb "orange"       lw 4
set style line 404  dt 4 pt 4 ps PS   lc rgb "#008080"      lw 4

set style line 11 dt 2 pt 1 ps PS    lc rgb "#D00040"        lw 4
set style line 22 dt 2 pt 6 ps PS    lc rgb "blue"           lw 4
set style line 33 dt 2 pt 2 ps 1    lc rgb "orange"          lw 4
set style line 44 dt 2 pt 4 ps 1    lc rgb "#008080"         lw 4


set style line 9   lt 5 dt 2   lc rgb "black"       lw 2
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1


set output 'fig_dzeta_1434.pdf'

set multiplot

LM1=0.07;  RM1=0.47
LM2=0.58;  RM2=0.98


BM1=0.10;  TM1=0.97

set key samplen 4

set xlabel font "Times-Italic"
set ylabel font "Times-Italic"
set key font "Times-Italic"

#==========================================================================

set tmargin at screen TM1
set bmargin at screen BM1
set lmargin at screen LM1
set rmargin at screen RM1

set xlabel 'q' offset 0, 1.5
set ylabel '{/Symbol D}' offset 1,2
set key t l spacing 1.1
set ytics 0.2
set xtics 4

d3a = 1.01564
d3b = 1.98765
d5a = 1.00437880472105
d5b = 1.99889437107828
d14 = 1.25113641083553
d34 = 1.74823849249817

#d3a = 1
#d3b = 2
#d5a = 1
#d5b = 2
#d14 = 1.25
#d34 = 1.75

plot [0:12][] \
   'dzeta_data.txt' u ($1):($1*d5a/3. - $5)  i 2 w lp ls 100 t '{/Symbol a} = 0   ', \
   'dzeta_data.txt' u ($1):($1*d5b/3. - $5)  i 3 w lp ls 200 t '{/Symbol a} = 1   ', \
   'dzeta_data.txt' u ($1):($1*d14/3. - $3)  i 4 w lp ls 300 t '{/Symbol a} = 1/4', \
   'dzeta_data.txt' u ($1):($1*d34/3. - $3)  i 5 w lp ls 400 t '{/Symbol a} = 3/4', \
0 w l ls 0 notitle
#   'dzeta_data.txt' u ($1):($1*d3a/3. - $3)  i 2 w lp ls 101 notitle, \
#   'dzeta_data.txt' u ($1):($1*d3b/3. - $3)  i 3 w lp ls 202 notitle, \
   


#--------------------

set tmargin at screen TM1
set bmargin at screen BM1
set lmargin at screen LM2
set rmargin at screen RM2

set xlabel 'i' offset 0, 1.3
set ylabel 'S_i' offset  1,0.5    #3,0
set key b c spacing 1.1
set ytics 0.02
set xtics 20


# slope = dzeta_1 * log(1.618) = 0.01 * log(1.618) = 0.005

# alpha = 0:   dzeta_1 = 0.010      slope = 0.00481
# alpha = 1:   dzeta_1 = 0.008      slope = 0.00385
# alpha = 1/4: dzeta_1 = 0.002      slope = 0.00096
# alpha = 3/4: dzeta_1 = 0.002      slope = 0.00096

unset key


plot [0:60][-0.15:-0.04]   \
   '../runsQ/MIx32x1/mi_S_Q8.txt' \
     i 0  u ($1):($8-$19) w lp ls 100 t '{/Symbol a} = 0', \
  '' i 1  u ($1):($8-$19) w lp ls 200 t '{/Symbol a} = 1', \
     '../runsQ/POSTFrisch/mi_S_1434.txt' \
    i 0  u ($1):($8-$19) w lp ls 300 t '{/Symbol a} = 1/4', \
  '' i 1  u ($1):($8-$19) w lp ls 400 t '{/Symbol a} = 3/4', \
   -0.00385*(x-10)+0.02 w l ls 202 t '\~ 0.00385 i', \
    0.00096*(x-55)-0.025 w l ls 303 t '\~ 0.00096 i', \
   -0.00385/2*(x-30)-0.078 w l ls 22 t '\~ 0.00385 i', \
    0.00096/2*(x-30)-0.055 w l ls 33 t '\~ 0.00096 i'

#   -0.00096*(x-10)-0.025 w l ls 404 t '\~ 0.00096 i', \
#   -0.00096/2*(x-30)-0.055 w l ls 44 t '\~ 0.00096 i', \
#    0.00481*(x-55)+0.05 w l ls 101 t '\~ 0.00481 i', \
#    0.00481/2*(x-30)-0.085 w l ls 11 t '\~ 0.00481 i', \


 
#--------------------


set tmargin at screen 0.52
set bmargin at screen 0.17

set lmargin at screen 0.12
set rmargin at screen 0.30

unset xlabel
unset ylabel
unset key

set xtics 1 font "Helvetica, 7" offset 0, 0.4
set ytics 0.01 font "Helvetica, 7" center offset -1, 0 
set grid

plot [0:4][-0.02:0.02] \
   'dzeta_data.txt' u ($1):($1*d5a/3. - $5)  i 2 w lp ls 100 t '{/Symbol a} = 0', \
   'dzeta_data.txt' u ($1):($1*d5b/3. - $5)  i 3 w lp ls 200 t '{/Symbol a} = 1', \
   'dzeta_data.txt' u ($1):($1*d14/3. - $3)  i 4 w lp ls 300 t '{/Symbol a} = 1/4', \
   'dzeta_data.txt' u ($1):($1*d34/3. - $3)  i 5 w lp ls 400 t '{/Symbol a} = 3/4'
 

#   'dzeta_data.txt' u ($1):($1*d3a/3. - $3)  i 2 w lp ls 101 t '{/Symbol a} = 0', \
#   'dzeta_data.txt' u ($1):($1*d3b/3. - $3)  i 3 w lp ls 202 t '{/Symbol a} = 1', \



#==========================================================================
unset multiplot

set output



