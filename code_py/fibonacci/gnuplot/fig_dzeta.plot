set term pdf fsize 7 enhanced color dashed lw 1 size 4, 2

set datafile commentschars '%#'


PS = 1.5

set style line 1 lt 1 pt 1 ps PS    lc rgb "#D00040"         lw 2
set style line 2 lt 1 pt 6 ps PS    lc rgb "#FF4000"         lw 2
set style line 3 lt 1 pt 2 ps 1    lc rgb "orange"          lw 2
set style line 4 lt 1 pt 4 ps 1    lc rgb "#60B000"         lw 2
set style line 5 lt 1 pt 6 ps 1    lc rgb "#008080"         lw 2
set style line 6 lt 1 pt 12 ps 1   lc rgb  "blue"           lw 2
set style line 7 lt 1 pt 4 ps 1    lc rgb "#800080"         lw 2

set style line 10   dt 2 pt 7 ps 0.5  lc rgb "black"        lw 2
set style line 100  lt 3 pt 6 ps PS    lc rgb "#D00040"      lw 3
set style line 200  lt 3 pt 4 ps PS    lc rgb "blue"         lw 3
set style line 101  dt 3 pt 6 ps 2    lc rgb "#D00040"      lw 4
set style line 202  dt 3 pt 6 ps 2    lc rgb "blue"         lw 4
set style line 1000 lt 4 pt 6 ps 1    lc rgb "black"        lw 2

set style line 9   lt 5 dt 2   lc rgb "black"       lw 2
set style line 99  lt 1 dt 4   lc rgb "black"       lw 1


set output 'fig_dzeta.pdf'

set multiplot

LM1=0.06;  RM1=0.48
LM2=0.57;  RM2=0.98


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
set ylabel '{/Symbol D}' offset 2,0
set key t l spacing 1.1
set ytics 0.5
set xtics 4

d3a=1.01564
d3b=1.98765

sa =  0.0029  #      s = (d1 - d2/2)*log(1.618)
sb = -0.0028

plot [0:12][] \
   'dzeta_data.txt' u ($1):($1*d3a/3. - $3)  i 2 w lp ls 100 t '{/Symbol a} = 0', \
   'dzeta_data.txt' u ($1):($1*d3b/3. - $3)  i 3 w lp ls 200 t '{/Symbol a} = 1'
   


#--------------------

set tmargin at screen TM1
set bmargin at screen BM1
set lmargin at screen LM2
set rmargin at screen RM2

set xlabel 'i' offset 0, 1.3
set ylabel 'S_i' offset 3,0
set key b c spacing 1.1
set ytics 0.05
set xtics 20


plot [0:60][-0.25:0]   '../runsQ/MIx32x1/mi_S_Q8.txt' \
     i 0  u ($1):($8-$19) w lp ls 100 t '{/Symbol a} = 0', \
  '' i 1  u ($1):($8-$19) w lp ls 200 t '{/Symbol a} = 1', \
    0.003*(x-55) w l ls 101 t '\~ 0.003 i', \
   -0.003*(x-10) w l ls 202 t '\~ (-0.003) i'

#--------------------


set tmargin at screen 0.52
set bmargin at screen 0.17

set lmargin at screen 0.12
set rmargin at screen 0.30

unset xlabel
unset ylabel
unset key

set xtics 1 font "Helvetica, 6" offset 0, 0.4
set ytics 0.04 font "Helvetica, 6" center offset -1, 0 


plot [0:4][-0.04:0.04] \
   'dzeta_data.txt' u ($1):($1*d3a/3. - $3)  i 2 w lp ls 100 t '{/Symbol a} = 0', \
   'dzeta_data.txt' u ($1):($1*d3b/3. - $3)  i 3 w lp ls 200 t '{/Symbol a} = 1' #, \
  -0.003*x w l ls 202 t '- 0.003 p', \
   0.003*x w l ls 101 t '  0.003 p'


#==========================================================================
unset multiplot

set output



