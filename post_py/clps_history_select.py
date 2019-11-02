#!/usr/bin/env python
import sys

infile  = open('zoombumps02.txt', 'r')
outfile = open('zoombumps08a.txt', 'w')

copy = 1

for line in infile:

    if line[0] == "#":
        
        c  = line.split()

        #  output =    0 :  hmax =   3.5655 at tmax =  10.7706 ( n =  113 )
        #0 1      2    3 4  5    6   7


        h = float(c[7])
        
        if h>8:
            outfile.write(line)
            copy = 1
        else:
            copy = 0

    else:
        if copy == 1:

            if len(line) > 40: 
                c  = line.split()
                if float(c[0]) <= 0:
                    outfile.write(line)
            else:
                outfile.write(line)


infile.close()
outfile.close()


