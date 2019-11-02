#!/usr/bin/env python
import sys, commands

#-----------------------------------------

def process_file(fname):

    infile  = open(fname, 'r')
    outfile = open(fname + '.sort1', 'a')

    for line in infile:

        c  = line.split()

        if len(line) <= 1:
            continue

        if c[0] == "NaN":
            continue

        if (c[0][0] == "#" and  not c[0] == "#collapse"):
            continue

        if c[0] == "#collapse":

            outfile.close()

            hmax = float(c[10])

            if hmax > 40:
                outfile = open(fname + ".sort4", "a")
            elif hmax > 20:
                outfile = open(fname + ".sort3", "a")
            elif hmax > 10:
                outfile = open(fname + ".sort2", "a")
            else:
                outfile = open(fname + ".sort1", "a")
        
            outfile.write("\n\n")
            outfile.write(line)
            
        else:

            outfile.write(line)

    infile.close()
    outfile.close()


#-----------------------------------------

for n in range(1,15,1):

    s = "%02d" % n
    fname = 'seed1b.post.' + s

    process_file(fname)

#-----------------------------------------
