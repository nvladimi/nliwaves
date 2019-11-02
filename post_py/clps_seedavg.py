#!/usr/bin/env python
import sys

infile  = open('eps1e-2_4096.clps.14', 'r')
outfile = open('time_rescaled.txt.14c', 'w')

copy = 0

for line in infile:

    c  = line.split()

    if len(line) < 4:
        continue

    if c[0] == "#_1.t":
        continue

    if c[0] == "NaN":
        continue

    if c[0] == "#":
        continue

    if c[0] == "#output":

        tmax = float(c[11])
        hmax = float(c[14])

        if hmax > 40:
            outfile.write("\n\n")
            outfile.write(line)
            copy = 1
        else:
            copy = 0

    else:
        if copy == 1:
            t = (float(c[0]) - tmax)*hmax*hmax
            h = float(c[1])/hmax

            s = " %12.6e  %12.6e\n" % (t,h)

            outfile.write(s)

infile.close()
outfile.close()


