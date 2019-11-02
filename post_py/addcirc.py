#!/usr/bin/env python
import sys, math, commands


#-----------------------------------------

def process_one(n):

    L  = 25.6;
    N  = 512;
    dx = 25.6/512;

    s =  "%04d" % n

    circfname = 'circ_data/eps1e-2_4096.rad.' + s + '.txt'
    img1fname = 'IMG_Nb03g/eps1e-2_4096.Nblur0300.' + s + '.ppm'
    img2fname = 'IMG_Nb03_circ/Nc.' + s + '.jpg'

    cfile = open(circfname, 'r')

    cmd = 'convert ' + img1fname

    for line in cfile:

        if line[0] == "#":
            continue
    
        c  = line.split()
    
        x  = float(c[1])/dx
        y  = float(c[2])/dx
        rN = float(c[4])/dx
        rH = float(c[5])/dx

        cmd = cmd + ' -fill none -stroke yellow -draw "circle '
        cmd = cmd + str(x) + ',' + str(y) + ' ' + str(x) + ',' + str(y+rH) + '" '

        cmd = cmd + ' -fill none -stroke cyan -draw "circle '
        cmd = cmd + str(x) + ',' + str(y) + ' ' + str(x) + ',' + str(y+rN) + '" '

    cfile.close();

    cmd = cmd + img2fname

    #print cmd

    output = commands.getstatusoutput(cmd)


#-----------------------------------------

nmax = 1121

for n in range(1,1+nmax,1):
    process_one(n)

#-----------------------------------------
