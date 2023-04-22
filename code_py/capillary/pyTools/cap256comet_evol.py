import sys

fbase      = sys.argv[1]
fname_out =  sys.argv[2]
istart    =  int(sys.argv[3])    # first file to process
iend      =  int(sys.argv[4])    # stop before reading it (even non-existing)
every     =  int(sys.argv[5])    # process one in every so many slices
dt        =  float(sys.argv[6])  # time between snapshots


No     = 256         # size of original simulation
Nklo   =  80         # size of reduced data

fbase     = sys.argv[1]
fbaseout =  sys.argv[2]


#-------------------------

import numpy as np
import os
import nlwTools as nlw


nslices = 0
nfiles  = 0

for ifile in np.arange(istart, iend):

    fname = fbase + '.klo.' + str(ifile).zfill(4)
    #print(fname)

    if os.path.isfile(fname) :

        klo = nlw.ReadKlo(fname, Nklo, No, spectrum = False)
        t0 = (ifile - 1) * dt * klo.shape[0]

        nlw.KloEvolution(klo, t0, dt, every, fname_out)

        nslices += klo.shape[0]
        nfiles += 1
        msg = "Last file read: {}".format(fname)


print("{}. Number of slices in {} files read: {}".format(msg, nfiles, nslices))


#-------------------------
