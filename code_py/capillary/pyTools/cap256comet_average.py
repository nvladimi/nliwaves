import sys

fbase     =  sys.argv[1]
fbaseout  =  sys.argv[2] 
istart    =  int(sys.argv[3])  # first file to process
iend      =  int(sys.argv[4])  # stop before reading it (even non-existing)

dNk    = 1     # spacing in |k| for polar grid
Nphi   = 120   # number of angles in polar grid

No     = 256   # size of original simulation
Nklo   = 80    # size of reduced data

#----------------------

import os
import numpy as np
import nlwTools as nlw

nfiles = 0
nslices  = 0

nlw.KloAverageClear()

for ifile in np.arange(istart, iend):

    fname = fbase + '.klo.' + str(ifile).zfill(4)
    #print(fname)

    if os.path.isfile(fname):

        klo = nlw.ReadKlo(fname, Nklo, No, spectrum = False)

        nlw.KloAverageAdd(klo, dNk=dNk, Nphi=Nphi)

        nslices += klo.shape[0]
        nfiles += 1             
        msg = "Last file read: {}".format(fname)

print("{}. Total {} slices in {} files".format(msg, nslices, nfiles))

a, b, count = nlw.KloAverageGet()

np.save(fbaseout+"_kxky.npy", a)
np.save(fbaseout+"_kphi.npy", b)

print("Number of slices in {} files read: {} ".format(nfiles,  nslices) )

#------------------------
