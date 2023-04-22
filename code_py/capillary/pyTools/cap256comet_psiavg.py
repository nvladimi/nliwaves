import sys

fbase     =  sys.argv[1]
fbaseout  =  sys.argv[2] 
istart    =  int(sys.argv[3])  # first file to process
iend      =  int(sys.argv[4])  # stop before reading it (even non-existing)


No     = 256   # size of original simulation

#----------------------

import os
import numpy as np
import nlwTools as nlw

nk = np.zeros((No,No))
nfiles = 0


for ifile in np.arange(istart, iend):

    fname = fbase + '.psi.' + str(ifile).zfill(4)
    #print(fname)

    if os.path.isfile(fname):

        nk += nlw.ReadPsi(fname, spectrum = True)
        nfiles += 1             
        msg = "Last file read: {}".format(fname)

print("{}. Number files read: {}".format(msg, nfiles))

nk = nk/nfiles

    
np.save(fbaseout+"_nk.npy", nk)

print("Total number of files read: {} ".format(nfiles) )

#------------------------
