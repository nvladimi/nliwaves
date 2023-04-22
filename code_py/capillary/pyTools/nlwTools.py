#import importlib ; import nlwTools 
#importlib.reload(nlwTools) ; nlwTools.polartest()


#==============================================================================

#==============================================================================

def polartest():
    
    import numpy as np
    import matplotlib.pyplot as plt

    N    = 40  # original grid
    Nphi = 16  # number of angles in polar grid
    dNk  =  2  # radial spacing in polar grid

    def testfun(kx, ky):
        
        #ff = 2*kx**2
        #ff = kx**2 + 0.5*(ky+5)**2
        #ff = kx**2 + ky**2
        ff = (kx - 4)**2 + (ky + 3)**2
        #ff = (kx - ky)**2 + 0.5*(kx + ky)**2
        #ff = (kx - ky)**2 + 0.5*(kx + ky - 10)**2
        
        return np.exp( -8*ff/N/N) 
            
    #-- create test data -- 

    k = np.fft.ifftshift(np.arange(N) - N/2)

    kx, ky = np.meshgrid(k,k)

    a = testfun(kx, ky)

    qkx = np.sum(a * kx)
    qky = np.sum(a * ky)
    phi0 = np.arctan2(qky, qkx)

    print("qkx, qky, phi0/pi: ", qkx, qky, phi0/np.pi)

    
    #-- setup polar mesh --
    
    phi = 2*np.pi*np.arange(Nphi)/Nphi + phi0; 
    r   = np.arange(dNk, N/2, dNk);
    Nk  = r.shape[0]
    pp,rr  = np.meshgrid(phi, r)

    x = rr * np.cos(pp)
    y = rr * np.sin(pp)
  
    #-- interpolation --

    ix1 = np.floor(x).astype(int)
    ix2 = np.floor(x+1).astype(int)
    iy1 = np.floor(y).astype(int)
    iy2 = np.floor(y+1).astype(int)
        
    b = (iy2 - y)*(ix2 - x)*a[iy1,ix1] + (iy2 - y)*(x - ix1)*a[iy1,ix2] + \
        (y - iy1)*(ix2 - x)*a[iy2,ix1] + (y - iy1)*(x - ix1)*a[iy2,ix2]
    
    #-- image --
    
    fig, ax = plt.subplots(ncols=3, nrows=1, constrained_layout=True, figsize=(12,4))
    
    ax[0].imshow(np.fft.fftshift(a), cmap='jet', origin = 'lower');
    print("original     min, max: ", np.min(a), np.max(a))

    ax[1].imshow(b, cmap='jet', origin = 'lower');
    print("interpolated min, max: ", np.min(b), np.max(b))
    
    z = np.linspace(-N/2, N/2, 200)
    ax[2].plot(k, a[:,0], 'or', markersize = 3)
    ax[2].plot(r, b[0:,0], 'ob', markersize = 3)
    ax[2].plot(z, testfun(z,0), '-k', linewidth = 0.5)

    
    plt.show()


    
    
#==============================================================================

def KloAverageAdd(klo, dNk=4, Nphi=64):

    import numpy as np

    global klo_avg_aak_orig
    global klo_avg_aak_polar
    global klo_avg_count

    #-- initiate domains --
    
    n, Nx, Ny = klo.shape

    if Nx == Ny:
        N = Nx
        Nk = int(N/2)
    else:
        print("Non-square domain")
        return(-1)

    #-- setup polar mesh --
   
    phi = 2*np.pi*np.arange(Nphi)/Nphi; 
    r   = np.arange(dNk, N/2, dNk);
    Nk  = r.shape[0]

    #-- setup global data arrays to store averages --
    
    if "klo_avg_count" in globals():
        pass
    else:
        klo_avg_count = 0
        klo_avg_aak_orig  = np.zeros((N, N))
        klo_avg_aak_polar = np.zeros((Nk, Nphi))
        
        
    #-- setup cartesian mesh --
        
    k = np.fft.ifftshift(np.arange(N) - N/2)

    kx, ky = np.meshgrid(k,k)
    kk = kx*kx + ky*ky
    qx = kx * kk * np.sqrt(kk)
    qy = ky * kk * np.sqrt(kk)


    #-- process all slices in "klo" array --
      
    for i in np.arange(0, n):

        #-- get slice --
        
        ak   = klo[i, :, :]
        aak  = np.real( ak * np.conjugate(ak) )       
        a    = aak * kk * kk

        #-- compute anisotropy vector --
        
        qkx  = np.sum(aak * qx)
        qky  = np.sum(aak * qy)
        phi0 = np.arctan2(qky, qkx)

        #-- shift polar mesh --
        
        pp,rr  = np.meshgrid(phi+phi0, r)
       
        x = rr * np.cos(pp)
        y = rr * np.sin(pp)

        #-- interpolate --

        ix1 = np.floor(x).astype(int)
        ix2 = np.floor(x+1).astype(int)
        iy1 = np.floor(y).astype(int)
        iy2 = np.floor(y+1).astype(int)

        b = (iy2 - y)*(ix2 - x)*a[iy1,ix1] + (iy2 - y)*(x - ix1)*a[iy1,ix2] + \
            (y - iy1)*(ix2 - x)*a[iy2,ix1] + (y - iy1)*(x - ix1)*a[iy2,ix2]

        #-- add to averages --
        
        klo_avg_aak_orig  += a
        klo_avg_aak_polar += b
        klo_avg_count += 1
    
    return

    
#==============================================================================

def nkInterpolate(aa, dNk=1, Nphi=120, AngleAvg = True, iShift = False):

    import numpy as np

    N = aa.shape[0]

    if iShift == True:
        aa = np.fft.ifftshift(aa)
    
    #-- setup polar mesh --
   
    phi = 2*np.pi*np.arange(Nphi)/Nphi; 
    r   = np.arange(dNk, N/2, dNk);
  
    pp,rr  = np.meshgrid(phi, r)
       
    x = rr * np.cos(pp)
    y = rr * np.sin(pp)

    #-- interpolate --

    ix1 = np.floor(x).astype(int)
    ix2 = np.floor(x+1).astype(int)
    iy1 = np.floor(y).astype(int)
    iy2 = np.floor(y+1).astype(int)

    bb = (iy2 - y)*(ix2 - x)*aa[iy1,ix1] + (iy2 - y)*(x - ix1)*aa[iy1,ix2] + \
         (y - iy1)*(ix2 - x)*aa[iy2,ix1] + (y - iy1)*(x - ix1)*aa[iy2,ix2]

    if AngleAvg == True:
        return np.average(bb,1), r
    else:
        return bb


#==============================================================================

def KloAverageGet():

    import numpy as np
    
    global klo_avg_aak_orig
    global klo_avg_aak_polar
    global klo_avg_count
    
    akk_orig  = klo_avg_aak_orig/klo_avg_count
    akk_polar = klo_avg_aak_polar/klo_avg_count
    
    return np.fft.fftshift(akk_orig), akk_polar, klo_avg_count



#==============================================================================

def KloAverageClear():

    global klo_avg_aak_orig
    global klo_avg_aak_polar
    global klo_avg_count

    if "klo_avg_count" in globals():

        del klo_avg_aak_orig
        del klo_avg_aak_polar
        del klo_avg_count

    
#==============================================================================


def KloEvolution(klo, t0, dt, every, fname_out):

    import numpy as np
    
    nfields = 7   # (t, nk, nq, pkx, pky, qkx, qky)
    
    n, Nx, Ny = klo.shape

    if Nx == Ny:
        N = Nx
    else:
        print("Non-square domain")
        return(-1)

    data = np.zeros((int(n/every), nfields))
    
    k = np.fft.ifftshift(np.arange(N)-N/2)

    kx, ky = np.meshgrid(k,k)
    kk = kx*kx + ky*ky
    qx = kx * kk * np.sqrt(kk)
    qy = ky * kk * np.sqrt(kk)
    
    
    for i in np.arange(0, data.shape[0]):

        t  = t0 + dt*i*every
        ak = klo[i*every, :, :]
        aak = np.real(ak * np.conjugate(ak))

        nk  = np.sum(aak)
        nq  = np.sum(aak * kk * kk)
        
        pkx = np.sum(aak * kx)
        pky = np.sum(aak * ky)
        qkx = np.sum(aak * qx)
        qky = np.sum(aak * qy)
        
        data[i,:] = (t, nk, nq, pkx, pky, qkx, qky)
        
    fid=open(fname_out,'ab')
    np.savetxt(fid,data)
    fid.close()

    return


#==============================================================================

def AvgSpectrumPsi(fbase, seeds, istart = 1,  iend = -1):
#
#  reads a sequence of *.psi files and returns 2D array with average spectrum
#  file name in the form "fbase_s00.psi.0000" is assumed
#
#  fbase   string file base
#  seeds   iterable sequence of integer seeds to process
#  istart  first file for each seed
#  iend    first file not to process; last existing if iend=-1 
#
    import numpy as np

    nreads = 0

    for s in seeds:

        fbase_seed =  fbase + '_s' + str(s).zfill(2)  
        i = istart
        print(fbase_seed)
        msg = "For seed {} no files read".format(s)       
        
        while True:

            fname = fbase_seed + '.psi.' + str(i).zfill(4)
            print(fname)
        
            if (os.path.isfile(fname) &  ((iend < 0) | (i < iend)) ) :

                spc0 = ReadPsi(fname, spectrum = True)
                
                if (len(psi) > 0):
                    if "spc" in locals():
                        spc += spc0
                    else:
                        spc = spc0
                    nreads += 1
                    msg = "For seed {} last file read: {}".format(s, fname)
                    
                i += 1
                
            else:
                print(msg)
                break

    print("Number of reads:  " +  str(nreads) )

    if nreads > 0:
        spc = spc/nreads
        return(spc)
    else:
        return(-1)
    
#==============================================================================

def AvgSpectrumKlo(fbase, N, No, seeds, istart = 1,  iend = -1):
#
#  reads a sequence of *.klo files and returns 2D array with average spectrum
#  file name in the form "fbase00/fbase00.klo.0000" is assumed
#
#  fbase   string file base
#  N       gridpoints in each direction in klo file
#  No      gridpoints in each direction in original simulation 
#  seeds   iterable sequence of integer seeds to process
#  istart  first file for each seed
#  iend    first file not to process; last existing if iend=-1 
#   
    import numpy as np
    import os
    
    nreads = 0

    for s in seeds:

        fbase_seed =  fbase + str(s).zfill(2)  
        i = istart
        
        while True:

            fname = fbase_seed + '/' + fbase_seed + '.klo.' + str(i).zfill(4)
        
            if (os.path.isfile(fname) &  ((iend < 0) | (i < iend)) ) :

                spc0, n = ReadKlo(fname, N, No, spectrum = True)
                
                if (len(spc0) > 0):
                    if "spc" in locals():
                        spc += spc0
                    else:
                        spc = spc0
                    nreads += 1
                    
                i += 1
                
            else:
                print("For seed {} last file read: {}".format(s, fname) )
                break

    print("Number of reads:  " +  str(nreads) )

    if nreads > 0:
        spc = spc/nreads
        return(spc)
    else:
        return(-1)
#==============================================================================

def CorrelatorsKlo(fbaseout, fbase, modes, N, No, seeds, istart = 1,  iend = -1, absval=False):
#
#  reads a sequence of *.klo files and returns 2D array with average spectrum
#  file name in the form "fbase00/fbase00.klo.0000" is assumed
#
#  fbase   string file base
#  modes   tuple of modes
#  N       gridpoints in each direction in klo file
#  No      gridpoints in each direction in original simulation 
#  seeds   iterable sequence of integer seeds to process
#  istart  first file for each seed
#  iend    first file not to process; last existing if iend=-1 
#   
    import numpy as np
    import os

    corr = np.zeros((N, N), dtype="complex")
    nk   = np.zeros((N, N))
    nc   = np.zeros(len(modes))

    nreads = 0

    for s in seeds:

        fbase_seed =  fbase + str(s).zfill(2)  
        i = istart
        
        while True:

            fname = fbase_seed + '/' + fbase_seed + '.klo.' + str(i).zfill(4)
        
            if (os.path.isfile(fname) &  ((iend < 0) | (i < iend)) ) :

                fk = ReadKlo(fname, N, No, spectrum = False)
                n = fk.shape[0]
                                
                for nn in range(n):
                    c = 1+0j
                    for mm in range(len(modes)):
                        m = modes[mm]
                        fc = fk[nn, m[1], m[0]]
                        c  = c * np.conj(fc)
                        nc[mm] += np.real( fc*np.conj(fc) )
                        
                    if absval:
                        corr += np.abs(fk[nn,:,:] * c)
                    else:
                        corr += fk[nn,:,:] * c
                        
                    nk   += np.real(fk[nn,:,:] * np.conj(fk[nn,:,:]))
 
                nreads += n
                i += 1
                
            else:
                print("For seed {} last file read: {}".format(s, fname) )
                break

    print("Number of reads:  " +  str(nreads) )

    if nreads > 0:
        corr = corr/nreads
        nk = nk/nreads
        nc = nc/nreads
        for mm in range(len(modes)):
            nk = nk*nc[mm]
        corr = corr/np.sqrt(nk)
        np.save(fbaseout + ".npy", corr)
        return(corr)        
    else:
        return(-1)

#==============================================================================

#==============================================================================

def ReadPsi(fname, spectrum = False):
# if "spectrum = true" return real 2D array 
# else returns complex 2D array

    import numpy as np
  
    dat = np.fromfile(fname, 'float64')
    
    N = int(np.sqrt(len(dat)/2))
    if len(dat)%2%N%N > 0:
        print("Wrong file size: ", fname, len(dat))
        return(0)

    dat = dat.reshape(N,N,2)
    psi = dat[:,:,0] + 1j*dat[:,:,1]  # transpose?

    if spectrum:
        NN = psi.shape[0] * psi.shape[1]                  
        fk = np.fft.fft2(psi)/NN
        spc = np.real(fk * np.conjugate(fk))
        return(spc)
    else:
        return(psi)
 
#===============================================================================


def ReadKlo(fname, N, No, spectrum = False):
# N  - gridpoints in each direction in "klo" file  (N = klow)
# No - gridpoints in each direction in original simulation 
# if "spectrum = true" return real 2D array of averaged data
# else returns complex 3D array with all snapshots
     
    import numpy as np
    
    dat = np.fromfile(fname, 'float64')

    n = int(len(dat)/2/N/N)    
    if len(dat)%2%n%N%N > 0:
        print("Wrong file size: ", fname)
        return(0)
           
    dat = dat.reshape(n,N,N,2) * No
    fk  = dat[:,:,:,0] + 1j*dat[:,:,:,1]  # transpose?

    if spectrum:
        spc = np.average( np.real( fk*np.conjugate(fk) ), 0)
        return spc, n
    else:
        return(fk)
 
#===============================================================================
    
def datasets(fbase):
# not used

    import numpy as np
    import os
    
    sets = list()
    
    for file in os.listdir("."):
        if file.startswith(fbase):
            sets.append(file)
            #print(os.path.join("/mydir", file))

    return(sets)
    
#===============================================================================

def CleanTime(t, debug = False):
#
# input:   t - 1D array of time from *.dat file
# output:  1D array of indexes for consecutive time
#
    import numpy as np

    #-- find left edges of intervals  --
    
    i1 = np.asarray(np.nonzero(t[:-1] > t[1:])).squeeze() + 1

    if debug:
        print("i1 = ", i1, ",   shape = ", i1.shape)

    if i1.shape == (0,):
        return(np.arange(len(t)))

        
    #-- find rigth edges of intervals --
        
    i2 = []

    if i1.shape == ():
        try:
            q = np.max(np.nonzero(t<t[i1]))
        except:
            q = 0
        i2.append(q)

    else:
        for i in i1:
            try:
                q = np.max(np.nonzero(t<t[i]))
            except:
                q = 0
            i2.append(q)

    if debug:
        print("i2 = ", i2)

    #-- find intervals between jumps --
        
    if i2[0] == 0:
        i2 = i2[1:]
        i2.append(len(t)-1)
    else:
        i2.append(len(t)-1)
        i2 = np.asarray(i2) 
        i1 = np.hstack((np.asarray(0),i1))
    
    ind = np.asarray([],dtype=int)
    for i in range(len(i1)):
        ind = np.hstack((ind, np.arange(i1[i], i2[i]+1)))

    if debug:
        print("i1 = ", i1)
        print("i2 = ", i2)
        print(t[ind])

    return(ind)


#===============================================================================

def data_index():
    
    v = {
        "time"    :   0,
        "E_pot"   :   1,  
        "E_kin"   :   2,
        "E_nl"    :   3,
        "eta_rms" :   4,
        "eta_min" :   5,
        "eta_max" :   6,
        "sumA_sq" :   7,
        "sumAk_sq":   8,
        "Px"      :   9,
        "Py"      :  10
    }

    return(v)


#===============================================================================






