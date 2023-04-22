#import importlib ; import nlwTools  ; import nlwTests 
#importlib.reload(nlwTools) ; importlib.reload(nlwTests); nlwTests.Test()

#==============================================================================

def Test():

    #test_klo_vs_psi("test/a1",     80, fnum=1)
    #test_klo_vs_psi("test/tmp_m05", 6, fnum=2)
    
    #test_klo_readavg("test/a1",     80, 256,  fnum=1)
    #test_klo_readavg("test/a1s55", 160, 512,  fnum=124)
    #test_klo_readavg("test/a1s55", 160, 512,  fnum=125)    

    #test_interp_klo("test/a1s55", 160, 512,  fnum=124)

    #test_AvgSpectrumKlo("c1s", 60, 128, [14,15], istart = 20,  iend = -1)
    test_CorrelatorsKlo("c1s+03+05", "c1s", ((3,5),), 60, 128, [14,15], istart = 20,  iend = -1)


#==============================================================================


def test_CorrelatorsKlo(fbaseout, fbase, modes, N, No, seeds, istart = 1,  iend = -1):
#
#  reads a sequence of *.klo files and returns 2D array with average spectrum
#  file name in the form "fbase00/fbase00.klo.0000" is assumed
#
#  fbase   string file base
#  modes   tuple of modes
#  N       gridpoints in each direction
#  No      gridpoints in each direction in original simulation 
#  seeds   iterable sequence of integer seeds to process
#  istart  first file for each seed
#  iend    first file not to process; last existing if iend=-1 
#   

    import numpy as np
    import matplotlib.pyplot as plt
    import nlwTools as nlw

    corr = nlw.CorrelatorsKlo(fbaseout, fbase, modes, N, No, seeds, istart = istart,  iend = iend)

    #-- image --
    
    fig, ax = plt.subplots(ncols=1, nrows=1, constrained_layout=True, figsize=(4,4))

    #a=np.fft.fftshift(np.abs(corr));
    a=np.abs(corr);
    ax.imshow(a, cmap='jet', origin = 'lower', vmin=0, vmax=1);
    print("corr min,max: ",np.min(a), np.max(a))
    
    plt.show();
    
    return()


    
#==============================================================================


def test_AvgSpectrumKlo(fbase, N, No, seeds, istart = 1,  iend = -1):
#
#  reads a sequence of *.klo files and returns 2D array with average spectrum
#  file name in the form "fbase00/fbase00.klo.0000" is assumed
#
#  fbase   string file base
#  N       gridpoints in each direction
#  No      gridpoints in each direction in original simulation 
#  seeds   iterable sequence of integer seeds to process
#  istart  first file for each seed
#  iend    first file not to process; last existing if iend=-1 
#   

    import numpy as np
    import matplotlib.pyplot as plt
    import nlwTools as nlw

    spc = nlw.AvgSpectrumKlo(fbase, N, No, seeds, istart = istart,  iend = iend)

    #return()
    #-- image --
    
    fig, ax = plt.subplots(ncols=1, nrows=1, constrained_layout=True, figsize=(4,4))

    a=np.log10(np.fft.fftshift(spc));
    ax.imshow(a, cmap='jet', origin = 'lower', vmin=-6, vmax=-3);
    print("log10(klo) min,max: ",np.min(a), np.max(a))
    
    plt.show();
    
    return()


    
#==============================================================================

def test_interp_klo(fname, N, N0, fnum = 0):

    import numpy as np
    import matplotlib.pyplot as plt
    import nlwTools as nlw
    
    fname_klo = fname + ".klo." + str(fnum).zfill(4)
        
    klo = nlw.ReadKlo(fname_klo, N, N0, spectrum=False);
    n, N1, N2 = klo.shape

    nlw.KloAverageClear()
    
    nlw.KloAverageAdd(klo, dNk=4, Nphi=32)

    a, b, count = nlw.KloAverageGet()

    print("Spectrum from klo:  {} x {},  n = {} slices, count = {}".format(N1, N2, n, count))

    
    #-- image --
    
    fig, ax = plt.subplots(ncols=3, nrows=1, constrained_layout=True, figsize=(12,4))
    
    ax[0].imshow(a, cmap='jet', origin = 'lower', vmin=0, vmax=0.002);
    print("min, max: ", np.min(a), np.max(a))

    ax[1].imshow(b, cmap='jet', origin = 'lower');
    print("min, max: ", np.min(b), np.max(b))

    ax[2].plot(np.arange(b.shape[1]), np.sum(b,0), '-o');

    
    plt.show();
    
    return()

#==============================================================================


def test_klo_readavg(fname, N, N0, fnum = 0):

    import numpy as np
    import matplotlib.pyplot as plt
    import nlwTools as nlw
    
    n = 1  # number of snapshots in "klo" file, use n=1 for better comaprison 
    L = 2*np.pi

    fname_klo = fname + ".klo." + str(fnum).zfill(4)
        
    klo, n = nlw.ReadKlo(fname_klo, N, N0, spectrum=True);
    (N1, N2) = klo.shape
    print("Spectrum from klo:  {} x {}, sum = {},  n = {}".format(N1, N2, np.sum(klo), n))

    klo = np.fft.fftshift(klo)

    
    #-- image --
    
    fig, ax = plt.subplots(ncols=1, nrows=1, constrained_layout=True, figsize=(4,4))

    a=np.log10(klo);
    ax.imshow(a, cmap='jet', origin = 'lower', vmin=-6, vmax=-3);
    print("log10(klo) min,max: ",np.min(a), np.max(a))
    
    plt.show();
    
    return()

#==============================================================================

def test_klo_vs_psi(fname, Nklo, fnum = 0):
#
# On normalization of FFT.
#
# See "*.dat" file for:  8.sumA_sq  9.sumAk_sq
#
#    #8.sumA_sq   = (1/2) * sum(|psi|^2) dxdy
#    #9.sumAk_sq  = sum( |A_k|^2 ) 
#    (#9) = (2/L^2)*(#8)
#
#    From "test/tmp_m05.dat": (#8) = 1.577169407370e-01,  (#9) = 7.990033557959e-03
#    From "test/a1.dat":      (#8) = 2.825774223384e-01,  (#9) = 1.431553945097e-02
    
    import numpy as np
    import matplotlib.pyplot as plt
    import nlwTools as nlw

    
    L = 2*np.pi

    fname_psi = fname + ".psi." + str(fnum).zfill(4)
    fname_klo = fname + ".klo." + str(fnum).zfill(4)
     
    #-- compare integral quantities --

    psi=0
    
    psi = nlw.ReadPsi(fname_psi,spectrum=False);
    (N1, N2) = psi.shape
    sumA_sq = 0.5 * np.sum( np.real(psi * np.conjugate(psi)) ) * L/N1 * L/N2
    sumAk_sq = sumA_sq * 2/L/L 
    
    print("Integral of psi:  {} x {}, integral = {} ({})".format(
        N1, N2, sumA_sq, sumAk_sq))
    
    psi = nlw.ReadPsi(fname_psi, spectrum=True);
    (N1, N2) = psi.shape
    print("Spectrum from psi:  {} x {}, sum = {}".format(N1, N2, np.sum(psi)))

    N0 = N1
    
    klo, n = nlw.ReadKlo(fname_klo, Nklo, N0, spectrum=True);
    (N1, N2) = klo.shape
    print("Spectrum from klo:  {} x {}, sum = {}".format(N1, N2, np.sum(klo)))

          
    #-- recenter and trim "psi" to match "klo"

    klo = np.fft.fftshift(klo)
    psi = np.fft.fftshift(psi)
          
    i1 = int(N0/2 - N1/2)
    i2 = i1 + N1
          
    psi = psi[i1:i2, i1:i2]
    print("Trimmed psi:  {} x {}, sum = {}".format(N1, N1, np.sum(psi)))

    #-- images --
    
    fig, ax = plt.subplots(ncols=2, nrows=1, constrained_layout=True, figsize=(6,3))

    a=np.log10(klo);
    ax[0].imshow(a, cmap='jet', origin = 'lower', vmin=-6, vmax=-3);
    print("log10(klo) min,max: ",np.min(a), np.max(a))

    a=np.log10(psi);
    ax[1].imshow(a, cmap='jet', origin = 'lower', vmin=-6, vmax=-3);
    print("log10(psi) min,max: ",np.min(a), np.max(a))
    
    plt.show();
    
    return()


#==============================================================================



#===============================================================================






