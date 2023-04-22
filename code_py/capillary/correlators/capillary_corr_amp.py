
import matplotlib
#matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import numpy as np
from matplotlib.colors import LinearSegmentedColormap




T = [[  4,  -2],  [  4,   8]]

#Q = [[ 5,  5],  [5,  -5], [0, 5]]
Q = [[ 5,  5],  [5,  -5], [0, 8]]



fnameout = "capillary_corr_amp.pdf"
prefix = "../pyPost/"
N=60;

runs = ("c1", "c5")
nicknames = ("Low forcing", "High forcing")

cmax = (0.01, 0.005)


#------------------------------------------

ucolors = [(0.95, 0.95, 0.95), (0.9, 0.8, 0), (0.7, 0, 0)]  # R -> G -> B
ucmap = LinearSegmentedColormap.from_list("ucolors", ucolors, N=8)

mmax=12
i1 = int(N/2-mmax)
i2 = int(N/2+mmax+1)
extent=(-mmax-0.5, mmax+0.5,-mmax-0.5, mmax+0.5)

ms=7


#-------------------------------

#fig, ax = plt.subplots(ncols=2, nrows=2, constrained_layout=True, figsize=(8,8))
fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6,6))

for i in (0,1):

    #-- triple corrs: -------------------------
    iax = ax[0][i]
    run = runs[i]

    corr = "{:+03d}{:+03d}{:+03d}{:+03d}".format(T[0][0], T[0][1], T[1][0], T[1][1])
    Tsum = np.sum(T,0)

    #-- data --    
    a = np.load(prefix + run + corr + ".npy")
    c = np.fft.fftshift(a);
    c = c[i1:i2, i1:i2]
    Csum = np.abs(a[Tsum[1],Tsum[0]])

    im0=iax.imshow(np.abs(c),  cmap=ucmap, origin = 'lower', vmin=0, vmax=cmax[0], extent=extent);

    #-- annotations --

    iax.set_title(nicknames[i] + ": $C_{res} = $" + "{:.3f}".format(Csum) )
 
    kmodes=( (T[0][0], T[1][0], Tsum[0]), (T[0][1], T[1][1], Tsum[1]))

    iax.plot(kmodes[0], kmodes[1], 'sk', ms=ms, linewidth=0.5, mfc='none')
    iax.plot(0, 0, '+k', ms=ms, linewidth=0.5)
    iax.plot(Tsum[0],Tsum[1], 'xw', ms=ms, linewidth=0.5)


    #-- quad corrs: -------------------------
    
    iax = ax[1][i]

    corr = "{:+03d}{:+03d}{:+03d}{:+03d}{:+03d}{:+03d}".format(
        Q[0][0], Q[0][1], Q[1][0], Q[1][1], Q[2][0], Q[2][1])
    Qsum = np.sum(Q,0)

    #-- data --    
    a = np.load(prefix + run + corr + ".npy")
    c = np.fft.fftshift(a);
    c = c[i1:i2, i1:i2]
    Csum = np.abs(a[Qsum[1],Qsum[0]])

    im1=iax.imshow(np.abs(c),  cmap=ucmap, origin = 'lower', vmin=0, vmax=cmax[1], extent=extent);

    #-- annotations --

    iax.set_title(nicknames[i] + ": $C_{res} = $" + "{:.3f}".format(Csum) )

    kmodes=( (Q[0][0], Q[1][0], Q[2][0], Qsum[0]), (Q[0][1], Q[1][1], Q[2][1], Qsum[1]))

    iax.plot(kmodes[0], kmodes[1], 'sk', ms=ms, linewidth=0.5, mfc='none')
    iax.plot(0, 0, '+k', ms=ms, linewidth=0.5)
    iax.plot(Qsum[0],Qsum[1], 'xk', ms=ms, linewidth=0.5)


#------------------------------

Ttitle = "Triple correlators:   [{}, {}] + [{}, {}] = [{}, {}]".format(
    T[0][0], T[0][1], T[1][0], T[1][1], Tsum[0], Tsum[1])

Qtitle = "Quadruple correlators:   [{}, {}] + [{}, {}] + [{}, {}] = [{}, {}]".format(
     Q[0][0], Q[0][1], Q[1][0], Q[1][1],  Q[2][0], Q[2][1], Qsum[0], Qsum[1])


ax[0][0].text(-5,16, Ttitle, fontsize=12)
ax[1][0].text(-11.5,16, Qtitle, fontsize=12)


cbar0=fig.colorbar(im0, ax=ax[0,1], orientation='vertical', shrink=0.8)
#cbar0.set_ticks([0,cmax/2, cmax])
#cbar0.set_ticklabels([0,cmax/2, cmax])

cbar1=fig.colorbar(im1, ax=ax[1,1], orientation='vertical', shrink=0.8)
#cbar0.set_ticks([0,cmax/2, cmax])
#cbar0.set_ticklabels([0,cmax/2, cmax])



#plt.show()
#plt.close()


plt.subplots_adjust(left=0.02, right=0.96, top=0.92, bottom=0.04, wspace=0.06, hspace=0.37)
plt.savefig(fnameout, pad_inches=0)



exit()


