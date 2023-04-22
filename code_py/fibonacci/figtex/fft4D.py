import numpy as np
import fiboMI as mi
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import matplotlib.pylab as plt

#-------------


#fbase = "MIx32x1/s5m10_i180" ; title="$\\alpha=1/2$ (direct), $|i-d|$=20" ; outfile='fft4D_s5m10_i180.pdf'
#fbase = "MIx32x1/s5m10_i160" ; title="$\\alpha=1/2$ (direct), $|i-d|$=40" ; outfile='fft4D_s5m10_i160.pdf'
#fbase = "MIx32x1/s5m10_i140" ; title="$\\alpha=1/2$ (direct), $|i-d|$=60" ; outfile='fft4D_s5m10_i140.pdf'
#fbase = "MIx32x1/s5m10_i120" ; title="$\\alpha=1/2$ (direct), $|i-d|$=80" ; outfile='fft4D_s5m10_i120.pdf'
#fbase = "MIx32x1/s5m10_i80" ; title="$\\alpha=1/2$ (direct), $|i-d|$=120" ; outfile='fft4D_s5m10_i080.pdf'
#fbase = "MIx32x1/s5m10_i40" ; title="$\\alpha=1/2$ (direct), $|i-d|$=160" ; outfile='fft4D_s5m10_i040.pdf'
#fbase = "MIx32x1/s5m10_i20" ; title="$\\alpha=1/2$ (direct), $|i-d|$=180" ; outfile='fft4D_s5m10_i020.pdf'



#fbase = "MIx32x1/s5m90_i20" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=20" ; outfile='fft4D_s5m90_i020.pdf' 
#fbase = "MIx32x1/s5m90_i40" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=40" ; outfile='fft4D_s5m90_i040.pdf' 
#fbase = "MIx32x1/s5m90_i60" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=60" ; outfile='fft4D_s5m90_i060.pdf' 
#fbase = "MIx32x1/s5m90_i80" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=80" ; outfile='fft4D_s5m90_i080.pdf' 
#fbase = "MIx32x1/s5m90_i120" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=120" ; outfile='fft4D_s5m90_i120.pdf' 
#fbase = "MIx32x1/s5m90_i160" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=160" ; outfile='fft4D_s5m90_i160.pdf' 
fbase = "MIx32x1/s5m90_i180" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=180" ; outfile='fft4D_s5m90_i180.pdf' 


#-------------

pi = np.pi

imode, numbins, numbinsPhi, binsize, (navg1, navg2, navg3), ncount =  mi.ReadParamMI(fbase)

dn1  = binsize #*navg1
dn2  = binsize #*navg2
dn3  = binsize #*navg3
dphi = 2*pi/numbinsPhi
nnn = navg1*navg2*navg3

P123 = np.fromfile(fbase + "_P123.dat", 'int32')

#-- normalize and reshape probabilities --

P123 = P123 / (ncount *dn1*dn2*dn3*dphi)
P123 = P123.reshape(numbins, numbins, numbins, numbinsPhi)
#print(np.sum(P123)*dn1*dn2*dn3*dphi)

lnP123 = np.log(P123)

n1 = (np.arange(numbins) + 0.5)*binsize;  #*navg1
n2 = (np.arange(numbins) + 0.5)*binsize;  #*navg2
n3 = (np.arange(numbins) + 0.5)*binsize;  #*navg2

phi =  (np.arange(-numbinsPhi/2, numbinsPhi/2) + 0.5 )/numbinsPhi*2

nn1, nn2, nn3, pp = np.meshgrid(n1, n2, n3, phi, sparse=False)


dPeff = np.log( 2/binsize * np.sinh(binsize/2) )
eqP = - (nn1 + nn2 + nn3) + 3*dPeff - np.log(2*pi)
dP = lnP123 - eqP

fP = np.fft.ifft(dP, axis=3);
a0 =   np.abs(fP[:,:,:,0])
a1 = 2*np.abs(fP[:,:,:,1])
a2 = 2*np.abs(fP[:,:,:,2])
a3 = 2*np.abs(fP[:,:,:,3])

nn1, nn2, nn3 = np.meshgrid(n1, n2, n3, sparse=False)
n = np.sqrt(nn1*nn2*nn3)

n  = np.reshape(n,  numbins**3)
a0 = np.reshape(a0, numbins**3)
a1 = np.reshape(a1, numbins**3)
a2 = np.reshape(a2, numbins**3)
a3 = np.reshape(a3, numbins**3)


dn=1
nmax=12
imax = np.int(nmax/dn)

a0avg = np.zeros(imax)
a1avg = np.zeros(imax)
a2avg = np.zeros(imax)
a3avg = np.zeros(imax)

for i in np.arange(0,imax):
    ind = np.logical_and(n>i*dn, n<(i+1)*dn)
    ind  = np.logical_and(ind, np.isfinite(a0))
    a0avg[i] = np.average(a0[ind])
    a1avg[i] = np.average(a1[ind])
    a2avg[i] = np.average(a2[ind])
    a3avg[i] = np.average(a3[ind])


#-- plot data --
MS = 1
LW = 1


fig, ax = plt.subplots(ncols=1, nrows=4, figsize=(3.5,14))

ax[0].plot(n, a0, 'o', mfc='none', ms=MS,  lw=LW)
ax[1].plot(n, a1, 'o', mfc='none', ms=MS,  lw=LW)
ax[2].plot(n, a2, 'o', mfc='none', ms=MS,  lw=LW)
ax[3].plot(n, a3, 'o', mfc='none', ms=MS,  lw=LW)

ni = (np.arange(0,imax) + 0.5)*dn


ax[0].plot(ni, a0avg, 'o-r', mfc='none', ms=MS*3,  lw=LW)
ax[1].plot(ni, a1avg, 'o-r', mfc='none', ms=MS*3,  lw=LW)
ax[2].plot(ni, a2avg, 'o-r', mfc='none', ms=MS*3,  lw=LW)
ax[3].plot(ni, a3avg, 'o-r', mfc='none', ms=MS*3,  lw=LW)


for i in (0,1,2,3):
    #ax[i].set_yscale('log')
    #ax[i].set_xscale('log')
    ax[i].set_xlim(0,12)
    ax[i].set(xlabel='$\sqrt{x_1 x_2 x_3}$')

ax[0].set_title(title + ",  " + '$\\ln {\cal P} / {\cal P}_{eq}$')
ax[0].set(ylabel='0th harmonic')
ax[1].set(ylabel='1st harmonic')
ax[2].set(ylabel='2nd harmonic')
ax[3].set(ylabel='3rd harmonic')

#-----------------------------------------------


#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.15, wspace=0.3)
plt.subplots_adjust(left=0.2, right=0.95, top=0.97, bottom=0.05, hspace=0.16)

plt.savefig(outfile) #pad_inches=0

plt.close()


#-------------------------------------------------



#import matplotlib.pyplot as plt
#x = np.arange(-5, 5, 0.1)
#y = np.arange(-5, 5, 0.1)
#xx, yy = np.meshgrid(x, y, sparse=True)
#z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)

#import matplotlib.pyplot as plt
#h = plt.contourf(x, y, z)
#plt.axis('scaled')
#plt.show()




