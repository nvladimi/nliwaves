
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
outfile = "plotJpdf.pdf"
binsize = 1
numbins = 20000
nmodes  = 200


fig, ax = plt.subplots(ncols=5, nrows=1, figsize=(15,3))
fbase = 'Post01/s5m10binA'

Pre  = np.fromfile(fbase+"_JPre.dat", 'int32').reshape([numbins, nmodes])
Pim  = np.fromfile(fbase+"_JPim.dat", 'int32').reshape([numbins, nmodes])

x = (np.arange(-numbins/2,numbins/2) + 0.5) * binsize;

for i in (50, 100, 150):

    y =  Pre[:,i]/np.sum(Pre[:,i])/binsize
    ax[0].plot(x,y, '-',  mfc='none', ms=1,  lw=0.5,  label="i = {}".format(i) )
    ax[1].plot(x,y, '-',  mfc='none', ms=1,  lw=0.5)
    ax[2].plot(np.abs(x),y, '-',  mfc='none', ms=1,  lw=0.5)
    print("Real, i = ", i," :  ", np.sum(y*x)*binsize, np.sum(y*x*x)*binsize)
 
    
    y =  Pim[:,i]/np.sum(Pim[:,i])/binsize
    #ax[0][0].plot(x,y, '-',  mfc='none', ms=1,  lw=0.5)
    #ax[1][1].plot(np.abs(x),y, '-',  mfc='none', ms=1,  lw=0.5)
    print("Imag, i = ", i," :  ", np.sum(y*x)*binsize, np.sum(y*x*x)*binsize)
 
    
re1=np.zeros(nmodes); re2=np.zeros(nmodes); im1=np.zeros(nmodes); im2=np.zeros(nmodes);

for i in range(nmodes):
    y =  Pre[:,i]/np.sum(Pre[:,i])/binsize
    re1[i] = np.sum(y*x)*binsize
    re2[i] = np.sum(y*x*x)*binsize
    y =  Pim[:,i]/np.sum(Pim[:,i])/binsize
    im1[i] = np.sum(y*x)*binsize
    im2[i] = np.sum(y*x*x)*binsize
    
i=np.arange(200);
ax[3].plot(i,re2/(200-i)**2, '-',  mfc='none', ms=1,  lw=0.5)
ax[3].plot(i,im2/(200-i)**2, '--',  mfc='none', ms=1,  lw=0.5)

ax[4].plot(i,re1, '-',  mfc='none', ms=1,  lw=0.5)
ax[4].plot(i,im1, '--',  mfc='none', ms=1,  lw=0.5)

    
#-- canvas options --




iax = ax[0]  
iax.set_xlim(-100,100)
#iax.set_ylim(0,8)
#iax.set(xlabel='$i$', ylabel='$n=C_0$')
#iax.grid(axis="y")
#iax.set_yscale('log')



iax = ax[1]  
iax.set_xlim(-2000,2000)
iax.set_ylim(1e-8,0.1)
iax.set_yscale('log')

iax = ax[2] 
iax.set_xlim(1,2000)
iax.set_ylim(1e-8,0.1)
iax.set_xscale('log')
iax.set_yscale('log')

iax = ax[3] 
iax.set_ylim(0,1)

iax = ax[4] 
iax.set_ylim(-0.8,0.8)


ax[0].legend(frameon=False, loc="upper left") # bbox_to_anchor=(1.02, 0.5))


#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.06, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
