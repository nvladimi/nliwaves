import numpy as np
import fiboMI as mi
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import matplotlib.pylab as plt

#-------------


#fbase = "MIx32x1/s4m10_i80" ; title="$\\alpha=1/2$ (direct), $|i-d|$=20" ; outfile='fit4D_s4m10_i80.pdf'
#fbase = "MIx32x1/s4m10_i60" ; title="$\\alpha=1/2$ (direct), $|i-d|$=40" ; outfile='fit4D_s4m10_i60.pdf'
#fbase = "MIx32x1/s4m10_i40" ; title="$\\alpha=1/2$ (direct), $|i-d|$=60" ; outfile='fit4D_s4m10_i40.pdf'
#fbase = "MIx32x1/s4m10_i20" ; title="$\\alpha=1/2$ (direct), $|i-d|$=80" ; outfile='fit4D_s4m10_i20.pdf'


#fbase = "MIx32x1/s4m90_i20" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=20" ; outfile='fit4D_s4m90_i20.pdf' 
#fbase = "MIx32x1/s4m90_i40" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=40" ; outfile='fit4D_s4m90_i40.pdf' 
#fbase = "MIx32x1/s4m90_i60" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=60" ; outfile='fit4D_s4m90_i60.pdf' 
#fbase = "MIx32x1/s4m90_i80" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=80" ; outfile='fit4D_s4m90_i80.pdf' 


#-------------


#fbase = "MIx32x1/s5m10_i180" ; title="$\\alpha=1/2$ (direct), $|i-d|$=20" ; outfile='fit4D_s5m10_i180.pdf'
#fbase = "MIx32x1/s5m10_i160" ; title="$\\alpha=1/2$ (direct), $|i-d|$=40" ; outfile='fit4D_s5m10_i160.pdf'
#fbase = "MIx32x1/s5m10_i140" ; title="$\\alpha=1/2$ (direct), $|i-d|$=60" ; outfile='fit4D_s5m10_i140.pdf'
#fbase = "MIx32x1/s5m10_i120" ; title="$\\alpha=1/2$ (direct), $|i-d|$=80" ; outfile='fit4D_s5m10_i120.pdf'
#fbase = "MIx32x1/s5m10_i80" ; title="$\\alpha=1/2$ (direct), $|i-d|$=120" ; outfile='fit4D_s5m10_i080.pdf'
#fbase = "MIx32x1/s5m10_i40" ; title="$\\alpha=1/2$ (direct), $|i-d|$=160" ; outfile='fit4D_s5m10_i040.pdf'
#fbase = "MIx32x1/s5m10_i20" ; title="$\\alpha=1/2$ (direct), $|i-d|$=180" ; outfile='fit4D_s5m10_i020.pdf'



#fbase = "MIx32x1/s5m90_i20" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=20" ; outfile='fit4D_s5m90_i020.pdf' 
#fbase = "MIx32x1/s5m90_i40" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=40" ; outfile='fit4D_s5m90_i040.pdf' 
#fbase = "MIx32x1/s5m90_i60" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=60" ; outfile='fit4D_s5m90_i060.pdf' 
#fbase = "MIx32x1/s5m90_i80" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=80" ; outfile='fit4D_s5m90_i080.pdf' 
#fbase = "MIx32x1/s5m90_i120" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=120" ; outfile='fit4D_s5m90_i120.pdf' 
#fbase = "MIx32x1/s5m90_i160" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=160" ; outfile='fit4D_s5m90_i160.pdf' 
fbase = "MIx32x1/s5m90_i180" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=180" ; outfile='fit4D_s5m90_i180.pdf' 


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

#print(lnP123.shape, eqP.shape)


#-- plot data --
MS = 1
LW = 1

def color(a):
    amin=2  #1.5
    amax=4.3
    ncolors = 100
    c = plt.cm.plasma((a-amin)/(amax-amin))
    #c = plt.cm.gist_rainbow((a-amin)/(amax-amin))
    return(c)


fig, ax = plt.subplots(ncols=1, nrows=3, figsize=(3.5,9.5))

print ("nnn:", navg1, navg2, navg3, nnn)

i123 = (
    ((1, 1, 2), (1, 2, 1), (2, 1, 1)),
    ((1, 2, 2), (2, 1, 2), (2, 2, 1)),
    ((1, 1, 4), (1, 4, 1), (4, 1, 1)),
    ((2,2,2),),
    ((1,2,4), (1,4,2), (2,1,4), (2,4,1), (4,1,2), (4,2,1))
    )

for iset in i123:
    q = np.zeros(numbinsPhi)
    for (i1, i2, i3)  in iset:
        a = np.sqrt(n1[i1]*n2[i2]*n3[i3])
        q +=  dP[i1,i2,i3,:];
    q = q/len(iset)    
    fq = np.fft.fft(q);
    fq[5:-4] = 0
    q= np.fft.ifft(fq);
    print(n1[i1], n2[i2], n3[i3], '   ',  a)
    ax[0].plot(phi, q, color=color(a), mfc='none', ms=MS,  lw=LW)
    ax[1].plot(phi, q/a, color=color(a), mfc='none', ms=MS,  lw=LW)
    ax[2].plot(phi, q*a, color=color(a), mfc='none', ms=MS,  lw=LW)



"""

for i1 in (1,2,4):
    for i2 in (1,2,4):
        for i3 in (1,2,4):
            a = np.sqrt(n1[i1]*n2[i2]*n3[i3])
            b = a/np.sqrt(nnn)
            if (a<5):
                q =   dP[i1,i2,i3,:];
                fq = np.fft.fft(q);
                fq[5:-4] = 0
                q= np.fft.ifft(fq);
                print(n1[i1], n2[i2], n3[i3], '   ',  a)
                ax[0].plot(phi, q, color=color(a), mfc='none', ms=MS,  lw=LW)
                ax[1].plot(phi, q/a, color=color(a), mfc='none', ms=MS,  lw=LW)
                ax[2].plot(phi, q*a, color=color(a), mfc='none', ms=MS,  lw=LW)
"""

#-----------------------------------------------


iax = ax[0]
iax.set_title(title)
iax.set_xlim(-1,1)
#iax.set_ylim(-0.4,0.4)
#iax.set_ylim(-0.01,0.07)
iax.set(ylabel='$\\ln {\cal P} - \ln {\cal P}_{eq}$')


iax = ax[1]  
iax.set_xlim(-1,1)
#iax.set_ylim(-2,2)
#iax.set_ylim(-1.0,5.0)
#iax.set(xlabel='$\\theta/\pi$')
iax.set(ylabel='$[\, \ln {\cal P} - \ln {\cal P}_{eq}\,] \, / \, \sqrt{x_1 x_2 x_3}$')
#iax.set_yticks(np.arange(-0.1, 0.11, step=0.05)) 
#iax.set_yticks(np.arange(-0, 0.2, step=0.05)) 

iax = ax[2]  
iax.set_xlim(-1,1)
#iax.set_ylim(-2,2)
#iax.set_ylim(-1.0,5.0)
iax.set(xlabel='$\\theta/\pi$')
iax.set(ylabel='$[\, \ln {\cal P} - \ln {\cal P}_{eq}\,]  \, sqrt{x_1 x_2 x_3}$')
#iax.set_yticks(np.arange(-0.1, 0.11, step=0.05)) 
#iax.set_yticks(np.arange(-0, 0.2, step=0.05)) 






#iax.set_ylim(0,2.0)
#iax.grid(axis="y")
#iax.legend(frameon=False)
#iax.text(2, 1.50, "direct")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_xscale('log')
#iax.set_yticks(np.arange(0, 2.1, step=0.5)) 


#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.15, wspace=0.3)
plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.07, hspace=0.14)

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




