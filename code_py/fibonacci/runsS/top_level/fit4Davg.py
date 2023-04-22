import numpy as np
import fiboMI as mi
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

#-------------

#a=-1
#fbase = "MIx32x1/s4m10_i80" ; title="$\\alpha=1/2$ (direct), $|i-d|$=20" ; outfile="fit4Do_s4m10_i80.pdf"
#fbase = "MIx32x1/s4m10_i60" ; title="$\\alpha=1/2$ (direct), $|i-d|$=40" ; outfile="fit4Do_s4m10_i60.pdf" 
#fbase = "MIx32x1/s4m10_i40" ; title="$\\alpha=1/2$ (direct), $|i-d|$=60" ; outfile="fit4Do_s4m10_i40.pdf"
#fbase = "MIx32x1/s4m10_i20" ; title="$\\alpha=1/2$ (direct), $|i-d|$=80" ; outfile="fit4Do_s4m10_i20.pdf"

#a=1
#fbase = "MIx32x1/s4m90_i20" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=20" ; outfile="fit4Do_s4m90_i20.pdf"
#fbase = "MIx32x1/s4m90_i40" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=40" ; outfile="fit4Do_s4m90_i40.pdf"
#fbase = "MIx32x1/s4m90_i60" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=60" ; outfile="fit4Do_s4m90_i60.pdf"
#fbase = "MIx32x1/s4m90_i80" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=80" ; outfile="fit4Do_s4m90_i80.pdf"


#-------------

#a=-1
#fbase = "MIx32x1/s5m10_i180" ; title="$\\alpha=1/2$ (direct), $|i-d|$=20" ; outfile='fit4Do_s5m10_i180.pdf'
#fbase = "MIx32x1/s5m10_i160" ; title="$\\alpha=1/2$ (direct), $|i-d|$=40" ; outfile='fit4Do_s5m10_i160.pdf'
#fbase = "MIx32x1/s5m10_i140" ; title="$\\alpha=1/2$ (direct), $|i-d|$=60" ; outfile='fit4Do_s5m10_i140.pdf'
#fbase = "MIx32x1/s5m10_i120" ; title="$\\alpha=1/2$ (direct), $|i-d|$=80" ; outfile='fit4Do_s5m10_i120.pdf'
#fbase = "MIx32x1/s5m10_i80" ; title="$\\alpha=1/2$ (direct), $|i-d|$=120" ; outfile='fit4Do_s5m10_i080.pdf'
#fbase = "MIx32x1/s5m10_i40" ; title="$\\alpha=1/2$ (direct), $|i-d|$=160" ; outfile='fit4Do_s5m10_i040.pdf'
#fbase = "MIx32x1/s5m10_i20" ; title="$\\alpha=1/2$ (direct), $|i-d|$=180" ; outfile='fit4Do_s5m10_i020.pdf'

a=1
#fbase = "MIx32x1/s5m90_i20" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=20" ; outfile='fit4Do_s5m90_i020.pdf' 
#fbase = "MIx32x1/s5m90_i40" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=40" ; outfile='fit4Do_s5m90_i040.pdf' 
#fbase = "MIx32x1/s5m90_i60" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=60" ; outfile='fit4Do_s5m90_i060.pdf' 
#fbase = "MIx32x1/s5m90_i80" ; title="$\\alpha=1/2$ (inverse), $|i-d|$=80" ; outfile='fit4Do_s5m90_i080.pdf' 
#fbase = "MIx32x1/s5m90_i120"; title="$\\alpha=1/2$ (inverse), $|i-d|$=120"; outfile='fit4Do_s5m90_i120.pdf' 
#fbase = "MIx32x1/s5m90_i160"; title="$\\alpha=1/2$ (inverse), $|i-d|$=160"; outfile='fit4Do_s5m90_i160.pdf' 
#fbase = "MIx32x1/s5m90_i180"; title="$\\alpha=1/2$ (inverse), $|i-d|$=180"; outfile='fit4Do_s5m90_i180.pdf' 



#-------------


#outfile='fit4Davg.pdf'


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
c = np.sqrt(nn1*nn2*nn3/nnn)
dPc = dP/c
dPc[P123 == 0]=0



Ps = dPc; Cs = P123; 
avg = np.copy(np.sum(Ps*Cs, axis=(0,1,2))/np.sum(Cs, axis=(0,1,2) ))
h=np.average(avg)

q0 = np.copy(dPc[0,0,0,:]);
h1=np.average(q0)


Ps = dPc[:2, :2, :2, :];  Cs = P123[:2, :2, :2, :]; 
q1 = np.copy(np.sum(Ps*Cs, axis=(0,1,2))/np.sum(Cs, axis=(0,1,2) ))

Ps = dPc[:3, :3, :3, :];  Cs = P123[:3, :3, :3, :]; 
q2 = np.copy(np.sum(Ps*Cs, axis=(0,1,2))/np.sum(Cs, axis=(0,1,2) ))

Ps = dPc[:4, :4, :4, :];  Cs = P123[:4, :4, :4, :]; 
q3 = np.copy(np.sum(Ps*Cs, axis=(0,1,2))/np.sum(Cs, axis=(0,1,2) ))




Ps = dPc[1:, 1:, 1:, :];  Cs = P123[1:, 1:, 1:, :]; 
p1 = np.copy(np.sum(Ps*Cs, axis=(0,1,2))/np.sum(Cs, axis=(0,1,2) ))

Ps = dPc[2:, 2:, 2:, :];  Cs = P123[2:, 2:, 2:, :];
p2 = np.copy(np.sum(Ps*Cs, axis=(0,1,2))/np.sum(Cs, axis=(0,1,2) ))

Ps = dPc[3:, 3:, 3:, :];  Cs = P123[3:, 3:, 3:, :];
p3 = np.copy(np.sum(Ps*Cs, axis=(0,1,2))/np.sum(Cs, axis=(0,1,2) ))



Ps = dPc;  Ps[0,0,0, :]=0;  Cs = P123;  Cs[0,0,0,:]=0;
p0 = np.sum(Ps*Cs, axis=(0,1,2))/np.sum(Cs, axis=(0,1,2) )
h2=np.average(p0)


#print(lnP123.shape, eqP.shape)


#-- plot data --
MS = 1
LW = 1


fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(3.5,6.5))

iax=ax[0]

iax.plot(phi, q0, '-r', label="$x<1$")
iax.plot(phi, q1, '-b', label="$x<2$")
iax.plot(phi, q2, '-g', label="$x<3$")
iax.plot(phi, q3, '-y', label="$x<4$")
iax.plot(phi, a*np.sin(phi*pi)+h, '-k', lw=0.5)
iax.plot(phi, avg, '--k')
iax.plot(phi, a*np.sin(phi*pi)+h1, '-k', lw=0.5)


iax=ax[1]

iax.plot(phi, p0, '-r', label="no $x<1$")
iax.plot(phi, p1, '-b', label="$x>1$")
iax.plot(phi, p2, '-g', label="$x>2$")
iax.plot(phi, p3, '-y', label="$x>3$")
iax.plot(phi, avg, '--k')
iax.plot(phi, a*np.sin(phi*pi)+h2, '-k', lw=0.5)



#-----------------------------------------------

ax[0].set_title(title)
ax[1].set(xlabel = "$\\theta/\pi$")
ax[0].set_ylim(-1,5)
ax[1].set_ylim(-2,4)


for i in (0,1):
    iax=ax[i]
    iax.set_xlim(-1,1)
    iax.set(ylabel='$[\, \ln {\cal P} - \ln {\cal P}_{eq}\,] \, / \, \sqrt{n_1 n_2 n_3 / x_1 x_2 x_3}$')
    iax.legend(frameon=False)


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




