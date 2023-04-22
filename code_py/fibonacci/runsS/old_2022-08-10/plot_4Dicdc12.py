import numpy as np
import fiboMI as mi
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

#------------------------------

runsInverse= {"MIx32x1/s4m90": (("20", "$|i-d| = 20$", "-r"), 
                                ("40", "$|i-d| = 40$", "-b"), 
                                ("60", "$|i-d| = 60$", "-g"), 
                                ("80", "$|i-d| = 80$", "-y")),
              "MIx32x1/s5m90": (("20", "$|i-d| = 20$", "--r"), 
                                ("40", "$|i-d| = 40$", "--b"), 
                                ("60", "$|i-d| = 50$", "--g"), 
                                ("80", "$|i-d| = 80$", "--y"), 
                                ("120","$|i-d| = 120$", ":m"), 
                                ("160","$|i-d| = 160$", ":c"), 
                                ("180","$|i-d| = 180$", ":y")), 
} 
 


runsDirect = {"MIx32x1/s4m10": (("80", "$|i-d| = 20$", "-r"), 
                                ("60", "$|i-d| = 40$", "-b"), 
                                ("40", "$|i-d| = 60$", "-g"), 
                                ("20", "$|i-d| = 80$", "-y")),
              "MIx32x1/s5m10": (("180", "$|i-d| = 20$", "--r"), 
                                ("160", "$|i-d| = 40$", "--b"), 
                                ("140", "$|i-d| = 50$", "--g"), 
                                ("120", "$|i-d| = 80$", "--y"), 
                                ("80", "$|i-d| = 120$", ":m"), 
                                ("40", "$|i-d| = 160$", ":c"), 
                                ("20", "$|i-d| = 180$", ":y")), 
} 
 



#runs=runsDirect;  outfile='plot_4Ddir12.pdf' ; ylim=(-1.5,1); title="$\\alpha = 1/2$ (direct)"

runs=runsInverse;  outfile='plot_4Dinv12.pdf' ; ylim=(-1,1.5); title="$\\alpha = 1/2$ (inverse)"


MS = 1
LW = 1
pi = np.pi

fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(14,3.5))


#-----------------------------

x=runs.keys()
print(x)


#run="MIx32x1/s4m10"

for run in runs.keys():

    print(run)

    for r in runs[run]:

        print(r)

        fbase = run + "_i" + r[0]
        label = r[1]
        style = r[2]
        print(fbase)

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
        avg = avg - np.average(avg)

        q0 = np.copy(dPc[0,0,0,:]);
        q0 = q0 - np.average(q0)

        Ps = dPc[1:, 1:, 1:, :];  Cs = P123[1:, 1:, 1:, :]; 
        q1 = np.copy(np.sum(Ps*Cs, axis=(0,1,2))/np.sum(Cs, axis=(0,1,2) ))
        q1 = q1 - np.average(q1)

        Ps = dPc[2:, 2:, 2:, :];  Cs = P123[2:, 2:, 2:, :];
        q2 = np.copy(np.sum(Ps*Cs, axis=(0,1,2))/np.sum(Cs, axis=(0,1,2) ))
        q2 = q2 - np.average(q2)


        #-- plot data --

        ax[0].plot(phi, avg, style)
        ax[1].plot(phi, q0, style)
        ax[2].plot(phi, q1, style)
        ax[3].plot(phi, q2, style, label=label)


#-----------------------------------------------

#ax[0].set_title(title)

ax[0].set_ylim(-1.25,1.25)
ax[1].set_ylim(-1.25,1.25)
ax[2].set_ylim(ylim)
ax[3].set_ylim(ylim)

ax[0].set_title(title + ", average")
ax[1].set_title(title + ", $x < 1$")
ax[2].set_title(title + ", $x > 1$")
ax[3].set_title(title + ", $x > 2$")


ax[0].set(ylabel='$[\, \ln {\cal P} - \ln  {\cal P}_{eq} \,] \, / \, \sqrt{n_1 n_2 n_3 / x_1 x_2 x_3}$ - c')
ax[3].legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))


for i in (0,1,2,3):
    iax=ax[i]
    iax.set_xlim(-1,1)
    iax.set(xlabel = "$\\theta/\pi$")


#iax.set_ylim(0,2.0)
#iax.grid(axis="y")
#iax.legend(frameon=False)
#iax.text(2, 1.50, "direct")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_xscale('log')
#iax.set_yticks(np.arange(0, 2.1, step=0.5)) 


#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.15, wspace=0.3)
plt.subplots_adjust(left=0.05, right=0.90, top=0.9, bottom=0.15, hspace=0.14)

plt.savefig(outfile) #pad_inches=0

plt.close()


#-------------------------------------------------




