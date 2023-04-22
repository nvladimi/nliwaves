import numpy as np
import fiboMI as mi
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

#------------------------------

outfile='plot_phase_alpha12.pdf' 

runs= ("MIx32x1/s5m90", "MIx32x1/s5m10")


modes= ((("20", "$\,20$", "-r"), 
         ("40", "$\,40$", "-b"), 
         ("80", "$\,80$", "--c"), 
         ("120","$120$","--y"), 
         ("160","$160$", ":m"), 
         ("180","$180$", ":g")), 
        (("180", " $20$", "-r"), 
         ("160", " $40$", "-b"), 
         ("120", " $80$", "--c"), 
         ("80", "$120$", "--y"), 
         ("40", "$160$", ":m"), 
         ("20", "$180$", ":g")) 
)
 


MS = 1
LW = 1
pi = np.pi

fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,2.8))


#-----------------------------


for i in (0,1):

    iax = ax[i]

    run = runs[i]

    for r in modes[i]:

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


        #-- plot data --

        iax.plot(phi, avg, style, label=label)


#-----------------------------------------------

#ax[0].set_title(title)


ax[0].set(ylabel='$[\, \ln {\cal P} - \ln  {\cal P}_{eq} \,] \, / \, \sqrt{n_1 n_2 n_3 / x_1 x_2 x_3} - c$')
ax[0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0.02, 0.94), labelspacing = 0.4)
ax[0].text(-0.65, 1.05, "$|i-d|$")
ax[0].text(0.4, -1.05, "inverse", fontsize=12)
ax[1].text(-0.85, -1.05, "direct", fontsize=12)

for i in (0,1):
    iax=ax[i]
    iax.set_xlim(-1,1)
    iax.set(xlabel = "$\\theta/\pi$")
    iax.set_ylim(-1.25,1.25)
    #iax.set_title(title[i])


plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.15, wspace=0.23)

plt.savefig(outfile) #pad_inches=0

plt.close()


#-------------------------------------------------




