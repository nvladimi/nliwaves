import numpy as np
import fiboMI as mi
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

#------------------------------

outfile='plot_phase_alpha12simple.pdf' 

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

        P123 = np.fromfile(fbase + "_P123.dat", 'int32')

        #-- normalize and reshape probabilities --

        P123 = P123.reshape(numbins, numbins, numbins, numbinsPhi)
        
        phi =  (np.arange(-numbinsPhi/2, numbinsPhi/2) + 0.5 )/numbinsPhi*2

        Pphi = np.copy(np.sum(P123, axis=(0,1,2))) / ncount * numbinsPhi -1

        #-- plot data --

        iax.plot(phi, Pphi, style, label=label)


#-----------------------------------------------

#ax[0].set_title(title)


ax[0].set(ylabel='${\cal P} / {\cal P}_{eq} - 1$')
ax[0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0.02, 0.96), labelspacing = 0.2)
ax[0].text(-0.85, 0.045, "$|i-d|$")
ax[0].text(0.4, -0.035, "inverse", fontsize=12)
ax[1].text(-0.85, -0.035, "direct", fontsize=12)

for i in (0,1):
    iax=ax[i]
    iax.set_xlim(-1,1)
    iax.set(xlabel = "$\\theta/\pi$")
ax[0].set_ylim(-0.04,0.05)
ax[1].set_ylim(-0.04,0.05)
    #iax.set_title(title[i])


plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.15, wspace=0.23)

plt.savefig(outfile) #pad_inches=0

plt.close()


#-------------------------------------------------




