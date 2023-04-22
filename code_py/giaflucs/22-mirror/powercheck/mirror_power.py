import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


outfile = "mirror_power.pdf"

#-- common parameters --

MS = 0.5
LW = 0.5

runs = (
    ("01", 'r', "$\\nu =  1$"),
    ("02", 'y', "$\\nu =  2$"),
    ("04", 'g', "$\\nu =  4$"),
    ("08", 'b', "$\\nu =  8$"),
)

dx = 0.27
dd = dx*dx
w0 = 1.5
        
#---------------------------------



fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(12,3))


for rset in runs:
    run = rset[0]
    label = rset[2]
    color = rset[1]
    
    fname = "data/lega" + run + "_last.txt" 
    print(fname)
    
    dat = np.loadtxt(fname,  comments='%')
        
    r = dat[:,5]
    p = 1-dat[:,6]/dat[-1,6]
    
    
    ax[0].plot(r,p, color=color, mfc='none', ms=MS,  lw=LW, label=label);
    ax[1].plot(r,p, color=color, mfc='none', ms=MS,  lw=LW, label=label);

    
    fname = "data/zero" + run + "_last.txt" 
    print(fname)
    
    dat = np.loadtxt(fname,  comments='%')
        
    r = dat[:,5]
    p = 1-dat[:,6]/dat[-1,6]
    
    ax[2].plot(r,p, color=color,  mfc='none', ms=MS,  lw=LW,  label=label);
    ax[3].plot(r,p, color=color,  mfc='none', ms=MS,  lw=LW,  label=label);

    

x=np.arange(120)


#-- canvas options --

for j in (0,2):
    ax[j].set_xlim(0,5)   
    ax[j].set_ylim(0,1)
    ax[j].hlines(y=0.5, xmin=0, xmax=5, colors='k', ls=':', lw=LW)
    ax[j].hlines(y=0.1, xmin=0, xmax=5, colors='k', ls=':', lw=LW)
    ax[j].plot(r, np.exp(-2*(r/w0)**2), '--k', lw=LW) 
    ax[j].grid(axis='x')    
    ax[j].set_xlabel('r, cm')
    
for j in (1,3):
    ax[j].set_yscale('log')
    ax[j].set_xlim(0,100)
    ax[j].set_ylim(0.01,1)
    ax[j].hlines(y=0.5, xmin=0, xmax=100, colors='k', ls=':', lw=LW)
    ax[j].hlines(y=0.1, xmin=0, xmax=100, colors='k', ls=':', lw=LW)
    ax[j].grid(axis='x')    
    ax[j].set_xlabel('r, cm')

ax[2].legend(handlelength=3, framealpha=1)
#ax[2].legend(frameon=False, handlelength=2, framealpha=1) #facecolor='y', framealpha=1)
ax[0].set_ylabel("power outside circle")

ax[0].set_title("atmosphere, 1.5cm, lin scale")
ax[1].set_title("atmosphere, 1.5cm, log scale")
ax[2].set_title("vaccuum, 1.5cm, lin scale")
ax[3].set_title("vaccuum, 1.5cm, log scale")

#iax.set_yticks(np.arange(-20, 1, step=5))
#ax[1].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.94), labelspacing = 0.4)


    
#-----------------------

plt.subplots_adjust(left=0.09, right=0.98, top=0.9, bottom=0.15, wspace=0.25, hspace=0.25)


plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
