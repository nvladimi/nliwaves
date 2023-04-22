import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


outfile = "mirror-returned.pdf"

#-- common parameters --

MS = 2
LW = 1
  
fig, ax = plt.subplots(ncols=4, nrows=2, figsize=(10,5))

ls1='-or'
ls2='-og'
ls3='-ob'
#------------------------------------------------------


fname = "../2022-10-15_phi1p5/data/lega.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]
y2 = dat[:,9]
y3 = dat[:,10]
y5a = dat[:,15]
y5b = dat[:,16]

ax[0][0].plot(x,y2, ls1, label = '1.5 cm',  mfc='none', ms=MS,  lw=LW);
ax[0][1].plot(x,y3, ls1, label = '1.5 cm',  mfc='none', ms=MS,  lw=LW);
ax[0][2].plot(x,y5a,ls1, label = '1.5 cm',  mfc='none', ms=MS,  lw=LW);
ax[0][3].plot(x,y5b,ls1, label = '1.5 cm',  mfc='none', ms=MS,  lw=LW);


fname = "../2022-10-15_phi5/data/lega.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]
y2 = dat[:,9]
y3 = dat[:,10]
y5a = dat[:,15]
y5b = dat[:,16]

ax[0][0].plot(x,y2, ls2, label = '5 cm',  mfc='none', ms=MS,  lw=LW);
ax[0][1].plot(x,y3, ls2, label = '5 cm',  mfc='none', ms=MS,  lw=LW);
ax[0][2].plot(x,y5a,ls2, label = '5 cm',  mfc='none', ms=MS,  lw=LW);
ax[0][3].plot(x,y5b,ls2, label = '5 cm',  mfc='none', ms=MS,  lw=LW);


fname = "../2022-10-15_phi15/data/lega.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]
y2 = dat[:,9]
y3 = dat[:,10]
y5a = dat[:,15]
y5b = dat[:,16]

ax[0][0].plot(x,y2, ls3, label = '15 cm',  mfc='none', ms=MS,  lw=LW);
ax[0][1].plot(x,y3, ls3, label = '15 cm',  mfc='none', ms=MS,  lw=LW);
ax[0][2].plot(x,y5a,ls3, label = '15 cm',  mfc='none', ms=MS,  lw=LW);
ax[0][3].plot(x,y5b,ls3, label = '15 cm',  mfc='none', ms=MS,  lw=LW);



fname = "../2022-10-15_phi1p5/data/zero.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]
y2 = dat[:,9]
y3 = dat[:,10]
y5a = dat[:,15]
y5b = dat[:,16]

ax[1][0].plot(x,y2, ls1, label = '1.5 cm',  mfc='none', ms=MS,  lw=LW);
ax[1][1].plot(x,y3, ls1, label = '1.5 cm',  mfc='none', ms=MS,  lw=LW);
ax[1][2].plot(x,y5a,ls1, label = '1.5 cm',  mfc='none', ms=MS,  lw=LW);
ax[1][3].plot(x,y5b,ls1, label = '1.5 cm',  mfc='none', ms=MS,  lw=LW);

fname = "../2022-10-15_phi5/data/zero.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]
y2 = dat[:,9]
y3 = dat[:,10]
y5a = dat[:,15]
y5b = dat[:,16]

ax[1][0].plot(x,y2, ls2, label = '5 cm',  mfc='none', ms=MS,  lw=LW);
ax[1][1].plot(x,y3, ls2, label = '5 cm',  mfc='none', ms=MS,  lw=LW);
ax[1][2].plot(x,y5a,ls2, label = '5 cm',  mfc='none', ms=MS,  lw=LW);
ax[1][3].plot(x,y5b,ls2, label = '5 cm',  mfc='none', ms=MS,  lw=LW);


fname = "../2022-10-15_phi15/data/zero.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]
y2 = dat[:,9]
y3 = dat[:,10]
y5a = dat[:,15]
y5b = dat[:,16]

ax[1][0].plot(x,y2, ls3, label = '15 cm',  mfc='none', ms=MS,  lw=LW);
ax[1][1].plot(x,y3, ls3, label = '15 cm',  mfc='none', ms=MS,  lw=LW);
ax[1][2].plot(x,y5a,ls3, label = '15 cm',  mfc='none', ms=MS,  lw=LW);
ax[1][3].plot(x,y5b,ls3, label = '15 cm',  mfc='none', ms=MS,  lw=LW);



        
#-- canvas options --

ax[0][0].text(1.3, 0.93, "atmosphere", fontsize=12)
ax[1][0].text(1.3, 0.9, "vacuum",  fontsize=12)


for iax in (ax[0][0], ax[0][1], ax[0][2], ax[0][3], ax[1][0], ax[1][1], ax[1][2], ax[1][3]):
    iax.set_xlim(1,100)
    iax.set_ylim(0,1.0)
    iax.set(xlabel='$\\nu$')
    #iax.legend(frameon=False, labelspacing = 0.3, loc="upper right")
    iax.set_xscale('log')
    #iax.set_yscale('log')
    iax.set_xticks((1, 2, 4, 8, 16, 32, 64))
    iax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    iax.minorticks_off()
    
    
for i in (0,1):
    ax[i][0].set(ylabel='$|\Delta I|$')
    ax[i][1].set(ylabel='$P_{\\rm diff}/P$')
    ax[i][2].set(ylabel='$r(0.5P), \, cm$')
    ax[i][3].set(ylabel='$r(0.9P), \, cm$')

ax[0][2].set_ylim(0,30)
ax[0][3].set_ylim(0,100)
ax[1][2].set_ylim(0,20)
ax[1][3].set_ylim(0,60)


def r0(w0, p):
    return(w0 * np.sqrt( -0.5*np.log(1-p) ))

y0=r0(1.5,0.5)
ax[0][2].hlines(y=y0, xmin=1, xmax=100, ls="--", color="r", lw=1)
ax[1][2].hlines(y=y0, xmin=1, xmax=100, ls="--", color="r", lw=1)

y0=r0(5,0.5)
ax[0][2].hlines(y=y0, xmin=1, xmax=100, ls="--", color="g", lw=1)
ax[1][2].hlines(y=y0, xmin=1, xmax=100, ls="--", color="g", lw=1)

y0=r0(15,0.5)
ax[0][2].hlines(y=y0, xmin=1, xmax=100, ls="--", color="b", lw=1)
ax[1][2].hlines(y=y0, xmin=1, xmax=100, ls="--", color="b", lw=1)


y0=r0(1.5,0.9)
ax[0][3].hlines(y=y0, xmin=1, xmax=100, ls="--", color="r", lw=1)
ax[1][3].hlines(y=y0, xmin=1, xmax=100, ls="--", color="r", lw=1)

y0=r0(5,0.9)
ax[0][3].hlines(y=y0, xmin=1, xmax=100, ls="--", color="g", lw=1)
ax[1][3].hlines(y=y0, xmin=1, xmax=100, ls="--", color="g", lw=1)

y0=r0(15,0.9)
ax[0][3].hlines(y=y0, xmin=1, xmax=100, ls="--", color="b", lw=1)
ax[1][3].hlines(y=y0, xmin=1, xmax=100, ls="--", color="b", lw=1)



#for iax in (ax[0][0], ax[1][1]):
ax[0][0].legend(frameon=False, labelspacing = 0.3, loc="lower right")

#ax[1][1].legend(frameon=False, labelspacing = 0.3)


#iax.set_yticks(np.arange(-20, 1, step=5))
#ax[1].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.94), labelspacing = 0.4)
    
#-----------------------

plt.subplots_adjust(left=0.06, right=0.99, top=0.98, bottom=0.10, wspace=0.30, hspace=0.26)


plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
