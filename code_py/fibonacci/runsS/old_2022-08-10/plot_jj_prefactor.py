import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

outfile = "plot_jj_prefactor.pdf"



#-- canvas options --

fig, iax = plt.subplots(ncols=1, nrows=1, figsize=(4,4))


dat = np.loadtxt("jj-prefactor.txt")
x = dat[:,0];
y = dat[:,1]
iax.plot(x, y, 'or:',  mfc='none', ms=5,  lw=1, label="$f(j)$")

x = np.arange(0,9,0.1);
iax.plot(x, 2.3*x**(-2), '-g',  mfc='none', ms=5,  lw=1,   label="$2.3 / j^2$")
iax.plot(x, 2*np.exp(-0.5*x), '-b',  mfc='none', ms=5,  lw=1, label="$2 \exp(-j/2)$")

iax.set_yscale('log')


iax.set_title("$\\alpha = 1/2: \quad D(i,j) = 1 + f(j)(i-d)$")

iax.set_xlim(0.5, 8.5)
iax.set_ylim(0.01, 10)
iax.set(xlabel='$j$', ylabel="$f(j)$")
iax.legend(frameon=False)
iax.text(1,0.02,
         "$D = - \phi^{2(1+\\alpha)} \langle b_{i-j-2} b_{i-j-1} b^*_{i-j} b_{i-2}b_{i-1}b^*_i \\rangle $")


#Axes.ticklabel_format(self, *, axis='both', style='', scilimits=None, 



#-----------------------

plt.subplots_adjust(left=0.15, right=0.97, top=0.9, bottom=0.12, hspace=0.25, wspace=0.28)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
