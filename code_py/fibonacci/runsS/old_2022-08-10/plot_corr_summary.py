import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

datafile="all_corr.txt" ; outfile = "plot_corr_summary.pdf"
datafile="all_corr_cherry.txt" ; outfile = "plot_corr_cherry.pdf"


# all_corr.txt:
#
#  0.name   1.2.3.4.s5m90(I)    5.6.7.8.s5m10(D)
#              9.10.s4v38(I)      11.12.s4v58(D)
#             13.14.s4v14(I)      15.16.s4v34(D)
#  17.power  18.sum  19.count
#
#  Data is recoderd in pairs slope/intercept.
#  For alpha=1/2 two pairs are provided with
#  fitted and theoretical slopes.


from scipy.special import gamma


fig, ax = plt.subplots(ncols=4, nrows=2, figsize=(12,6))

ax[0][2].axis("off")

lt1='or'
lt2='db'

#-- data --

dat = np.loadtxt(datafile)

m = dat[:,17]
s12i = dat[:,1]
s12d = dat[:,5]
s12t = dat[:,3]
r12i = dat[:,2]
r12d = dat[:,6]
r12it = dat[:,4]
r12dt = dat[:,8]

p12i =  dat[:,2] + np.log(190)*dat[:,1]
p12d =  dat[:,6] + np.log(190)*dat[:,5]


p18i =  dat[:,10] + 90*dat[:,9]
p18d = -dat[:,12] - 10*dat[:,11]
p14i =  dat[:,14] + 90*dat[:,13]
p14d = -dat[:,16] - 10*dat[:,15]

#-- theoretical --

mt=np.arange(3,10)

p18t = gamma(mt/2 + 1 ) * (19)**(mt/3) * (1/8)
p14t = gamma(mt/2 + 1 ) * (9)**(mt/3) * (1/4)
p34t = gamma(mt/2 + 1 ) * (10)**(mt/3) * (1/4)



#-- plots --

iax=ax[0][0]
iax.plot(m, s12i, lt1,  mfc='none', ms=5,  lw=1, label="inverse")
iax.plot(m, s12d, lt2,  mfc='none', ms=5,  lw=1, label="direct")
iax.plot(m, s12t, '-k',  mfc='none', lw=0.5, label="theory")
iax.set(xlabel='$m$')
#iax.legend(frameon=False)
iax.set(title='fit $\ln C_m(i)=s\ln i + r^\prime$ for $\\delta = 0$')
iax.set(ylabel='$s$')

iax=ax[0][1]
iax.plot(m, r12it, lt1,  mfc='none', ms=5,  lw=1, label="inverse")
iax.plot(m, r12dt, lt2,  mfc='none', ms=5,  lw=1, label="direct")
iax.set(xlabel='$m$')
#iax.legend(frameon=False)
iax.set(title='fit $\ln C_m(i)=(m/3 - 1)\ln i + r$')
iax.set(ylabel='$r$')

iax=ax[0][3]
iax.plot(m, p12i, lt1,  mfc='none', ms=5,  lw=1, label="inverse")
iax.plot(m, p12d, lt2,  mfc='none', ms=5,  lw=1, label="direct")
iax.set(xlabel='$m$')
#iax.legend(frameon=False)
iax.set(ylabel='$r$')
iax.set(ylabel='$\ln|C_m(p)|$ for $\\delta = 0$')





iax=ax[1][0]
iax.plot(m, p18i, lt1,  mfc='none', ms=5,  lw=1, label="inverse")
iax.plot(m, p18d, lt2,  mfc='none', ms=5,  lw=1, label="direct")
iax.plot(mt, p18t, '-k', lw=0.5, label="$C_m^{(th)}$")
iax.plot(mt, p18t/p18t[0], ':k', lw=1, label="$C_m^{(th)}/C_3^{(th)}$")
iax.set(xlabel='$m$')
#iax.legend(frameon=False)
iax.set(ylabel='$|C_m(p)|$ for $\\delta = \pm 1/8$')
iax.set_yscale('log')
iax.set_ylim(0.5, 5e5)

iax=ax[1][1]
iax.plot(m, p14i, lt1,  mfc='none', ms=5,  lw=1, label="inverse")
iax.plot(m, p14d, lt2,  mfc='none', ms=5,  lw=1, label="direct")
iax.plot(mt, p14t, '-k', lw=0.5, label="$C_m^{(th)}$")
iax.plot(mt, p34t, '-k', lw=0.5)
iax.plot(mt, p14t/p14t[0], ':k', lw=1, label="$C_m^{(th)}/C_3^{(th)}$")
iax.plot(mt, p34t/p34t[0], ':k', lw=1)
iax.set(xlabel='$m$')
#iax.legend(frameon=False)
iax.set(ylabel='$|C_m(p)|$ for $\\delta = \pm 1/4$')
iax.set_yscale('log')
iax.set_ylim(0.5, 5e4)

iax=ax[1][2]
iax.plot(p12i, np.log(p18i), lt1,  mfc='none', ms=5,  lw=1, label="inverse")
iax.plot(p12d, np.log(p18d), lt2,  mfc='none', ms=5,  lw=1, label="direct")
iax.set(xlabel='$\ln|C_m(p)|$ for $\\alpha=1/2$' )
iax.set(ylabel='$\ln|C_m(p)|$ for $\\delta = \pm 1/8$' )

iax=ax[1][3]
iax.plot(p12i, np.log(p14i), lt1,  mfc='none', ms=5,  lw=1, label="inverse")
iax.plot(p12i, np.log(p14d), lt2,  mfc='none', ms=5,  lw=1, label="direct")
iax.set(xlabel='$\ln|C_m(p)|$ for $\\alpha=1/2$' )
iax.set(ylabel='$\ln|C_m(p)|$ for $\\delta = \pm 1/4$' )



#---------------------------------







#Axes.ticklabel_format(self, *, axis='both', style='', scilimits=None, 



#-----------------------

plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.07, hspace=0.28, wspace=0.25)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
