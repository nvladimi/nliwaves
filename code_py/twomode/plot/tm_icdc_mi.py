
import numpy as np
import pickle

import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

fnameout = "tm_icdc_mi.pdf"

#xlim = (-3.5, 3.5)
LW  = 1
MS  = 6
LP = -2

#-------------------------------------
                
#-------------------------------------
    
fig, ax = plt.subplots(ncols=2, nrows=2, tight_layout=True, figsize=(8,8))
fig.tight_layout(pad=1)


s1o = 0 #  1/np.log(2)
s2o = 0 # 1/np.log(2) - 1;
s3o = np.log2(2*np.pi);

dc1 = pickle.load(open("../POSTMI/mi3_dc.pkl", 'rb'))
ic1 = pickle.load(open("../POSTMI/mi3_ic.pkl", 'rb'))
dc2 = pickle.load(open("../POSTMI/mi4_dc.pkl", 'rb'))
ic2 = pickle.load(open("../POSTMI/mi4_ic.pkl", 'rb'))
dc3 = pickle.load(open("../POSTMI/mi5_dc.pkl", 'rb'))
ic3 = pickle.load(open("../POSTMI/mi5_ic.pkl", 'rb'))


#-- S1 -------------------------------------


iax = ax[0,0]


iax.plot(dc1["chi"], dc1["S1"] - s1o, '--rs', mfc='none', lw=LW, ms=MS)
iax.plot(ic1["chi"], ic1["S1"] - s1o, '--bs', mfc='none', lw=LW, ms=MS)
iax.plot(dc2["chi"], dc2["S1"] - s1o, ':rd', mfc='none', lw=LW, ms=MS)
iax.plot(ic2["chi"], ic2["S1"] - s1o, ':bd', mfc='none', lw=LW, ms=MS)
iax.plot(dc3["chi"], dc3["S1"] - s1o, '-ro', mfc='none', lw=LW, ms=MS)
iax.plot(ic3["chi"], ic3["S1"] - s1o, '-bo', mfc='none', lw=LW, ms=MS)


#iax.set(ylim = (-1, 1), xlim = xlim)
iax.set(xlabel = "$\chi$", ylabel = "$\Delta S_1$")  
iax.axhline(lw=0.5, c="gray")
iax.yaxis.labelpad = LP  
iax.set_xscale('log')
iax.grid()

ticks = (1e-4, 1e-3, 1e-2, 1, 5, 25, 100, 500)
iax.set_xticks(ticks)
iax.set_xticklabels(["$%g$"%y for y in ticks])



#-- S2 -------------------------------------

iax = ax[0,1]

iax.plot(dc1["chi"], dc1["S2"] - s2o, '--rs', mfc='none', lw=LW, ms=MS)
iax.plot(ic1["chi"], ic1["S2"] - s2o, '--bs', mfc='none', lw=LW, ms=MS)
iax.plot(dc2["chi"], dc2["S2"] - s2o, ':rd', mfc='none', lw=LW, ms=MS)
iax.plot(ic2["chi"], ic2["S2"] - s2o, ':bd', mfc='none', lw=LW, ms=MS)
iax.plot(dc3["chi"], dc3["S2"] - s2o, '-ro', mfc='none', lw=LW, ms=MS)
iax.plot(ic3["chi"], ic3["S2"] - s2o, '-bo', mfc='none', lw=LW, ms=MS)


#iax.set(ylim = (-1, 1), xlim = xlim)
iax.set(xlabel = "$\chi$", ylabel = "$\Delta S_2$")  
iax.axhline(lw=0.5, c="gray")
iax.yaxis.labelpad = LP  
iax.set_xscale('log')
iax.grid()

ticks = (1e-4, 1e-3, 1e-2, 1, 5, 25, 100, 500)
iax.set_xticks(ticks)
iax.set_xticklabels(["$%g$"%y for y in ticks])


#-- S12 ---------------------------------------------------


iax = ax[1,0]


iax.plot(dc1["chi"], dc1["S12"] - s1o - s2o -s3o, '--rs', mfc='none', lw=LW, ms=MS)
iax.plot(ic1["chi"], ic1["S12"] - s1o - s2o -s3o, '--bs', mfc='none', lw=LW, ms=MS)
iax.plot(dc2["chi"], dc2["S12"] - s1o - s2o -s3o, ':rd', mfc='none', lw=LW, ms=MS)
iax.plot(ic2["chi"], ic2["S12"] - s1o - s2o -s3o, ':bd', mfc='none', lw=LW, ms=MS)
iax.plot(dc3["chi"], dc3["S12"] - s1o - s2o -s3o, '-ro', mfc='none', lw=LW, ms=MS)
iax.plot(ic3["chi"], ic3["S12"] - s1o - s2o -s3o, '-bo', mfc='none', lw=LW, ms=MS)

#iax.set(ylim = (-1, 1), xlim = xlim)
iax.set(xlabel = "$\chi$", ylabel = "$\Delta S_{12\\theta}$")  
iax.axhline(lw=0.5, c="gray")
#iax.yaxis.labelpad = LP  
iax.set_xscale('log')

iax.grid()

ticks = (1e-4, 1e-3, 1e-2, 1, 5, 25, 100, 500)
iax.set_xticks(ticks)
iax.set_xticklabels(["$%g$"%y for y in ticks])


#-- I12 ---------------------------------------------------



iax = ax[1,1]


iax.plot(dc1["chi"], dc1["I12"], '--rs', mfc='none', lw=LW, ms=MS, label = "DC: $\Delta n=1$, far eq")
iax.plot(ic1["chi"], ic1["I12"], '--bs', mfc='none', lw=LW, ms=MS, label = "IC:   $\Delta n=1$, far eq")
iax.plot(dc2["chi"], dc2["I12"], ':rd',  mfc='none', lw=LW, ms=MS, label = "DC: $\Delta n=0.5$, far eq" )
iax.plot(ic2["chi"], ic2["I12"], ':bd',  mfc='none', lw=LW, ms=MS, label = "IC:   $\Delta n=0.5$, far eq" )
iax.plot(dc3["chi"], dc3["I12"], '-ro',  mfc='none', lw=LW, ms=MS, label = "DC: $\Delta n=1$,  near eq" )
iax.plot(ic3["chi"], ic3["I12"], '-bo',  mfc='none', lw=LW, ms=MS, label = "IC:   $\Delta n=1$,  near eq" )

#iax.set(ylim = (-1, 1), xlim = xlim)
iax.set(xlabel = "$\chi$", ylabel = "$I_{12}$")  
iax.axhline(lw=0.5, c="gray")
#iax.yaxis.labelpad = LP  
iax.set_xscale('log')
iax.grid()

ticks = (1e-4, 1e-3, 1e-2, 1, 5, 25, 100, 500)
iax.set_xticks(ticks)
iax.set_xticklabels(["$%g$"%y for y in ticks])


#iax.legend(frameon=True, fontsize=7.5, labelspacing=0.3, handlelength=1.5, loc="upper left")
iax.legend()


#-------------------------------------------------------------

plt.savefig(fnameout, pad_inches=0)
plt.close()     
        


#plt.show()
