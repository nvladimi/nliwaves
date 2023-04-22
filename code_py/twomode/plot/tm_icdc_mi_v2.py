
import numpy as np
import pickle

import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

fnameout = "tm_icdc_mi_v2.pdf"

#xlim = (-3.5, 3.5)
LW  = 1
MS  = 5
LP = -2

#-------------------------------------
                
#-------------------------------------
    
fig, ax = plt.subplots(ncols=1, nrows=1, tight_layout=True, figsize=(2.5,2.55))
fig.tight_layout(pad=-1)



d = pickle.load(open("../POSTMI/mi3_dc.pkl", 'rb')); x1d=d["chi"]; y1d=d["I12"]
d = pickle.load(open("../POSTMI/mi3_ic.pkl", 'rb')); x1i=d["chi"]; y1i=d["I12"]
d = pickle.load(open("../POSTMI/mi4_dc.pkl", 'rb')); x2d=d["chi"]; y2d=d["I12"]
d = pickle.load(open("../POSTMI/mi4_ic.pkl", 'rb')); x2i=d["chi"]; y2i=d["I12"]
d = pickle.load(open("../POSTMI/mi5_dc.pkl", 'rb')); x3d=d["chi"]; y3d=d["I12"]
d = pickle.load(open("../POSTMI/mi5_ic.pkl", 'rb')); x3i=d["chi"]; y3i=d["I12"]


#-- combine "1" and "3" --

x13d = np.hstack((x1d, x3d))
y13d = np.hstack((y1d, y3d))

ind = np.argsort(x13d)
x13d = x13d[ind]
y13d = y13d[ind]


x13i = np.hstack((x1i, x3i))
y13i = np.hstack((y1i, y3i))

ind = np.argsort(x13i)
x13i = x13i[ind]
y13i = y13i[ind]


#-- take only upper half of "2" --

ind = x2d > 0.5
x20d = x2d[ind]
y20d = y2d[ind]

ind = x2i > 0.5
x20i = x2i[ind]
y20i = y2i[ind]




#-- I12 ---------------------------------------------------



iax = ax


iax.plot(x13d, y13d, '-ro',  mfc='none', lw=LW, ms=MS, label = "Direct:  \t $\;\Delta \\rho^2 /n  = 1$" )
iax.plot(x13i, y13i, '-bo',  mfc='none', lw=LW, ms=MS, label = "Inverse: \t $\Delta \\rho^2 /n =1$" )
iax.plot(x20d, y20d, ':rd',  mfc='none', lw=LW, ms=MS, label = "Direct: \t $\;\Delta \\rho^2 / n = 0.5$" )
iax.plot(x20i, y20i, ':bd',  mfc='none', lw=LW, ms=MS, label = "Inverse: \t $\Delta \\rho^2 /n = 0.5$" )



iax.set(ylim = (0, 4.5))
iax.set(xlabel = "$\chi$", ylabel = "$I_{12}$")  
iax.axhline(lw=0.5, c="gray")
#iax.yaxis.labelpad = LP  
iax.set_xscale('log')
#iax.grid()

#ticks = (1e-4, 1e-3, 1e-2, 1, 5, 25, 100, 500)
#iax.set_xticks(ticks)
#iax.set_xticklabels(["$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "$1$", "$5$", "$25$", "$100$", "$500$"])

ticks = (1e-4, 1e-2, 1, 100)
iax.set_xticks(ticks)
iax.set_xticklabels(["$10^{-4}$", "$10^{-2}$", "$1$",  "$10^{2}$"])


ticks = range(5)
iax.set_yticks(ticks)
iax.set_yticklabels(["$%g$"%y for y in ticks])


#iax.legend(frameon=True, fontsize=7.5, labelspacing=0.3, handlelength=1.5, loc="upper left")
#iax.legend(frameon=False)


#-------------------------------------------------------------

plt.subplots_adjust(left=0.03, right=0.99, top=0.98, bottom=0.09)
#plt.tight_layout()


plt.savefig(fnameout, pad_inches=0)
plt.close()     
        


#plt.show()
