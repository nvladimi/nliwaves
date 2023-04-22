
import numpy as np
import pickle

import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt


fnameout = "tm_near_eq.pdf"
datpkl =   "POST/post05_mi.pkl"

with open(datpkl, 'rb') as handle:
    alldict = pickle.load(handle)

print(alldict[0].keys())

plotchi = {0.108: "#FF0000",
           0.256: "blue",
           0.5:   "#60B000",
           2.048: "#FF8000",
           4.0:   "#800080",
           8.788: "#008080"}

xlim = (-3.5, 3.5)
LW  = 1
MS  = 3

#-------------------------------------

if False:
    for chi in plotchi:
        for d in alldict:
            if d["chi"] == chi:
                print(chi)
                q = np.vstack((d["dTT"], d["n1s"],  d["n2s"],  d["ntot"]))
                print(q.transpose())
                
#-------------------------------------
    
fig, ax = plt.subplots(ncols=4, nrows=4, tight_layout=True, figsize=(12.5,12))



#-- S1 --

iax = ax[0,0]
iax.grid()
iax.set(ylim = (-0.5, 2.5), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$", ylabel = "$S_1, \quad S(n_1)$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x  = d["dTT"]
            n1 = d["n1s"]
            s1 = d["S1"]
            y1 = 1/np.log(2) * (1 + np.log(n1));
            iax.plot(x, s1, 'o', mfc='none', color = plotchi[chi], linewidth=LW, markersize=MS)
            iax.plot(x, y1, '--', color = plotchi[chi],  linewidth=LW)
x = np.linspace(-3.5, 3.5,100);
y = (1 + np.log(1+x/4))/np.log(2)
iax.plot(x,y, '-k', linewidth=0.5)

iax = ax[1,0]
iax.grid()
iax.set(ylim = (-0.08, 0.02), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$", ylabel = "$S_1 - S(n_1)$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x  = d["dTT"]
            n1 = d["n1s"]
            s1 = d["S1"]
            y1 = 1/np.log(2) * (1 + np.log(n1));
            iax.plot(x, s1 - y1, '-o', mfc='none',  linewidth=LW,  markersize=MS,
                     color = plotchi[chi], label = "$\chi={}$".format(chi))
iax.legend()


iax = ax[2,0]
iax.grid()
iax.set_yscale('log')
iax.set(ylim = (3e-4, 0.1), xlim = xlim)
iax.set(xlabel = "$\Delta T/T$", ylabel = "$S(n_1) - S_1$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x = d["dTT"]
            n1 = d["n1s"]
            s1 = np.array(d["S1"])
            y1 = 1/np.log(2) * (1 + np.log(n1));
            iax.plot(x, y1-s1, '-o', mfc='none',  linewidth=LW,  markersize=MS,
                     color = plotchi[chi], label = "$\chi={}$".format(chi))

iax = ax[3,0]
iax.grid()
iax.set_yscale('log')
iax.set(ylim = (3e-4, 0.1), xlim = xlim)
iax.set(xlabel = "$\\frac{\Delta T}{T}(\\frac{\chi}{1+\chi})^{1/3}$", ylabel = "$S(n_1) - S_1$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x = d["dTT"] * np.power(chi/(1+chi),1/3.)
            n1 = d["n1s"]
            s1 = np.array(d["S1"])
            y1 = 1/np.log(2) * (1 + np.log(n1));
            iax.plot(x, y1-s1, '-o', mfc='none',  linewidth=LW,  markersize=MS,
                     color = plotchi[chi], label = "$\chi={}$".format(chi))


#-- S2 -------------------------------------

iax = ax[0,1]
iax.grid()
iax.set(ylim = (-1.5, 1.5), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$", ylabel = "$S_2, \quad S(n_2)$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x = d["dTT"]
            n2 = d["n2s"]
            s2 = d["S2"]
            y2 = 1/np.log(2) * (1 + np.log(n2));
            iax.plot(x, s2, 'o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)
            iax.plot(x, y2, '--', color = plotchi[chi],  linewidth=LW)
x = np.linspace(-3.5, 3.5,100);
y = (1 + np.log(1-x/4))/np.log(2) - 1
iax.plot(x,y, '-k', linewidth=0.5)


iax = ax[1,1]
iax.grid()
iax.set(ylim = (-0.08, 0.02), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$", ylabel = "$S_2 - S(n_2)$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x = d["dTT"]
            n2 = d["n2s"]
            s2 = d["S2"]
            y2 = 1/np.log(2) * (1 + np.log(n2));
            iax.plot(x, s2 - y2, '-o', mfc='none',  linewidth=LW,  markersize=MS,
                     color = plotchi[chi], label = "$\chi={}$".format(chi) )
iax.legend()

iax = ax[2,1]
iax.grid()
iax.set_yscale('log')
iax.set(ylim = (3e-4, 0.1), xlim = xlim)
iax.set(xlabel = "$\Delta T/T$", ylabel = "$S(n_2) - S_2$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x = d["dTT"]
            n2 = d["n2s"]
            s2 = np.array(d["S2"])
            y2 = 1/np.log(2) * (1 + np.log(n2));
            iax.plot(x, y2-s2, '-o', mfc='none',  linewidth=LW,  markersize=MS,
                     color = plotchi[chi], label = "$\chi={}$".format(chi))

iax = ax[3,1]
iax.grid()
iax.set_yscale('log')
iax.set(ylim = (3e-4, 0.1), xlim = xlim)
iax.set(xlabel = "$\\frac{\Delta T}{T}(\\frac{\chi}{1+\chi})^{1/2}$", ylabel = "$S(n_2) - S_2$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x = d["dTT"] * np.power(chi/(1+chi),1/2.)
            n2 = d["n2s"]
            s2 = np.array(d["S2"])
            y2 = 1/np.log(2) * (1 + np.log(n2));
            iax.plot(x, y2-s2, '-o', mfc='none',  linewidth=LW,  markersize=MS,
                     color = plotchi[chi], label = "$\chi={}$".format(chi))



#-- Sphi --------------------------------------------------

iax = ax[0,2]
iax.grid()
iax.set(ylim = (2.54, 2.66), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$", ylabel = "$S_\\theta$")

y3 = np.log2(2*np.pi);
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x   = d["dTT"]
            sphi = d["Sphi"]
            iax.plot(x, sphi, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)
iax.plot(x, x*0 + y3, '--', color = "black",  linewidth=LW)                

iax = ax[1,2]
iax.grid()
iax.set(ylim = (-0.1, 0.02), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$",
        ylabel = "$S_\\theta - \log(2\pi)$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x   = d["dTT"]
            sphi = d["Sphi"]
            iax.plot(x, sphi-y3, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)


iax = ax[2,2]
iax.set_yscale('log')
iax.set_xscale('log')
iax.grid()
iax.set(ylim = (5e-5, 0.2), xlim=(0.1,10) )
iax.set(xlabel = "$|\Delta T / T|$",
        ylabel = "$\log(2\pi) - S_\\theta$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x   = np.abs(d["dTT"])
            sphi = d["Sphi"]
            iax.plot(x, y3-sphi, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)


iax = ax[3,2]
iax.grid()
iax.set(ylim = (-0.1, 0.02), xlim=xlim )
iax.set(xlabel =  "$\\frac{\Delta T}{T}(\\frac{\chi}{1+\chi})^{1/3}$",
        ylabel = "$S_\\theta -  \log(2\pi)$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x   = d["dTT"]*np.power(chi/(1+chi),1/3.)
            sphi = d["Sphi"]
            iax.plot(x, sphi-y3, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)






#-- S12 ---------------------------------------------------

iax = ax[0,3]
iax.grid()
iax.set(ylim = (3.4, 4.6), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$", ylabel = "$S_{12}, \quad S(n_1)+S(n_2)+\log(2\pi) $")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x   = d["dTT"]
            s12 = d["S12"]
            s1  = d["S1"]
            s2  = d["S2"]
            sphi = d["Sphi"]
            n1  = d["n1s"]
            n2  = d["n2s"]
            y1  = 1/np.log(2) * (1 + np.log(n1)) ;
            y2  = 1/np.log(2) * (1 + np.log(n2)) ;
            y3  = np.log2(2*np.pi);
            iax.plot(x, s12, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)
            iax.plot(x, y1+y2+y3, '--', color = plotchi[chi],  linewidth=LW)                

iax = ax[1,3]
iax.grid()
iax.set(ylim = (0, 0.14), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$",
        ylabel = "$I_{12}, \quad S(n_1)+S(n_2)+\log(2\pi) - S_{12} $")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x   =  d["dTT"]
            s12 =  d["S12"]
            s1  =  d["S1"]
            s2  =  d["S2"]
            sphi = d["Sphi"]
            n1 =   d["n1s"]
            n2 =   d["n2s"]
            y1 = 1/np.log(2) * (1 + np.log(n1)) ;
            y2 = 1/np.log(2) * (1 + np.log(n2)) ;
            y3 = np.log2(2*np.pi);
            iax.plot(x, s1+s2+sphi-s12, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)
            iax.plot(x, y1+y2+y3-s12, '--', color = plotchi[chi],  linewidth=LW)                

iax = ax[2,3]
iax.grid()
iax.set_yscale('log')
iax.set_xscale('log')
iax.set(ylim = (1e-4, 1), xlim = (0.1,10))
iax.set(xlabel = "$|\Delta T / T|$",
        ylabel = "$I_{12} - 0.005$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x   =  d["dTT"]
            I12 =  d["I12"] - 0.005
            iax.plot(x, I12, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)


iax = ax[3,3]
iax.grid()
iax.set(ylim = (0, 0.12), xlim = xlim)
iax.set(xlabel = "$\\frac{\Delta T}{T}(\\frac{\chi}{1+\chi})^{1/3}$",
        ylabel = "$I_{12}$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x   =  d["dTT"] * np.power(chi/(1+chi), 1/3.)
            I12 =  d["I12"]
            iax.plot(x, I12, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)




plt.savefig(fnameout, pad_inches=0)


plt.close()     
        
