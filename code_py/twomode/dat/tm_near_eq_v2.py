
import numpy as np
import pickle

import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt


fnameout = "tm_near_eq_v2.pdf"
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
LP = -2

#-------------------------------------
                
#-------------------------------------
    
fig, ax = plt.subplots(ncols=4, nrows=2, tight_layout=True, figsize=(10,5.0))
fig.tight_layout(pad=1)


s1o = 1/np.log(2)
s2o = 1/np.log(2) - 1;
s3o = np.log2(2*np.pi);

#-- S1 -------------------------------------


iax = ax[0,0]
iax.text(2.6, -0.9, '(a)')

iax.set(ylim = (-1, 1), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$", ylabel = "$\Delta S_1$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x  = d["dTT"]
            s1 = d["S1"] - s1o
            iax.plot(x, s1, '-o', mfc='none', color = plotchi[chi], lw=LW, ms=MS, label="$\chi={}$".format(chi))
x = np.linspace(-3.5, 3.5,100);
y = (1 + np.log(1+x/4))/np.log(2) - s1o
iax.plot(x,y, '--k', linewidth=0.7)
iax.axhline(lw=0.5, c="gray")
iax.axvline(lw=0.5, c="gray")
iax.yaxis.labelpad = LP  

iax.legend(frameon=False, fontsize=7.5, labelspacing=0.3, handlelength=1.5, loc="upper left")



iax = ax[1,0]
iax.text(2.6, 4.5e-4, '(e)')

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
iax.axvline(lw=0.5, c="gray")
iax.yaxis.labelpad = LP  


#-- S2 -------------------------------------


iax = ax[0,1]
iax.text(-3.2, -0.9, '(b)')


iax.set(ylim = (-1, 1), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$", ylabel = "$\Delta S_2$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x = d["dTT"]
            s2 = d["S2"] - s2o
            iax.plot(x, s2, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)
x = np.linspace(-3.5, 3.5,100);
y = (1 + np.log(1-x/4))/np.log(2) - 1 - s2o
iax.plot(x,y, '--k', linewidth=0.7)
iax.axhline(lw=0.5, c="gray")
iax.axvline(lw=0.5, c="gray")
iax.yaxis.labelpad = LP  

iax = ax[1,1]
iax.text(2.6, 4.5e-4, '(f)')

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
iax.axvline(lw=0.5, c="gray")
iax.yaxis.labelpad = LP  



#-- S12 ---------------------------------------------------

iax = ax[0,2]
iax.text(2.6, -0.93, '(c)')

iax.set(ylim = (-1, 0.2), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$", ylabel = "$\Delta S_{12\\theta}$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x   = d["dTT"]
            s12 = d["S12"] - s1o -s2o -s3o
            iax.plot(x, s12, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)
iax.axhline(lw=0.5, c="gray")
iax.axvline(lw=0.5, c="gray")

ticks = (-1, -0.6, -0.2, 0, 0.2)
iax.set_yticks(ticks)
iax.set_yticklabels(["$%g$"%y for y in ticks])
iax.yaxis.labelpad = LP   


            
iax = ax[1,2]
iax.text(3.8, 1.8e-4, '(g)')

iax.set_yscale('log')
iax.set_xscale('log')
iax.set(ylim = (1e-4, 0.3), xlim = (0.2,6))
iax.set(xlabel = "$|\Delta T / T|$",
         ylabel = "$\log(2\pi) - S_\\theta$")  
for chi in (0.108, 0.256, 0.5, 4.0):
    for d in alldict:
        if d["chi"] == chi:
            x   =  d["dTT"]
            y   =  s3o - d["Sphi"]
            iax.plot(x, y, '-o', mfc='none', color = plotchi[chi],
                     lw=LW,  ms=MS,  label="$\chi={}$".format(chi))
            iax.plot(-x, y, '--v', mfc='none', color = plotchi[chi],
                     lw=LW,  ms=MS)
iax.yaxis.labelpad = 0  
ticks = (0.5, 1, 2, 4)
iax.set_xticks(ticks)
iax.set_xticklabels(["$%g$"%y for y in ticks])
iax.legend(frameon=False, fontsize=7.5, labelspacing=0.3, handlelength=1.5, loc="upper left")


#-- I12 ---------------------------------------------------

        
iax = ax[0,3]
iax.text(-3.2, 0.008, '(d)')

iax.set(ylim = (0, 0.15), xlim = xlim)
iax.set(xlabel = "$\Delta T / T$",
        ylabel = "$I_{12}$")  
for chi in plotchi:
    for d in alldict:
        if d["chi"] == chi:
            x   =  d["dTT"]
            I12 =  d["I12"]
            iax.plot(x, I12, '-o', mfc='none', color = plotchi[chi],  linewidth=LW,  markersize=MS)
iax.axvline(lw=0.5, c="gray")

ticks = np.arange(0,0.16,0.05)
iax.set_yticks(ticks)
iax.set_yticklabels(["$%g$"%y for y in ticks])
iax.yaxis.labelpad = 0  


            
iax = ax[1,3]
iax.text(3.8, 1.8e-4, '(h)')


iax.set_yscale('log')
iax.set_xscale('log')
iax.set(ylim = (1e-4, 0.3), xlim = (0.2,6))
iax.set(xlabel = "$|\Delta T / T|$",
        ylabel = "$I_{12} - 0.005$")  
for chi in (0.108, 0.256, 0.5, 4.0):
    for d in alldict:
        if d["chi"] == chi:
            x   =  d["dTT"]
            I12 =  d["I12"] - 0.005
            iax.plot(x, I12, '-o', mfc='none', color = plotchi[chi],
                     lw=LW,  ms=MS,  label="$\chi={}$".format(chi))
            iax.plot(-x, I12, '--v', mfc='none', color = plotchi[chi],
                     lw=LW,  ms=MS)
iax.yaxis.labelpad = 0  
ticks = (0.5, 1, 2, 4)
iax.set_xticks(ticks)
iax.set_xticklabels(["$%g$"%y for y in ticks])
iax.legend(frameon=False, fontsize=7.5, labelspacing=0.3, handlelength=1.5, loc="upper left")


#(0.108, 0.256, 0.5, 2.048, 4.0, 8.788)



#-------------------------------------------------------------

plt.savefig(fnameout, pad_inches=0)


plt.close()     
        
