
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import numpy as np
import nlwTools as nlw


fnames = ("../2020-01-06_klo/postdata/a1n256psi_nk.npy",
          "../2020-01-06_klo/postdata/a2n256psi_nk.npy",
          "../2020-01-06_klo/postdata/a5n256psi_nk.npy",
          "../wexac/a1n512_20210123/n512psi_nk.npy")

colors = ("green", "blue", "red", "green")
lines = ("", "solid", "solid", "solid")
markers = ("o", "o", "o", "")

lbl = ("$f = 0.001, \, 256^2$",
       "$f = 0.002, \, 256^2$",
       "$f = 0.005, \, 256^2$",
       "$f = 0.001, \, 512^2$")

#-- images --

fig, iax = plt.subplots(ncols=1, nrows=1, constrained_layout=True, figsize=(3.5,3.25))

for i in range(4):
    aa = np.load(fnames[i])
    nk, k =  nlw.nkInterpolate(aa, dNk=1, Nphi=120, AngleAvg = True, iShift = False)
    iax.plot(k, nk, marker=markers[i], color=colors[i], linestyle=lines[i],
             markersize=1.5, linewidth=0.7, label=lbl[i])



    
iax.plot(k, 0.01*k**(-4), '--k', linewidth=0.7, label="$k^{-4}$")
iax.set_xscale('log')
iax.set_yscale('log', basey=10)

iax.set_xlabel('$|k|$')
iax.set_ylabel('$n_k$')
iax.set_xlim((1,128))
iax.set_ylim((1e-12, 0.1))

ticks = 2**(np.arange(0,8))
iax.set_xticks(ticks)
iax.set_xticklabels(["$%g$"%y for y in ticks])


iax.legend(frameon=False, loc="upper right")

plt.savefig('figAvgSpc.pdf', pad_inches=0)



#plt.show()
plt.close()








