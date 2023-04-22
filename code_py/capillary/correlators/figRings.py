
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import numpy as np
fnameout = "figRings.pdf"

w1 = np.array([170, 272, 485, 680])
r1 = np.array([1, 1, 1, 1])
r2 = np.array([2.17, 1.82, 1.52, 1.38])
r3 = np.array([2.90, 2.38, 1.90, 1.70])
r4 = np.array([3.50, 2.85, 2.25, 1.95])

#-- images --

fig, iax = plt.subplots(ncols=1, nrows=1, constrained_layout=True, figsize=(3,8.5))

iax.plot(w1, w1,       'ok', label='1st circle')
iax.plot(w1, w1*r2**2, 'or', label='2nd circle')
iax.plot(w1, w1*r3**2, 'og', label='3rd circle')
iax.plot(w1, w1*r4**2, 'ob', label='4th circle')
iax.plot(680, 50, 'sk', label='inner circle')

w=np.arange(0,900,100)
iax.plot(w, w,      '--k', linewidth=1, label='$\omega_1$')
iax.plot(w, w+630,  '--r', linewidth=1, label='$\omega_1 + 630$')
iax.plot(w, w+1260, '--g', linewidth=1, label='$\omega_1 + 1260$')
iax.plot(w, w+1900, '--b', linewidth=1, label='$\omega_1 + 1900$')
iax.plot(w, w-630,  ':k',  linewidth=1, label='$\omega_1 - 630$')


print(w)

#iax.plot(k, nk, marker=markers[i], color=colors[i], linestyle=lines[i],
#             markersize=1.5, linewidth=0.7, label=lbl[i])


  
#iax.plot(k, 0.01*k**(-4), '--k', linewidth=0.7, label="$k^{-4}$")

iax.set_xlabel('$\omega_1$')
iax.set_ylabel('$\omega_1, \quad   \omega_2, \quad \omega_3, \quad \omega_4$')

iax.set_xlim((0,800))
iax.set_ylim((0,3200))

#ticks = 2**(np.arange(0,8))
#iax.set_xticks(ticks)
#iax.set_xticklabels(["$%g$"%y for y in ticks])


iax.legend(frameon=False) # loc="upper right")


plt.savefig(fnameout, pad_inches=0)



#plt.show()
plt.close()








