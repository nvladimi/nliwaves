# this is a comment

import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import kinwave as kw

outfile = "plotP5.pdf"

#-- input parameters -- 

q = np.array((-1, 2))
#r = np.array((0, 1)); c1=1.7; c2=28; 
r = np.array((0, 1e-5)); c1=0.4e-4; c2=5e-3; 


#-- create linear arrays for plotting --

P = np.arange(0.01,1.01,0.01)
L = np.zeros(P.shape)
Wpqkr = np.zeros(P.shape)
Zpqkr = np.zeros(P.shape)
Zqpkr = np.zeros(P.shape)
Zkrpq = np.zeros(P.shape)
Zrkpq = np.zeros(P.shape)

i=0
for p0 in P:
    p = np.array((p0, -1))
    k = kw.kfun(p,q,r)
    L[i] = kw.Lfun(p,q,k,r)
    
    Wpqkr[i] = kw.Wfun(p,q,k,r)
    Zpqkr[i] = kw.Zfun(p,q,k,r)
    Zqpkr[i] = kw.Zfun(q,p,k,r)
    Zkrpq[i] = kw.Zfun(k,r,p,q)
    Zrkpq[i] = kw.Zfun(r,k,p,q)

    i = i+1

#-- create log arrays for plotting --

aP =np.geomspace(1e-3, 1e5, num=100, endpoint=True)
aL = np.zeros(aP.shape)

aWpqkr = np.zeros(aP.shape)
aZpqkr = np.zeros(aP.shape)
aZqpkr = np.zeros(aP.shape)
aZkrpq = np.zeros(aP.shape)
aZrkpq = np.zeros(aP.shape)

i=0
for p0 in aP:
    p = np.array((p0, -1))
    k = kw.kfun(p,q,r)
    aL[i] = np.abs(kw.Lfun(p,q,k,r))

    aWpqkr[i] = np.abs(kw.Wfun(p,q,k,r))
    aZpqkr[i] = np.abs(kw.Zfun(p,q,k,r))
    aZqpkr[i] = np.abs(kw.Zfun(q,p,k,r))
    aZkrpq[i] = np.abs(kw.Zfun(k,r,p,q))
    aZrkpq[i] = np.abs(kw.Zfun(r,k,p,q))

    
    i = i+1


#-- plot L(R) in linear and log-log coordinates --    

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(8,8))

iax=ax[0,0]
iax.plot(P, L, '-r', mfc='none', ms=4, lw=1, label="$\lambda$")

iax=ax[0,1]
iax.plot(P, Wpqkr*4, '-r', mfc='none', ms=4, lw=1, label="$4W_{pqkr}$")
iax.plot(P, Zpqkr, '-y', mfc='none', ms=4, lw=1, label="$Z_{pqkr}$")
iax.plot(P, Zqpkr, '--b', mfc='none', ms=4, lw=1, label="$Z_{qpkr}$")
iax.plot(P, Zkrpq, '-g', mfc='none', ms=4, lw=1, label="$Z_{krpq}$")
iax.plot(P, Zrkpq, '-c', mfc='none', ms=4, lw=1, label="$Z_{rkpq}$")


iax=ax[1,0]
iax.plot(aP, aL, '-r', mfc='none', ms=4, lw=2, label="$|\lambda|$")
iax.plot(aP, c1*aP, ':k', mfc='none', ms=4, lw=2, label="{} ".format(c1)+"$b$")
iax.set_yscale('log')
iax.set_xscale('log')

iax=ax[1,1]
iax.plot(aP, aWpqkr*4, '-r', mfc='none', ms=4, lw=1, label="$|4W_{pqkr}|$")
iax.plot(aP, aZpqkr, '-y',   mfc='none', ms=4, lw=1, label="$|Z_{pqkr}|$")
iax.plot(aP, aZqpkr, '--b',  mfc='none', ms=4, lw=1, label="$|Z_{qpkr}|$")
iax.plot(aP, aZkrpq, '-g',   mfc='none', ms=4, lw=1, label="$|Z_{krpq}|$")
iax.plot(aP, aZrkpq, '-c',   mfc='none', ms=4, lw=1, label="$|Z_{rkpq}|$")
iax.plot(aP, c2*aP**1.5,   ':k', mfc='none', ms=4, lw=1, label="{} ".format(c2)+"$b^{3/2}$")
iax.set_yscale('log')
iax.set_xscale('log')


#-- garnish the plot --

#ax[0].hlines(y=0,       xmin=0, xmax=1, linewidth=0.2, color='k')
#ax[1].hlines(y=0,       xmin=0, xmax=1, linewidth=0.2, color='k')
ax[1,0].vlines(x=1,       ymin=1e-6, ymax=1e1, linewidth=0.2, color='k')
ax[0,0].set_xlim(0,1)

ax[1,0].set_ylim(1e-6,1e1)
ax[1,1].set_ylim(1e-6,1e6)

for iax in (ax[0,0], ax[0,1], ax[1,0], ax[1,1]):
    iax.set_xlabel("$a$")
    iax.legend(frameon=True)

ax[0,0].grid()
ax[0,1].grid()

title = "p = [a, -1],  \ \ \ q = [ {}, {} ], \ \ \  r = [0, 1e-5]".format( q[0], q[1])
ax[0,0].set_title(title)

#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.10, right=0.98, top=0.95, bottom=0.05, wspace=0.2,  hspace=0.2)

plt.savefig(outfile) #pad_inches=0

plt.close()

print("Done.")




#===========================================



#=========================================================
