
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
outfile = "tcorr_auto.pdf"


jj = [20, 40, 60, 80, 100, 120, 140, 160, 180]
#jj = [20, 40, 60, 80, 100]


#modes0=len(jj)
modes=len(jj)

ntimes=200
dt = 0.1

Q0 = np.fromfile("t200/s5m10_aaQ0.dat")
R0 = np.fromfile("t200/s5m10_aaR0.dat")
Q1 = np.fromfile("t200/s5m90_aaQ0.dat")
R1 = np.fromfile("t200/s5m90_aaR0.dat")


Q0 = Q0.reshape((ntimes, 9, 2)) 
Q1 = Q1.reshape((ntimes, 9, 2)) 
R0 = R0.reshape((ntimes, 9))
R1 = R1.reshape((ntimes, 9))

Q0 = Q0[:,:,0] + 1j*Q0[:,:,1]
Q1 = Q1[:,:,0] + 1j*Q1[:,:,1]

Q0 = np.abs(Q0)
Q1 = np.abs(Q1)



#print(Q0.shape, Q1.shape, R0.shape, R1.shape)


fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(6,9))

t=np.arange(ntimes)*dt

iax=ax[0]


j=0

iax=ax[0][0]
iax.plot(t, Q0[:,j], '-or', mfc='none', ms=1, lw=0.5,
         label='$|\langle a_j(t) \, a_j^*(t+\\tau)\\rangle |$')
iax.plot(t, R0[:,j], '-ob', mfc='none', ms=1, lw=0.5,
         label='$\langle |a_j(t)| \, |a_j^*(t+\\tau)|\\rangle$')
iax.set_title("direct,   j={}".format(jj[j]))

iax=ax[0][1]
iax.plot(t, Q1[:,j], '-or', mfc='none', ms=1, lw=0.5,
         label='$|\langle a_j(t) \, a_j^*(t+\\tau)\\rangle |$')
iax.plot(t, R1[:,j], '-ob', mfc='none', ms=1, lw=0.5,
         label='$\langle |a_j(t)| \, |a_j^*(t+\\tau)|\\rangle$')
iax.set_title("inverse,   j={}".format(jj[j]))



j=4

iax=ax[1][0]
iax.plot(t, Q0[:,j], '-or', mfc='none', ms=1, lw=0.5,
         label='$|\langle a_j(t) \, a_j^*(t+\\tau)\\rangle |$')
iax.plot(t, R0[:,j], '-ob', mfc='none', ms=1, lw=0.5,
         label='$\langle |a_j(t)| \, |a_j^*(t+\\tau)|\\rangle$')
iax.set_title("direct,   j={}".format(jj[j]))

iax=ax[1][1]
iax.plot(t, Q1[:,j], '-or', mfc='none', ms=1, lw=0.5,
         label='$|\langle a_j(t) \, a_j^*(t+\\tau)\\rangle |$')
iax.plot(t, R1[:,j], '-ob', mfc='none', ms=1, lw=0.5,
         label='$\langle |a_j(t)| \, |a_j^*(t+\\tau)|\\rangle$')
iax.set_title("inverse,   j={}".format(jj[j]))


j=8

iax=ax[2][0]
iax.plot(t, Q0[:,j], '-or', mfc='none', ms=1, lw=0.5,
         label='$|\langle a_j(t) \, a_j^*(t+\\tau)\\rangle |$')
iax.plot(t, R0[:,j], '-ob', mfc='none', ms=1, lw=0.5,
         label='$\langle |a_j(t)| \, |a_j^*(t+\\tau)|\\rangle$')
iax.set_title("direct,   j={}".format(jj[j]))

iax=ax[2][1]
iax.plot(t, Q1[:,j], '-or', mfc='none', ms=1, lw=0.5,
         label='$|\langle a_j(t) \, a_j^*(t+\\tau)\\rangle |$')
iax.plot(t, R1[:,j], '-ob', mfc='none', ms=1, lw=0.5,
         label='$\langle |a_j(t)| \, |a_j^*(t+\\tau)|\\rangle$')
iax.set_title("inverse,   j={}".format(jj[j]))




#-- canvas options --

for i in (0,1,2):
    for j in (0,1):
        ax[i][j].legend(frameon=False)
        ax[i][j].set_xlim(0,1)
        ax[i][j].set(xlabel='$\\tau$')

    
#iax = ax[0]  
#iax.set_xlim(0,0.1)
#iax.set_ylim(-5,35)
#iax.set(xlabel='$i$', ylabel='$n=C_0$')
#iax.grid(axis="y")
#iax.legend(frameon=False)


#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.06, right=0.98, top=0.96, bottom=0.06, wspace=0.3,  hspace=0.35)

plt.savefig(outfile) #pad_inches=0

plt.close()

print("Done.")




#===========================================



#=========================================================
