
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle

phase=True


if (1 == 0):
    
    outfile = "tcorr0400_dir"

    base0 = "t200/s5m10_aaQ"
    base1 = "t400/s5m10_aaQ"
    
    jj = [20, 40, 60, 80, 100, 120, 140, 160, 180]

    kk = [1,2,3,4,7]

    def J(i):
        return 6+i

if (1 == 0):
    
    outfile = "tcorr0400_inv"

    base0 = "t200/s5m90_aaQ"
    base1 = "t400/s5m90_aaQ"
    
    jj = [20, 40, 60, 80, 100, 120, 140, 160, 180]

    kk = [1,2,3,4,7]

    def J(i):
        return i


if (1 == 0):
    
    outfile = "tcorr5000_dir"

    base0 = "t5000/s5m10_aa2500_Q"
    base1 = "t5000/s5m10_aa5000_Q"
    
    jj = [100, 120, 140, 160, 180]

    kk = [1,2,3,4,7]

    def J(i):
        return 2+i



if (1 == 1):
    
    outfile = "tcorr5000_inv"

    base0 = "t5000/s5m90_aa2500_Q"
    base1 = "t5000/s5m90_aa5000_Q"
    
    jj = [20, 40, 60, 80, 100]

    kk = [1,2,3,4,7]

    def J(i):
        return i


if phase:
    outfile = outfile + "_phi.pdf"
else:
    outfile = outfile + "_amp.pdf"
    
outfile = "t100.pdf"


#-----------------------------------------------

modes=len(jj)

ntimes=200
dt = 0.1
t=np.arange(ntimes)*dt

fig, ax = plt.subplots(ncols=3, nrows=5, figsize=(15,15))


for i1 in range(5):
    k=kk[i1]

    Q0 = np.fromfile(base0 + str(k) + ".dat")
    Q1 = np.fromfile(base1 + str(k) + ".dat")

    Q0 = Q0.reshape((ntimes, modes, 2)) 
    Q1 = Q1.reshape((ntimes, modes, 2))

    Q0 = Q0[:,:,0] + 1j*Q0[:,:,1]
    Q1 = Q1[:,:,0] + 1j*Q1[:,:,1]

    P0 = np.angle(Q0)/np.pi
    P1 = np.angle(Q1)/np.pi

    Q0 = np.abs(Q0)
    Q1 = np.abs(Q1)


    for i0 in (0,1,2):

        j = J(i0)

        iax=ax[i1][i0]

        if phase:
            iax.plot(t, P0[:,j], ':r',  mfc='none', ms=1, lw=0.5)
            iax.plot(t, P1[:,j], ':ob', mfc='none', ms=1, lw=0.5)
        else:
            iax.plot(t, Q0[:,j], ':r',  mfc='none', ms=1, lw=0.5)
            iax.plot(t, Q1[:,j], ':ob', mfc='none', ms=1, lw=0.5)

        iax.set_title("k={},   j={}".format(k, jj[j]))
        iax.set(xlabel='$\\tau$')
        iax.set_xlim(0,20)

#-- canvas options --


    
#iax = ax[0]  
#iax.set_xlim(0,0.1)
#iax.set_ylim(-5,35)
#iax.set(xlabel='$i$', ylabel='$n=C_0$')
#iax.grid(axis="y")
#iax.legend(frameon=False)


iax = ax[1]  
#iax.set_xlim(0,0.1)
#iax.set_ylim(-5,35)
#iax.set(xlabel='$i$', ylabel='$J$')
#iax.grid(axis="y")
#iax.legend(frameon=False)


#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.04, right=0.99, top=0.97, bottom=0.04, wspace=0.2,  hspace=0.35)

plt.savefig(outfile) #pad_inches=0

plt.close()

print("Done.")




#===========================================



#=========================================================
