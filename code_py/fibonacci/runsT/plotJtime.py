
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
import fiboPost

outfile = "plotJtime.pdf"

nt1=100
nt2=100
dt = 0.01

#-- s = (1j)**(len(jRe) - len(jIm))

corr = "21-0" ;     s=(1j)**(2 - 1);  lbl1="$t_1$"; lbl2="$t_2$"; 
#corr = "322-0" ;    s=(1j)**(3 - 1);  lbl1="$t_2$"; lbl2="$t_3$"; 
#corr = "11-30" ;    s=(1j)**(2 - 2);  lbl1="$t_1$"; lbl2="$t_3$"; 
#corr = "222-40" ;   s=(1j)**(3 - 2);  lbl1="$t_2$"; lbl2="$t_4$"; 
#corr = "4111-00" ;  s=(1j)**(4 - 2);  lbl1="$t_1$"; lbl2="$t_4$"; 
#corr = "44333-0" ;  s=(1j)**(5 - 1);  lbl1="$t_3$"; lbl2="$t_4$"; 
#corr = "63333-0" ;  s=(1j)**(5 - 1);  lbl1="$t_3$"; lbl2="$t_6$"; 




fbase = "t1inv"; d=1
#fbase = "t1dir"; d=200

for j in (20, 40, 60, 80, 100, 120, 140, 160, 180):
#for j in (100,):
    
    Q = np.fromfile("Post/b100/" + fbase + "_" + corr + "_j" + str(j).zfill(3) +".dat")

    Q = Q.reshape((2*nt1, 2*nt2, 2))
    Q = np.real(  (Q[:,:,0] + (1j)*Q[:,:,1])*s )
    
    t1=np.arange(-nt1,nt1)*dt
    t2=np.arange(-nt2,nt2)*dt

    c00  = Q[nt1,nt2]
    c00m = fiboPost.Cmodel(corr, j, d)

    Q = Q/c00m

    #-- minumum and maximum --
    
    imax = np.unravel_index(np.argmax(Q, axis=None), Q.shape)
    tmax=( np.array(imax[0:2]) - np.array((nt1,nt2)) )*dt
    qmax=Q[imax]

    imin = np.unravel_index(np.argmin(Q, axis=None), Q.shape)
    tmin=(np.array(imin[0:2]) - np.array((nt1,nt2)) )*dt
    qmin=Q[imin]

    #-- extrema in 1st and 3rd quaters --
    
    Q1  =  Q[nt1:, nt2:]
    ie1 = np.unravel_index(np.argmax(np.abs(Q1), axis=None), Q1.shape)
    te1 = np.array(ie1[0:2])*dt
    Qe1  = Q1[ie1]

    Q3  = Q[0:nt1+1, 0:nt2+1]
    ie3 = np.unravel_index(np.argmax(np.abs(Q3), axis=None), Q3.shape)
    te3 = (np.array(ie3[0:2]) - np.array((nt1,nt2)))*dt
    Qe3 = Q3[ie3]

    #-- dervatives at the center --

    q1   = Q[:,nt2]
    d1q1 = ( q1[nt1-2]  -  8*q1[nt1-1]              +  8*q1[nt1+1]  - q1[nt1+2]) /12/dt
    d2q1 = (-q1[nt1-2]  + 16*q1[nt1-1] - 30*q1[nt1] + 16*q1[nt1+1]  - q1[nt1+2]) /12/dt/dt
    d3q1 = (-q1[nt1-2]  +  2*q1[nt1-1]              -  2*q1[nt1+1]  + q1[nt1+2]) /2/dt/dt/dt
    d4q1 = ( q1[nt1-2]  -  4*q1[nt1-1] +  6*q1[nt1] -  4*q1[nt1+1]  + q1[nt1+2]) /dt/dt/dt/dt
    
    q2   = Q[nt1,:]
    d1q2 = ( q2[nt2-2]  -  8*q2[nt2-1]              +  8*q2[nt2+1]  - q2[nt2+2]) /12/dt
    d2q2 = (-q2[nt2-2]  + 16*q2[nt2-1] - 30*q2[nt2] + 16*q2[nt2+1]  - q2[nt2+2]) /12/dt/dt
    d3q2 = (-q2[nt2-2]  +  2*q2[nt2-1]              -  2*q2[nt2+1]  + q2[nt2+2]) /2/dt/dt/dt
    d4q2 = ( q2[nt2-2]  -  4*q2[nt2-1] +  6*q2[nt2] -  4*q2[nt2+1]  + q2[nt2+2]) /dt/dt/dt/dt


# 3rd   −1/2   1   0   −1   1/2
# 3rd   ( -1   2   0   -2   1 )*(1/2)
# 4th      1  −4   6   −4   1

    
       
    outstr1="{:4.0f}   {:5.2f}  {:5.2f}  {:7.2f}   {:5.2f}  {:5.2f}  {:7.2f}  {:7.4f}  {:7.4f}"
    outstr2="    {:9.2e} {:9.2e}   {:9.2e} {:9.2e}   {:9.2e} {:.2e}   {:9.2e} {:9.2e}"

    
    print(outstr1.format(j, te1[0], te1[1], Qe1, te3[0], te3[1], Qe3,  c00,  c00m), \
          outstr2.format(d1q1, d1q2,  d2q1, d2q2, d3q1, d3q2, d4q1, d4q2))
    
    #print(outstr1.format(j, tmax[0], tmax[1], qmax, tmin[0], tmin[1], qmin,  c00,  c00m))

quit()

#-------------------------------

fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(9,3))
ymax=20

#print(imax, imin, tmax, tmin, Q[imax], Q[imin])

iax=ax[0]
#iax.plot(t1, Q[:,nt2,0], '-r', mfc='none', ms=1, lw=0.5)
iax.plot(t1, Q[:,nt2],     '-b', mfc='none', ms=1, lw=1, label=lbl2 + "=0")
iax.plot(t1, Q[:,imax[1]], '-r', mfc='none', ms=1, lw=1, label=lbl2 + "@ max")
iax.plot(t1, Q[:,imin[1]], '-g', mfc='none', ms=1, lw=1, label=lbl2 + "@ min")

iax.hlines(y=0,       xmin=-ymax, xmax=ymax, linewidth=0.2, color='k')
iax.vlines(x=0,       ymin=-ymax, ymax=ymax, linewidth=0.2, color='k')
iax.vlines(x=tmin[0], ymin=qmin,  ymax=0,    linewidth=0.3, color='g')
iax.vlines(x=tmax[0], ymin=0,     ymax=qmax, linewidth=0.3, color='r')

iax.set_xlim(-1,1)
#iax.set_ylim(-ymax,ymax)
iax.legend(frameon=False)



iax=ax[1]
#iax.plot(t2, Q[nt1,:,0], '-r', mfc='none', ms=1, lw=0.5)
iax.plot(t2, Q[nt1,:],     '-b', mfc='none', ms=1, lw=1, label=lbl1 + "=0")
iax.plot(t2, Q[imax[0],:], '-r', mfc='none', ms=1, lw=1, label=lbl1 + "@ max")
iax.plot(t2, Q[imin[0],:], '-g', mfc='none', ms=1, lw=1, label=lbl1 + "@ min")

iax.hlines(y=0,       xmin=-ymax, xmax=ymax, linewidth=0.2, color='k')
iax.vlines(x=0,       ymin=-ymax, ymax=ymax, linewidth=0.2, color='k')
iax.vlines(x=tmin[1], ymin=qmin,  ymax=0,    linewidth=0.3, color='g')
iax.vlines(x=tmax[1], ymin=0,     ymax=qmax, linewidth=0.3, color='r')

iax.set_xlim(-1,1)
#iax.set_ylim(-ymax,ymax)
iax.legend(frameon=False)



iax=ax[2]
extent=((-nt1-0.5)*dt, (nt1-0.5)*dt, (-nt2-0.5)*dt, (nt2-0.5)*dt)

vrange= max(qmax, -qmin)*1.2
im0=iax.imshow(np.transpose(Q[:,:]), cmap="jet", origin = 'lower', vmin=-vrange, vmax=vrange, extent=extent);

iax.hlines(y=0, xmin=-nt1*dt, xmax=nt1*dt, linewidth=0.5, color='k')
iax.vlines(x=0, ymin=-nt2*dt, ymax=nt2*dt, linewidth=0.5, color='k')
iax.set_xlim(-(nt1-0.5)*dt,(nt1-0.5)*dt)
iax.set_ylim(-(nt2-0.5)*dt,(nt2-0.5)*dt)

ax[0].set(xlabel=lbl1)
ax[1].set(xlabel=lbl2)
ax[2].set(xlabel=lbl1)
ax[2].set(ylabel=lbl2)





if corr == "21-0":
    title = "$\langle a_{j-2}(t-t_2) \, a_{j-1}(t-t_1) \, a^*_j(t) \\rangle / C_{21\\bar{0}}$"
    ax[0].text(0.1, -100, "inverse, $j=100$")

if corr == "322-0":
    title = "$\langle a_{j-3}(t-t_3) \, a_{j-2}^2(t-t_2) \, a^*_j(t) \\rangle / C_{322\\bar{0}}$"
    ax[0].text(-0.9, -20, "inverse, $j=100$")

if corr == "11-30":
    title="$\langle a^*_{j-3}(t-t_3) \, a_{j-1}^2(t-t_1) \, a^*_j(t) \\rangle / C_{\\bar{3}11\\bar{0}}$"
    ax[0].text(-0.9, -20, "inverse, $j=100$")

if corr == "222-40":
    title="$\langle a^*_{j-4}(t-t_4) \, a_{j-2}^3(t-t_2) \, a^*_j(t) \\rangle / C_{\\bar{4}222\\bar{0}}$"
    ax[0].text(-0.9, -350, "inverse, $j=100$")

if corr == "4111-00":
    title="$\langle a_{j-4}(t-t_4) \, a_{j-1}^3(t-t_1) \, a^{*2}_j(t) \\rangle / C_{4111\\bar{0}\\bar{0}}$"
    ax[0].text(0.1, -55, "inverse, $j=100$")

if corr == "44333-0":
    title="$\langle a^{2}_{j-4}(t-t_4) \, a_{j-3}^3(t-t_3) \, a^*_j(t) \\rangle / C_{44333\\bar{0}}$"
    ax[0].text(-0.9, -50, "inverse, $j=100$")

if corr == "63333-0":
    title="$\langle a_{j-6}(t-t_6) \, a_{j-3}^4(t-t_3) \, a^*_j(t) \\rangle / C_{63333\\bar{0}}$"
    ax[0].text(0.1, -200, "inverse, $j=100$")

ax[2].set_title(title)

    
#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.04, right=0.97, top=0.91, bottom=0.13, wspace=0.3,  hspace=0.35)

plt.savefig(outfile) #pad_inches=0

plt.close()

print("Done.")




#===========================================



#=========================================================
