
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
from scipy.special import gamma, factorial 

matplotlib.rc('xtick', labelsize=11) 
matplotlib.rc('ytick', labelsize=11) 
matplotlib.rc('axes', labelsize=13)


#corr="21-0"
#corr="11-30"
#corr="322-0"

#corr="222-40"

#corr="4111-00"
#corr="44333-0"
corr="63333-0"

fig, ax = plt.subplots(ncols=6, nrows=2, figsize=(18,6))

LW = 1.2
MS = 6


#--data --

m=200; dd=m; di=1;

j=np.array((20, 40, 60, 80, 100, 120, 140, 160, 180))
ji = np.abs(j-di)
jd = np.abs(j-dd)

j0 =np.array((0, 200))


Dd  = np.loadtxt("Post/maxcor_dir" + corr + "_b100.txt")
Di  = np.loadtxt("Post/maxcor_inv" + corr + "_b100.txt")
outfile = "plot_maxcorr_" + corr + ".pdf"

if corr == "21-0":
   lbl1="$t_1$";  lbl2="$t_2$"
   
elif corr == "11-30":
   lbl1="$t_1$";  lbl2="$t_3$"

elif corr == "322-0":
   lbl1="$t_2$";  lbl2="$t_3$"

elif corr == "222-40":
   lbl1="$t_2$";  lbl2="$t_4$"

elif corr == "4111-00":
   lbl1="$t_1$";  lbl2="$t_4$"

elif corr == "44333-0":
   lbl1="$t_3$";  lbl2="$t_4$"

elif corr == "63333-0":
   lbl1="$t_3$";  lbl2="$t_6$"
   
else:
   print("Unknown correlator")


#0.j  1.t1[0]  2.t1[1]  3.Q1  4.t3[0] 5.t3[1]  6.Q3  7.c00  8.c00m
# 9.d1q1  10.d1q2  11.ddq1  12.d1q2  13.d3q1  14.d3q2  15.d4q1  16.d4q2

T=0.1

#-- direct cascade --

iax=ax[0,0]
iax.plot(jd, Dd[:,3], '-^r', mfc='none', ms=MS, lw=LW/2,  label=lbl1+','+lbl2+"\,$>0$");
iax.plot(jd, Dd[:,6], '-vb', mfc='none', ms=MS, lw=LW/2,  label=lbl1+','+lbl2+"\,$<0$");

f=0.5*(np.abs(Dd[:,3]) + np.abs(Dd[:,6]))

iax=ax[0,1]
iax.plot(jd, Dd[:,1], '-or',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(jd, Dd[:,2], '--sr', mfc='none', ms=MS, lw=LW/2,  label=lbl2);
iax.plot(jd, Dd[:,4], '-ob',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(jd, Dd[:,5], '--sb', mfc='none', ms=MS, lw=LW/2,  label=lbl2);

iax=ax[0,2]
iax.plot(jd, Dd[:,9]/f*T,  '-om',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(jd, Dd[:,10]/f*T, '--sm', mfc='none', ms=MS, lw=LW/2,  label=lbl2);

iax=ax[0,3]
iax.plot(jd, Dd[:,11]/f*T**2, '-om',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(jd, Dd[:,12]/f*T**2, '--sm', mfc='none', ms=MS, lw=LW/2,  label=lbl2);

iax=ax[0,4]
iax.plot(jd, Dd[:,13]/f*T**3, '-om',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(jd, Dd[:,14]/f*T**3, '--sm', mfc='none', ms=MS, lw=LW/2,  label=lbl2);

iax=ax[0,5]
iax.plot(jd, Dd[:,15]/f*T**4, '-om',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(jd, Dd[:,16]/f*T**4, '--sm', mfc='none', ms=MS, lw=LW/2,  label=lbl2);




#-- inverse cascade --

iax=ax[1,0]
iax.plot(ji, Di[:,3], '-^r', mfc='none', ms=MS, lw=LW/2,  label=lbl1+','+lbl2+"$\,>0$");
iax.plot(ji, Di[:,6], '-vb', mfc='none', ms=MS, lw=LW/2,  label=lbl1+','+lbl2+"$\,<0$");

f=0.5*(np.abs(Di[:,3]) + np.abs(Di[:,6]))

iax=ax[1,1]
iax.plot(ji, Di[:,1], '-or',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(ji, Di[:,2], '--sr', mfc='none', ms=MS, lw=LW/2,  label=lbl2);
iax.plot(ji, Di[:,4], '-ob',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(ji, Di[:,5], '--sb', mfc='none', ms=MS, lw=LW/2,  label=lbl2);

iax=ax[1,2]
iax.plot(ji, Di[:,9]/f*T,  '-om',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(ji, Di[:,10]/f*T, '--sm', mfc='none', ms=MS, lw=LW/2,  label=lbl2);

iax=ax[1,3]
iax.plot(ji, Di[:,11]/f*T**2, '-om',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(ji, Di[:,12]/f*T**2, '--sm', mfc='none', ms=MS, lw=LW/2,  label=lbl2);

iax=ax[1,4]
iax.plot(ji, Di[:,13]/f*T**3, '-om',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(ji, Di[:,14]/f*T**3, '--sm', mfc='none', ms=MS, lw=LW/2,  label=lbl2);

iax=ax[1,5]
iax.plot(ji, Di[:,15]/f*T**4, '-om',  mfc='none', ms=MS, lw=LW/2,  label=lbl1);
iax.plot(ji, Di[:,16]/f*T**4, '--sm', mfc='none', ms=MS, lw=LW/2,  label=lbl2);



#-- garnish plots --

ax[0,0].set_title("value at extrema, $F=|f(ext)|$")
ax[0,1].set_title("times at extrema, typical $T=0.1$")
ax[0,2].set_title("$f'(0) \, T/F$")
ax[0,3].set_title("$f''(0) \, T^2/F$")
ax[0,4].set_title("$f'''(0) \, T^3/F$")
ax[0,5].set_title("$f''''(0) \, T^4/F$")

for iax in (ax[1,0], ax[1,1], ax[1,2], ax[1,3], ax[1,4], ax[1,5]):
   iax.set(xlabel='$|j-d|$')

for iax in (ax[0,0], ax[0,1], ax[0,2], ax[0,3], ax[0,4], ax[0,5],
            ax[1,0], ax[1,1], ax[1,2], ax[1,3], ax[1,4], ax[1,5]):
   iax.set_xlim(0,200)

for iax in (ax[0,1], ax[1,1]):
   iax.grid("on")
   
   
if corr == "21-0":
   ax[0,0].set_ylim(-200,200)
   ax[1,0].set_ylim(-200,200)
   ax[0,0].plot(j0,  1.1*j0, '--k', lw=LW/2,  label="$1.1 \, |j-d|$");
   ax[0,0].plot(j0, -1.1*j0, '--k', lw=LW/2,  label="$-1.1 \,|j-d|$");
   ax[1,0].plot(j0,  1.0*j0, '--k', lw=LW/2,  label="$1.0 \, j-d|$");
   ax[1,0].plot(j0, -1.0*j0, '--k', lw=LW/2,  label="$-1.0 \, |j-d|$");
   ax[0,0].text(15, 160, '"'+corr+'"')
   ax[0,0].text(20, 130, 'direct')
   ax[1,0].text(15, 160, '"'+corr+'"')
   ax[1,0].text(20, 130, 'inverse')
  
elif corr == "11-30":
   ax[0,0].set_ylim(-120,0)
   ax[1,0].set_ylim(0,120)
   ax[0,0].plot(j0, -0.7*j0, '--k', lw=LW/2,  label="$-0.7\,|j-d|$");
   ax[1,0].plot(j0,  0.7*j0, '--k', lw=LW/2,  label="$0.7\,|j-d|$");
   ax[0,0].text(15, -100, '"'+corr+'"')
   ax[0,0].text(20, -110, 'direct')
   ax[1,0].text(145, 20, '"'+corr+'"')
   ax[1,0].text(150, 10, 'inverse')

elif corr == "322-0":
   ax[0,0].set_ylim(-250,0)
   ax[1,0].set_ylim(0,250)
   ax[0,0].plot(j0, -1.4*j0, '--k', lw=LW/2,  label="$-1.4\,|j-d|$");
   ax[1,0].plot(j0,  1.4*j0, '--k', lw=LW/2,  label="$1.4\,|j-d|$");
   ax[0,0].text(15, -210, '"'+corr+'"')
   ax[0,0].text(20, -230, 'direct')
   ax[1,0].text(145, 40, '"'+corr+'"')
   ax[1,0].text(150, 20, 'inverse')
   
elif corr == "222-40":
   ax[0,0].set_ylim(-750,750)
   ax[1,0].set_ylim(-750,750)
   ax[0,0].plot(j0,  3.8*j0, '--k', lw=LW/2,  label="$3.8 \, |j-d|$");
   ax[0,0].plot(j0, -3.8*j0, '--k', lw=LW/2,  label="$-3.8 \,|j-d|$");
   ax[1,0].plot(j0,  3.7*j0, '--k', lw=LW/2,  label="$3.7 \, j-d|$");
   ax[1,0].plot(j0, -3.7*j0, '--k', lw=LW/2,  label="$-3.7 \, |j-d|$");
   ax[0,0].text(15, 620, '"'+corr+'"')
   ax[0,0].text(20, 500, 'direct')
   ax[1,0].text(15, 620, '"'+corr+'"')
   ax[1,0].text(20, 500, 'inverse')

elif corr == "4111-00":
   ax[0,0].set_ylim(-175,0)
   ax[1,0].set_ylim(0,175)
   ax[0,0].plot(j0, -0.9*j0, '--k', lw=LW/2,  label="$-0.9\,|j-d|$");
   ax[1,0].plot(j0,  0.9*j0, '--k', lw=LW/2,  label="$0.9\,|j-d|$");
   ax[0,0].text(15, -150, '"'+corr+'"')
   ax[0,0].text(20, -165, 'direct')
   ax[1,0].text(145, 25, '"'+corr+'"')
   ax[1,0].text(150, 10, 'inverse')

elif corr == "44333-0":
   ax[0,0].set_ylim(-500,0)
   ax[1,0].set_ylim(0,500)
   ax[0,0].plot(j0, -2.8*j0, '--k', lw=LW/2,  label="$-2.8\,|j-d|$");
   ax[1,0].plot(j0,  2.6*j0, '--k', lw=LW/2,  label="$2.6\,|j-d|$");
   ax[0,0].text(15, -440, '"'+corr+'"')
   ax[0,0].text(20, -480, 'direct')
   ax[1,0].text(145, 80, '"'+corr+'"')
   ax[1,0].text(150, 40, 'inverse')


elif corr == "63333-0":
   ax[0,0].set_ylim(-400,400)
   ax[1,0].set_ylim(-400,400)
   ax[0,0].plot(j0, -2.0*j0, '--k', lw=LW/2,  label="$-2.0\,|j-d|$");
   ax[1,0].plot(j0,  2.0*j0, '--k', lw=LW/2,  label="$2.0\,|j-d|$");
   ax[0,0].text(15, -300, '"'+corr+'"')
   ax[0,0].text(20, -350, 'direct')
   ax[1,0].text(15, 340, '"'+corr+'"')
   ax[1,0].text(20, 290, 'inverse')

   
else:
   print("Unknown correlator " + corr )


for iax in (ax[0,0], ax[0,2], ax[0,3], ax[0,4], ax[0,5],
            ax[1,0], ax[1,2], ax[1,3], ax[1,4], ax[1,5]):
   iax.legend(frameon=False)
 

#---------------------------

#-----------------------

plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.10, wspace=0.30, hspace=0.25)

plt.savefig(outfile) #pad_inches=0

plt.close()


print("Done.")


#===========================================



#=========================================================
