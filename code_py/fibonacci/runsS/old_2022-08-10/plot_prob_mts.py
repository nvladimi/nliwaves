
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
from scipy.special import gamma, factorial 

matplotlib.rc('xtick', labelsize=11) 
matplotlib.rc('ytick', labelsize=11) 
matplotlib.rc('axes', labelsize=13)



    
outfile = "plot_prob_mts.pdf"

direct=True

fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,2.9))

LW = 1.2
mycolors=['#d00000',  '#40a000', '#0040ff', '#e07000', '#00a0f0',  '#c000f0']


#---------------------------------

iax=ax[0]

if direct:
   run="MIx32x1/s5m10"; irun=1
else:
   run="MIx32x1/s5m90"; irun=0

modes= ((("20", "$\,20$", "-",   3), 
         ("40", "$\,40$", "-",   1), 
         ("80", "$\,80$", "--",  2), 
#        ("120","$120$","--",    3), 
         ("160","$160$", "--",   1), 
#        ("180","$180$", ":",    5)
         ), 
        (
         ("180", " $20$", "-",   3), 
         ("160", " $40$", "-",   1), 
         ("120", " $80$", "--",  2),
#        ("80", "$120$", "--",   ), 
         ("40", "$160$", "--",    0), 
#        ("20", "$180$", ":",    5)
        ) 
)
 


for r in modes[irun]:

   fname = run + "_i" + r[0] + "_prob.txt"

   label = r[1]
   style = r[2]
   color = mycolors[r[3]]
   dat = np.loadtxt(fname)

   x = dat[:,0]
   y = np.log(dat[:,3])
   iax.plot(x,y, style, mfc='none', lw=LW, color=color,  label=label);
        
x=np.arange(0,21)
iax.plot(x,-x, '-k', lw=LW/3)

   
iax.set_xlim(0,18)
iax.set_ylim(-18,0)
iax.set(xlabel='$|b_i|^2 / \langle |b_i|^2 \\rangle $')
iax.set_yticks(np.arange(-18, 1, step=6)) 
iax.yaxis.set_label_coords(-0.10,0.5)

iax.text(1.5, -11, "$|i-d|$", fontsize=12)
iax.legend(frameon=False, loc="lower left", handlelength=3)
#bbox_to_anchor=(1.0, 0.94), labelspacing = 0.4)

iax.set(ylabel='$\ln {\cal P}$')

  
   
#---------------------------
iax=ax[1]

#MS=0.5

#nrange = (1,2,3,6,9,12)
nrange = (1,4,8,12)

m=200

if direct:
   datM  = np.loadtxt("Post/s5m10_mts.txt"); d=m
else:
   datM  = np.loadtxt("Post/s5m90_mts.txt"); d=1


x = np.abs(d-np.arange(1,m+1))

ic=3
iax.set_prop_cycle(None)
for n in nrange:
   y = datM[:,n+1] / gamma(n/2+1) / (datM[:,2+1])**(n/2)
   iax.plot(x,y, mfc='none',  lw=LW,  label="{}".format(n), color=mycolors[ic]);
   ic = ic-1
   
iax.set_xlim(0,200)
iax.set_ylim(0.75,3)
iax.set(xlabel='$|i-d|$')

iax.set(ylabel='$A_m  A_2^{-m/2} / \Gamma(\\frac{m}{2}+1)$')
iax.yaxis.set_label_coords(-0.13,0.5)


#---------------------------

inax = ax[1].inset_axes([0.37, 0.37, 0.58, 0.58])

inax.set_prop_cycle(None)

ic=3
for n in nrange:
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   inax.plot(x,y, mfc='none',  lw=LW, label="${}$".format(n), color=mycolors[ic]);
   ic = ic - 1
   
inax.plot(x,1.13*x, ":k") # label="$1.13|i-d|$")
inax.axhline( y=28, xmin=0.29, xmax=0.42, ls=":", color="k")


inax.set_xlim(0,200)
inax.set_ylim(0,250)
inax.set(xlabel='$|i-d|$')
inax.legend(frameon=False, handlelength=1)

inax.legend(frameon=False, labelspacing = 0.32, handlelength=1,
            loc="upper left", bbox_to_anchor=(0.02, 0.93))
inax.text(30, 225, "$m$", fontsize=11)



inax.text(95, 20, "$1.13|i-d|$", fontsize=11)

inax.set(ylabel='$[ \\frac{A_m}{\Gamma(m/2+1)} ]^{\\frac{3}{m}}$')
inax.yaxis.set_label_coords(-0.25,0.5)
inax.yaxis.label.set_fontsize(15)
inax.xaxis.label.set_fontsize(11)




#-----------------------

plt.subplots_adjust(left=0.07, right=0.98, top=0.97, bottom=0.16, wspace=0.27)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
