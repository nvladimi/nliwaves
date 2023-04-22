
import numpy as np

import matplotlib.pyplot as plt
#from matplotlib import rc

import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

fnameout = "tm_burst.pdf"


#-- prepare data --

gamma= 0.1
P = 0.000002
R = 0.001
V = 1.0
nu1 = 4*R*R/gamma
r_1 = 0.0110
r_2=2.453*gamma/(2*V)
t_min = 101540.0
t_max = 101600.0
N_t = 100


step_t = (t_max - t_min)/N_t

for_plot = np.zeros([3, N_t+1])


for j in range(0,N_t+1):
    t = t_min + j*step_t
    for_plot[0,j] = t
    for_plot[2,j] = gamma/(2.0*V) + np.sqrt((r_2 - gamma/(2*V))**2 + r_1**2/2) * \
        np.tanh(-2*V*(t-t_min) * np.sqrt((r_2-gamma/(2*V))**2+r_1**2/2) + \
                0.5*np.log((np.sqrt((r_2-gamma/(2*V))**2 + r_1**2/2) +  r_2 - gamma/(2*V)) \
                           /(np.sqrt((r_2-gamma/(2*V))**2 + r_1**2/2)-r_2+gamma/(2*V))))
    
    for_plot[1,j] = np.sqrt(r_1**2-2*(for_plot[2,j]-gamma/(2*V))**2+2*(r_2-gamma/(2*V))**2)



a,b,c,d,e = np.loadtxt('burst.txt', skiprows=1, unpack=True)


#-- plot --

fig,ax = plt.subplots(ncols=1, nrows=1, figsize=(3.5, 3.25), tight_layout=True)
fig.tight_layout(pad=0)

iax = ax


iax.plot(a/gamma, (c**2+b**2), '-',  color='#FFC000', lw = 2, label = '$\\rho_1^2$ (num)')
iax.plot(a/gamma, (d**2+e**2), '-',  color='#00C0D0', lw = 2, label = '$\\rho_2^2$ (num)')


iax.plot(for_plot[0], for_plot[1]**2, '--', color='#C00000', lw=1.5, label = '$\\rho_1^2$ (theory)')
iax.plot(for_plot[0], for_plot[2]**2, '--', color='#0000FF', lw=1.5, label = '$\\rho_2^2$ (theory)')

         
iax.plot([101500,101600], [(gamma/(2*V))**2, (gamma/(2*V))**2], 'k--', linewidth=0.7)


iax.set(xlim=(101540., 101580.))
iax.set(ylim=(0, 0.0160))


ticks = (0, 0.0025, 0.005, 0.01, 0.015)
iax.set_yticks(ticks)
iax.set_yticklabels(["0", '$\\frac{\gamma^2}{4|V|^2}$', "0.005", "0.01", "0.015"])


iax.yaxis.labelpad = -2

iax.legend(frameon=False, numpoints=4, labelspacing=0.4, loc='best')


iax.set_ylabel('$\\rho_1^2$,  $\\rho_2^2$')
iax.set_xlabel('$t$')

######################################################3
plt.savefig(fnameout, pad_inches=0)
plt.close()     



#plt.show()    
    
