#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 15:40:46 2020

@author: belan
"""


import math 
#import sympy.mpmath as mp
#from mpmath import *
from numpy import *
import numpy as np
import scipy

#from matplotlib import rc
'''rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)'''
#
'''
from matplotlib import rc
import matplotlib.pylab
import os
#rc('font', **{ 'family':'serif' })
rc('text', usetex=True)
rc('text.latex', unicode = True )
rc('text.latex', preamble = '\usepackage[utf8]{inputenc}')
rc('text.latex', preamble = '\usepackage[russian]{babel}')'''
#
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.pyplot import figure, show, savefig
from matplotlib import cm, colors
import time

font = {'style' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

rc('font', **font)
'''
rc('xtick', labelsize=12) 
rc('ytick', labelsize=12) '''






gamma= 0.1
P = 0.000002
R = 0.001
V = 1.0
nu1 = 4*R*R/gamma
print nu1

r_1 = 0.0110
r_2=2.453*gamma/(2*V)
print gamma/(2*V)
print r_2


t_min = 101540.0
t_max = 101600.0
N_t = 100
step_t = (t_max - t_min)/N_t

for_plot = np.zeros([3, N_t+1])


    

for j in range(0,N_t+1):
    t = t_min + j*step_t
    for_plot[0,j] = t
    for_plot[2,j] = gamma/(2.0*V)+np.sqrt((r_2-gamma/(2*V))**2+r_1**2/2)*np.tanh(-2*V*(t-t_min)*np.sqrt((r_2-gamma/(2*V))**2+r_1**2/2)+0.5*np.log((np.sqrt((r_2-gamma/(2*V))**2+r_1**2/2)+r_2-gamma/(2*V))/(np.sqrt((r_2-gamma/(2*V))**2+r_1**2/2)-r_2+gamma/(2*V))))
    for_plot[1,j] = np.sqrt(r_1**2-2*(for_plot[2,j]-gamma/(2*V))**2+2*(r_2-gamma/(2*V))**2)


print for_plot[2,0]  
#
#label_size = 10
#plt.rcParams['xtick.labelsize'] = label_size 
#plt.rcParams['ytick.labelsize'] = label_size








#plt.rc('font', **font)

#plt.rc('lines', linewidth=20, color='r')


#plt.rc('font', size=12)          # controls default text sizes
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
##plt.rc('figure', titlesize=30)  # fontsize of the figure title
##rc('text', usetex=True)
#
#'''
#gs = gridspec.GridSpec(4, 2,
#                       width_ratios=[1, 1],
#                       height_ratios=[1,2.5,1,2.5]
#                       )
#'''
#
fig = plt.figure(figsize=(8, 10), dpi=80)
##fig = plt.figure(figsize=(7, 14),dpi=80)

######################################################3
#plt.subplot(221)



#ax1 = plt.subplot(gs[2])
#fig, ax1 = plt.subplots()
#t = np.arange(0.01, 10.0, 0.01)
#s1 = np.exp(t)

#fig3 = ax1.plot(for_plot_2[0], for_plot_2[1], '#777777', linestyle='None', marker='o', markersize=8, label = '$p_r$' )
#fig3333=ax1.plot(fp[0], fp[1], '#777777', marker='o', markersize=8, linestyle='None', label = '$p_r$')
#plt.locator_params(axis = 'x', nbins=6)
#plt.locator_params(axis = 'y', nbins=5)
#plt.annotate(r'$p_r$', xy=(2., 0.29), xytext=(2, 0.3))
#plt.annotate(r'optimal', xy=(1., 0.43), xytext=(0.98, 0.423))
#plt.annotate(r'Poisson', xy=(1., 0.43), xytext=(0.98, 0.415))
#plt.annotate(r'restart', xy=(1., 0.43), xytext=(0.98, 0.407))
#plt.annotate('', xy=(0.78,0.402), xytext=(0.92, 0.41),
 #           arrowprops=dict(facecolor='black', frac=0.4, shrink=0.01),
 #           )
# Make the y-axis label, ticks and tick labels match the line color.

#ax1.set_xlabel('$r$')

#fp = np.genfromtxt("burst1.txt")

#f = open('burst1.txt', 'r')






a,b,c,d,e = np.loadtxt('burst.txt', skiprows=1, unpack=True)

#with open('burst1.txt', 'r') as f:
#    lines = f.readlines()
#    x = [float(line.split()[0]) for line in lines]
#    y = [float(line.split()[1]) for line in lines]



#x = [row.split(' ')[0] for row in dat]
#y = [row.split(' ')[1] for row in dat]


fig = plt.figure(figsize=(6, 5), dpi=120)

#print rho[1]

#plt.plot(a, (c**2+b**2)/(2*P**2), '#348abd', marker='o', markersize=1, linestyle='None', label = '$\\rho_1$')
#plt.plot(a, 10*(d**2+e**2)/(2*P**2), '#8eba42', marker='o', markersize=1, linestyle='None', label = '$\\rho_2$')

plt.plot(for_plot[0], for_plot[1]**2, '#e24a33', linestyle='--', label = '$\\rho_1^2$ (theory)')
plt.plot(for_plot[0], for_plot[2]**2, '#988ed5', linestyle='--', label = '$\\rho_2^2$ (theory)')




plt.plot(a/gamma, (c**2+b**2), '#fbc15e', marker='.', markersize=1, linestyle='None', label = '$\\rho_1^2$ (numerics)')
plt.plot(a/gamma, (d**2+e**2), '#348abd', marker='.', markersize=1, linestyle='None', label = '$\\rho_2^2$ (numerics)')

         
plt.plot([101500,101600], [(gamma/(2*V))**2, (gamma/(2*V))**2], 'k--', linewidth=1., label = '$\\frac{\gamma^2}{4|V|^2}$')


plt.xlim(101540., 101580.)
#plt.ylim(0., 1.1*r_2)
#s2 = np.sin(2 * np.pi * t)

'''
lns = fig22+fig2+fig3333  #+fig3+fig222+fig0+fig1
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, labelspacing=0.08, prop={'size':10}, loc=4)
'''

#plt.locator_params(axis = 'y', nbins=7)
#plt.locator_params(axis = 'x', nbins=4)
#ax2.tick_params('y', colors='r')
plt.legend(numpoints=4, labelspacing=0.1, loc='best')




plt.ylabel('amplitude, $\\rho_i^2$')
plt.xlabel('time, $t$')

######################################################3
plt.show()    
    