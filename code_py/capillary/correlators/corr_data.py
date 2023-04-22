#import importlib ; import plotCap256spectrum
#importlib.reload(plotCap256spectrum) ; plotCap256spectrum.Plot()

import matplotlib
#matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import numpy as np


#------------------------------------------

prefix = "../pyPost/c5";

N = 60
mmax = 6

corrname="+10+06-03+05"
filename = prefix + corrname + ".npy"

#-----------------------------------------

def ind2mode(ind):
    
    i0 = ind[0]
    i1 = ind[1]
    if i0 > N/2:
        i0 = i0 - N
    if i1 > N/2:
        i1 = i1 - N
    return(i1,i0)

c = np.load(filename)

cabs = np.abs(c)
cphi = np.angle(c)/np.pi
cavg = np.average(cabs)


print("")
print("Corr = {},  ({}), C is time average of complex correlator".format(corrname,  filename))
print("Average over modes <|C|> = {:9.3e}".format(cavg))
print("First {} modes with highest |C|:".format(mmax))

print("     1.mode         2.|C|   3.|C|/<|C|>  4.phase/pi" )


m=0
while m<mmax: 
    cmax = np.max(cabs)
    ii   = np.argmax(cabs)
    ind  = np.unravel_index(ii, cabs.shape)
    mode = ind2mode(ind)
    phi  = cphi[ind] 

    print("    [{:+3d} {:+3d}]     {:9.3e}   {:8.2f}   {:6.2f}".format(mode[0], mode[1], cmax, cmax/cavg, phi))
    
    cabs[ind] = 0
    m = m+1
    

print("")

#corrnames=("+02+11", "+06+03-04+08", "+07+11", "+10+06-03+05", "+05+03-03+05+12-03", "+07+00+00-07+05+05")
#corrtitles=("(2,11)", "(6,3) + (-4,8)", "(7,11)",  "(10,6) + (-3,5)", "(5,3)+(-3,5)+(12,-3)", "(7,0) + (0,-7) + (5,5)")




exit()




#------------------------------------------
