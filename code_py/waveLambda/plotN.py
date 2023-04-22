# this is a comment

import kinwave as kw
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

#-------------------------------------------

kw.write_N2D_pi2(2, 2)

#--------------------------------------------

def plot_N2D():

    fname = "sumN_pi2_q02r02.npy"; outfile = "plotN.pdf"

    
    #fname = "2023-04-08/sumN_pi2_q02r02.npy"; outfile = "2023-04-08/sumN_pi2_q02r02.pdf"
    #fname = "2023-04-08/sumN_pi2_q02r04.npy"; outfile = "2023-04-08/sumN_pi2_q02r04.pdf"
    #fname = "2023-04-08/sumN_pi2_q04r02.npy"; outfile = "2023-04-08/sumN_pi2_q04r02.pdf"
    #fname = "2023-04-08/sumN_pi2_q03r03.npy"; outfile = "2023-04-08/sumN_pi2_q03r03.pdf"
    #fname = "2023-04-08/sumN_pi2_q04r04.npy"; outfile = "2023-04-08/sumN_pi2_q04r04.pdf"

    #fname = "2023-04-08/resN_pi2_n02m04.npy"; outfile = "2023-04-08/rezN_pi2_n02m04.pdf"
    #fname = "2023-04-08/resN_pi2_n03m06.npy"; outfile = "2023-04-08/rezN_pi2_n03m06.pdf"
    #fname = "2023-04-08/resN_pi2_n04m08.npy"; outfile = "2023-04-08/rezN_pi2_n04m08.pdf"
    #fname = "2023-04-08/resN_pi2_n08m16.npy"; outfile = "2023-04-08/rezN_pi2_n08m16.pdf"
    #fname = "2023-04-08/resN_pi2_n16m32.npy"; outfile = "2023-04-08/rezN_pi2_n16m32.pdf"
    #fname = "2023-04-08/resN_pi2_n32m64.npy"; outfile = "2023-04-08/rezN_pi2_n32m64.pdf"

    #fname = "2023-04-08/resN_pi2_n02m06.npy"; outfile = "2023-04-08/rezN_pi2_n02m06.pdf"
    #fname = "2023-04-08/resN_pi2_n04m12.npy"; outfile = "2023-04-08/rezN_pi2_n04m12.pdf"
    #fname = "2023-04-08/resN_pi2_n08m24.npy"; outfile = "2023-04-08/rezN_pi2_n08m24.pdf"
    #fname = "2023-04-08/resN_pi2_n16m48.npy"; outfile = "2023-04-08/rezN_pi2_n16m48.pdf"

    
    
    
    with open(fname, 'rb') as f:
        p   = np.load(f)
        p1  = np.load(f)
        (mq, dq, mr, dr) = np.load(f) 
        N2D = np.load(f)
        m = int ( (N2D.shape[0]-1)/2 )
    
            
    print(m,p,p1, np.shape(N2D))
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(8,4)) 

    vrange = np.max( [np.abs(np.max(N2D)), np.abs(np.min(N2D))] ) /10
    print(vrange)
    vrange = 10000
    
    iax = ax[0]
    extent=(-m-0.5,m+0.5,-m-0.5,m+0.5)
    
    im0=iax.imshow(np.transpose(N2D), cmap="jet", origin = 'lower', vmin=-vrange, vmax=vrange, extent=extent)


    plt.subplots_adjust(left=0.05, right=0.95, top=0.99, bottom=0.05)


    plt.savefig(outfile) #pad_inches=0

    plt.close()


    return() 

plot_N2D()

print("Done.")




#===========================================



#=========================================================
