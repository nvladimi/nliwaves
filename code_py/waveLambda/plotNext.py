# this is a comment

import kinwave as kw
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

#--------------------------------------------


outfile = "extN_pi2.pdf";  fbase = "2023-04-08/extN_pi2_"
#outfile = "extN_pi4.pdf";  fbase = "2023-04-08/extN_pi4_"



runs = (("n02m04", "n02m08"),
        ("n04m08", "n04m16")
)

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6.2,6.2)) 

vrange = 10000
 
for i in (0,1):
    for j in (0,1):

        fname = fbase + runs[i][j] + ".npy"
        with open(fname, 'rb') as f:
            p   = np.load(f)
            p1  = np.load(f)
            (mq, dq, mr, dr) = np.load(f) 
            N2D = np.load(f)
            m = int ( (N2D.shape[0]-1)/2 )
    
            
  
        iax = ax[i,j]
        
        extent=( (-m-0.5)*dq, (m+0.5)*dq, (-m-0.5)*dq, (m+0.5)*dq )
    
        im0=iax.imshow(np.transpose(N2D), cmap="jet", origin = 'lower',
                       vmin=-vrange, vmax=vrange, extent=extent)

        iax.hlines(y=0, xmin=-dq/2, xmax=dq/2, linewidth=0.25, color='k')
        iax.vlines(x=0, ymin=-dq/2, ymax=dq/2, linewidth=0.25, color='k')

        iax.set_title( "$m_q = {}, \quad \Delta q = 1/{}$".format(int(mq), int(1/dq)) )

plt.subplots_adjust(left=0.04, right=0.99, top=0.95, bottom=0.05, hspace=0.2, wspace=0.08)

plt.savefig(outfile)

plt.close()

print("Done.")


#=========================================================
