# this is a comment

import kinwave as kw
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

#--------------------------------------------


outfile = "rezN_pi2.pdf";  fbase = "2023-04-08/resN_pi2_"
#outfile = "rezN_pi4.pdf";  fbase = "2023-04-08/resN_pi4_"



runs = (("n04m08", "n08m16", "n16m32", "n32m64" ),
        ("n04m12", "n08m24", "n16m48", "n32m96" )
)

fig, ax = plt.subplots(ncols=4, nrows=2, figsize=(11.2,6.2)) 

vrange = 10000
 
for i in (0,1):
    for j in (0,1,2,3):

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

        iax.set_title( "$m_r = {}, \quad \Delta r = 1/{}$".format(int(mr), int(1/dr)) )

plt.subplots_adjust(left=0.02, right=0.99, top=0.95, bottom=0.05, hspace=0.2, wspace=0.15)

plt.savefig(outfile)

plt.close()

print("Done.")


#=========================================================
