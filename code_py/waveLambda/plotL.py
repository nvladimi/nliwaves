# this is a comment

import kinwave as kw
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

#-------------------------------------------

def write_L2D():

    #m=10; dq=0.1; p0=1; fname = "lambda_m10_dq01_p1.npy"
    #m=10; dq=0.1; p0=10; fname = "lambda_m10_dq01_p10.npy"
    #m=10; dq=0.1; p0=2; fname = "lambda_m10_dq01_p2.npy"
    #m=10; dq=0.1; p0=0.5; fname = "lambda_m10_dq01_p05.npy"
    m=10; dq=0.1; p0=0.1; fname = "lambda_m10_dq01_p01.npy"

    pp = np.array(
        ( [ 1.0,  0.0], [ 1.0,  0.4],  [ 1.0,  1.0], [ 0.4,  1.0],
          [ 0.1,  1.0], [-0.4,  1.0],  [-1.0,  1.0], [-1.0,  0.4],
          [-1.0,  0.1], [-1.0, -0.4],  [-1.0, -1.0], [-0.4, -1.0],
          [ 0.1, -1.0], [ 0.4, -1.0],  [ 1.0, -1.0], [ 1.0, -0.4]) )


    kw.L2Dfun(m, dq, pp*p0, fname)


#write_L2D()



#--------------------------------------------

def write_L4d():

    m=4; mr=2;  fname="lambda_04o02.npy"
    #m=8; mr=4; fname="lambda_08o04.npy"
    #m=16; mr=8; fname="lambda_16o08.npy"

    kw.L4Dfun(m, mr, fname)

    with open(fname, 'rb') as f:
        [m, mr] = np.load(f)
        L4D = np.load(f)

    print(m, mr, L4D.shape)

#--------------------------------------------

def plot_L2D():

    fname = "lambda_m10_dq01_p01.npy"
    outfile = "plotL.pdf"

    with open(fname, 'rb') as f:
        [m, dq] = np.load(f)
        pp = np.load(f)
        L2D = np.load(f)
        m = int(m)
        
    #print(m, dq, pp, L2D.shape)
    vrange = np.max( [np.abs(np.max(L2D)), np.abs(np.min(L2D))] )
    print(vrange)

    fig, ax = plt.subplots(ncols=5, nrows=5, figsize=(10,10)) 

    place =   ( [ 2,  4], [ 1, 4],  [ 0, 4], [ 0, 3],
                [ 0,  2], [ 0, 1],  [ 0, 0], [ 1, 0],
                [ 2,  0], [ 3, 0],  [ 4, 0], [ 4, 1],
                [ 4,  2], [ 4, 3],  [ 4, 4], [ 3, 4] )

    empty =  ([1,1], [1,2], [1,3],   [2,1], [2,2], [2,3],  [3,1], [3,2], [3,3])
              
    

    for i in range(16):

        Q = L2D[i,:,:]
        #Q=np.sign(Q); vrange=2

        iax = ax[place[i][0], place[i][1]]
       
        im0=iax.imshow(np.transpose(Q), cmap="jet", origin = 'lower', vmin=-vrange, vmax=vrange)
        #extent=extent);
        iax.set_xticks([])
        iax.set_yticks([])
        iax.hlines(y=m, xmin=0, xmax=2*m, linewidth=0.25, color='k')
        iax.vlines(x=m, ymin=0, ymax=2*m, linewidth=0.25, color='k')


    #-- information in the middle --
        
    for i in range(9):
        iax = ax[empty[i][0], empty[i][1]]
        iax.axis('off')

    iax = ax[1,2]

    iax.text(0,0.95, "$\mathbf{r} = [0,1]$")
    iax.text(0,0.80, "${} < q_x, q_y < {}$ in each square".format(-m*dq,m*dq) )
    iax.text(0,0.65, "$\lambda$ is shown for $\mathbf{p}$: "  )

    iax.text(0, 0.5,
             str(pp[6]) + "  " +
             str(pp[5]) + "  " +
             str(pp[4]) + "  " +
             str(pp[3]) + "  " +
             str(pp[2]) )

    iax.text(0,   0.4,   str(pp[7]))
    iax.text(1.3, 0.4,   str(pp[1]))

    iax.text(0,   0.3,   str(pp[8]))
    iax.text(1.3, 0.3  , str(pp[0]))

    iax.text(0,   0.2,   str(pp[9]))
    iax.text(1.3, 0.2  , str(pp[15]))


    iax.text(0, 0.1,
             str(pp[10]) + "  " +
             str(pp[11]) + "  " +
             str(pp[12]) + "  " +
             str(pp[13]) + "  " +
             str(pp[14]) )


    iax.text(0, -0.1, "colormap limits = $\pm {}$".format(vrange) )

    #title = "p = [a, -1],  \ \ \ q = [ {}, {} ], \ \ \  r = [0, 1e-5]".format( q[0], q[1])
    #ax[0,0].set_title(title)

    #-----------------------

    #plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)


    plt.savefig(outfile) #pad_inches=0

    plt.close()


plot_L2D()

print("Done.")




#===========================================



#=========================================================
