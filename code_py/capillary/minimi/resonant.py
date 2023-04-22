#import importlib ; import resonant
#importlib.reload(resonant); resonant.Test()

def Test(kmax=10, fnameout='test.pdf', debug=False, showplot=False):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.rcParams['text.usetex'] = True

    krange = np.arange(kmax)
    numTriangles  = np.zeros(kmax)

    for k in krange:
        allT = Triangles(k, debug=False)
        numTriangles[k] = len(allT)

    if showplot:
        fig, iax = plt.subplots(ncols=1, nrows=1, figsize=(6,4))

        iax.set_title("Number of Triangles in box N x N");
        iax.set_xscale('log')
        iax.set_yscale('log')
        
        iax.set_xlabel('$N$')
        iax.set_ylabel('triangles')
        
        iax.set_xlim(2,100)
        iax.set_ylim(1,1e5)
 
        n = krange*2+1
        
        iax.plot(n, numTriangles, '-or', markersize=3)
        iax.plot(n, 0.2*n**3, '--k', label='$N^3$')
        iax.plot(n, n**(7/3.), '--b', label='$N^{7/3}$')

        iax.legend()

        
        #plt.savefig('disk_zoom_select.png', pad_inches=0, dpi=figDPI)
        plt.tight_layout()
        plt.show()
  


#============================================================

def Modes(kmax, m=(0,0), fnameout='test.pdf', debug=False, showplot=False):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm

    N = 2*kmax+1
    
    All  = np.zeros([N,N], dtype=int)
    Diag = np.zeros([N,N], dtype=int)
 
    allT = Triangles(kmax, debug=False)

    for T in allT:

        for v in T:
            All[v[0], v[1]] += 1           
        i,j = hypothenuse(T)
        Diag[i][j] += 1

        if debug:
            if m in T:
                print(T)

    All  = np.fft.fftshift(All)
    Diag = np.fft.fftshift(Diag)
    Side = All - Diag

    if debug:
    
        print("\n Diagonals:")
        print(Diag)
    
        print("\n Sides:")
        print(Side)
    
        print("\n All:")  
        print(All)
    
    print("\n Totals (triangles, diagonals, sides, all):",
          len(allT), np.sum(Diag), np.sum(Side), np.sum(All))
        

    #-- plot ---

    if showplot:
    
        fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(10,3))
        extent=[-kmax-0.5, kmax+0.5, -kmax-0.5, kmax+0.5]

        vmax = np.max(All)+1 

        jet = cm.get_cmap('jet', vmax)

        im = ax[0].imshow(Side, cmap=jet, vmin=0, vmax=vmax, origin='lower', extent=extent);
        im = ax[1].imshow(Diag, cmap=jet, vmin=0, vmax=vmax, origin='lower', extent=extent);
        im = ax[2].imshow(All,  cmap=jet, vmin=0, vmax=vmax, origin='lower', extent=extent);

        ax[0].set_title("Sides");
        ax[1].set_title("Diagonals");
        ax[2].set_title("Overall");

        cbar = fig.colorbar(im, extend='both', shrink=1,  ax=ax[2])
    
        #plt.savefig('disk_zoom_select.png', pad_inches=0, dpi=figDPI)
        plt.tight_layout()
        plt.show()

    return(Diag, Side)
    
        

#============================================================
def hypothenuse(T):

    x1 = T[0][0]
    y1 = T[0][1]

    x2 = T[1][0]
    y2 = T[1][1]
    
    x3 = T[2][0]
    y3 = T[2][1]

    rr1 = x1*x1 + y1*y1
    rr2 = x2*x2 + y2*y2
    rr3 = x3*x3 + y3*y3

    if (rr1 > rr2) and (rr1 > rr3):
        return(x1,y1)
    if (rr2 > rr1) and (rr2 > rr3):
        return(x2,y2)
    if (rr3 > rr1) and (rr3 > rr2):
        return(x3,y3)

    return(0)


    
#============================================================
def Triangles(kmax, debug=False):

    import numpy as np

    k = np.arange(-kmax, kmax+1)

    kx, ky = np.meshgrid(k,k)

    k1 = kx + 1.j*ky
    k2 = k1
    
    allT = set({})
   
    for q in k1:
        for q1 in q:

            if q1 == 0:
                continue
            
            c = np.real(k2 * np.conj(q1))
            ind = ( (c == 0) & (k2 != 0))
           
            k3 = k2[ind] + q1
            for q3 in k3:

                x3 = int (np.real(q3))
                y3 = int (np.imag(q3))
            
                if ( x3 >= -kmax and  x3 <=  kmax and  y3 >= -kmax and y3 <=  kmax):

                    x1 = int(np.real(q1))
                    y1 = int(np.imag(q1))
                               
                    T=[(x1, y1), (x3-x1, y3-y1), (x3, y3)]
                    T.sort()
                    #T.append(int(q3*np.conj(q3)))
                    allT.add(tuple(T))

    if debug:
        print("\n Triangles:")
        for T in allT:
            
            print(T)
        
    print("\n Box {} x {} total triangles: {}".format(2*kmax+1, 2*kmax+1, len(allT)) )

    return(allT)

#============================================================

