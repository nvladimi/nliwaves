
import matplotlib
#matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import numpy as np
from matplotlib.colors import LinearSegmentedColormap


Ta = [[
[[  1,   4],  [  4,  -1]],
[[  1,  -4],  [  4,   1]],
[[ -1,   4],  [  4,   1]],
[[ -4,  -1],  [  1,  -4]],
],[
[[ -4,   1],  [  1,   4]],
[[ -4,  -1],  [ -1,   4]],
[[ -1,  -4],  [  4,  -1]],
[[ -4,   1],  [ -1,  -4]],
],[
[[ -3,  -3],  [ -3,   3]],
[[ -3,  -3],  [  3,  -3]],
[[ -3,   3],  [  3,   3]],
[[  3,  -3],  [  3,   3]],
],[
[[  0,   7],  [  7,   0]],
[[ -7,   0],  [  0,   7]],
[[  0,  -7],  [  7,   0]],
[[ -7,   0],  [  0,  -7]],
]]


Tb = [[
[[ -1,   7],  [  7,   1]],
[[ -7,  -1],  [  1,  -7]],
[[ -7,   1],  [ -1,  -7]],
[[ -7,  -1],  [ -1,   7]],
],[
[[  1,   7],  [  7,  -1]],
[[  1,  -7],  [  7,   1]],
[[ -7,   1],  [  1,   7]],
[[ -1,  -7],  [  7,  -1]],
],[
[[  5,  -5],  [  5,   5]],
[[ -5,  -5],  [ -5,   5]],
[[ -5,  -5],  [  5,  -5]],
[[ -5,   5],  [  5,   5]],
],[
[[  2,  -4],  [ -8,  -4]],
[[  2,   4],  [ -8,   4]],
[[ -2,   4],  [  8,   4]],
[[ -4,   2],  [ -4,  -8]],
],[
[[  4,  -2],  [  4,   8]],
[[  4,   2],  [  4,  -8]],
[[ -4,  -2],  [ -4,   8]],
[[ -2,  -4],  [  8,  -4]],
]]




#prefix = "../pyPost/c1"; fnameout = "img_triple_c1a.pdf"; cmax=0.010; N=60; T=Ta
#prefix = "../pyPost/c1"; fnameout = "img_triple_c1b.pdf"; cmax=0.010; N=60; T=Tb
#prefix = "../pyPost/c5"; fnameout = "img_triple_c5a.pdf"; cmax=0.010; N=60; T=Ta
#prefix = "../pyPost/c5"; fnameout = "img_triple_c5b.pdf"; cmax=0.010; N=60; T=Tb
#prefix = "../pyPost/e1"; fnameout = "img_triple_e1a.pdf"; cmax=0.020; N=40; T=Ta
prefix = "../pyPost/e1"; fnameout = "img_triple_e1b.pdf"; cmax=0.020; N=40; T=Tb



#------------------------------------------

ucolors = [(0.95, 0.95, 0.95), (0.9, 0.8, 0), (0.7, 0, 0)]  # R -> G -> B
ucmap = LinearSegmentedColormap.from_list("ucolors", ucolors, N=4)

mmax=12
i1 = int(N/2-mmax)
i2 = int(N/2+mmax+1)
extent=(-mmax-0.5, mmax+0.5,-mmax-0.5, mmax+0.5)

nrows=len(T)
ncols=len(T[0])
fsize=4
ms=10
#-------------------------------

fig, ax = plt.subplots(ncols=ncols, nrows=nrows, constrained_layout=True, figsize=(ncols*fsize,nrows*fsize))

for j in range(nrows):
    for i in range(ncols):

        iax = ax[j][i]
        C = T[j][i]

        corr = "{:+03d}{:+03d}{:+03d}{:+03d}".format(C[0][0], C[0][1], C[1][0], C[1][1])
        Csum = np.sum(C,0)

        #-- data --    
        a = np.load(prefix + corr + ".npy")
        c = np.fft.fftshift(a);
        c = c[i1:i2, i1:i2]

        im0=iax.imshow(np.abs(c),  cmap=ucmap, origin = 'lower', vmin=0,     vmax=cmax, extent=extent);

        #-- annotations --

        iax.set_title("[{}, {}] + [{}, {}] = [{}, {}]: {:.3f}".format(
            C[0][0], C[0][1], C[1][0], C[1][1],
            Csum[0], Csum[1], np.abs(a[Csum[1],Csum[0]]) ) )

        kmodes=( (C[0][0], C[1][0], Csum[0]), (C[0][1], C[1][1], Csum[1]))

        iax.plot(kmodes[0], kmodes[1], 'sk', ms=ms, linewidth=0.5, mfc='none')
        iax.plot(0, 0, 'xk', ms=ms, linewidth=0.5)


#------------------------------

"""    
cbar0=fig.colorbar(im0, ax=ax[j,0], orientation='horizontal', shrink=0.9)
cbar0.set_ticks([0,cmax/2, cmax])
cbar0.set_ticklabels([0,cmax/2, cmax])
"""
plt.savefig(fnameout, pad_inches=0)

#plt.show()
#plt.close()

exit()




#------------------------------------------
