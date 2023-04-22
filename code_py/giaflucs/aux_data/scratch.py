
import os
import matplotlib.pyplot as plt
import numpy as np
import pickle


import importlib
import giaflucs ; giaflucs.ringImax()
importlib.reload(giaflucs)


with open("./slice_3ring.pkl", "rb") as f: data = pickle.load(f)

type(data)
#  <class 'list'>

len(data)
# 2

type(data[0])
# <class 'numpy.ndarray'>

type(data[1])
#  <class 'numpy.ndarray'>

len(data[0])
#  57518

len(data[1])
#  57518

>>> data[0].dtype
dtype('float64')

>>> data[0].shape
(57518, 3)

>>> data[0].ndim
2

>>> data[0].size
172554


>>> data[0].imag
array([[0., 0., 0.],
       [0., 0., 0.],
       [0., 0., 0.],
       ...,
       [0., 0., 0.],
       [0., 0., 0.],
       [0., 0., 0.]])

>>> data[0].real
array([[1.71050000e+04, 1.80650000e+04, 8.09010491e-05],
       [1.71490000e+04, 1.88560000e+04, 8.09031914e-05],
       [1.70950000e+04, 1.78730000e+04, 8.09037891e-05],
       ...,
       [1.72250000e+04, 1.89380000e+04, 5.39578718e-04],
       [1.72260000e+04, 1.89380000e+04, 5.44251946e-04],
       [1.72260000e+04, 1.89390000e+04, 5.60422337e-04]])


>>> data[1].dtype
dtype('int64')

>>> data[1].size
57518

>>> data[1].shape
(57518,)


#-----------------

disk.npy
disk_3.npy
disk_3full.npy

painting.ipynb

scratch.py

slice16.pkl
slice32.pkl
slice8.pkl
slice_3ring.pkl

#------------------


>>> actual_arr = np.load("disk_3.npy")
>>> actual_arr.shape
(3025, 3025)

>>> actual_arr = np.load("disk_3full.npy")
>>> actual_arr.shape
(3025, 3025)

>>> actual_arr = np.load("disk.npy")
>>> actual_arr.shape
(1493, 1493)

>>> plt.imshow(np.log(actual_arr), cmap='jet', vmin=-12, vmax=-7)
<matplotlib.image.AxesImage object at 0x7ff618044bd0>


#------------------


>>> import matplotlib.pyplot as plt
>>> import numpy as np
>>> ig, ax = plt.subplots()
>>> ax.plot([1, 2, 3, 4], [1, 4, 2, 3]) 
[<matplotlib.lines.Line2D object at 0x7ff9b0178790>]

>>> 
>>> plt.plot([1, 2, 3, 4])
[<matplotlib.lines.Line2D object at 0x7ff9d02fc210>]
>>> plt.ylabel('some numbers')
Text(0, 0.5, 'some numbers')

>>> plt.show()




>>> import os
>>> import matplotlib.pyplot as plt
>>> import numpy as np
>>> actual_arr = np.load("disk.npy")
>>> plt.imshow(np.log(actual_arr), cmap='jet')
<matplotlib.image.AxesImage object at 0x7fed28b8f2d0>
>>> 
>>> plt.show()



from numpy.random import default_rng
rng = default_rng()
a = abs(rng.standard_normal(25).reshape(5,5));
b=np.ones((5,5))
c = b*(a>1)


ind = np.nonzero(a>1)
b = np.zeros((5,5))
b[ind]=1


q=(a>1)

b = np.round(a*100)

q = np.random.rand(3,4)


q=a[:,0:-1]


q=a[:,1:]


b0 = b[1:-1, 1:-1]
bE = b[1:-1, 2:]
bW = b[1:-1, 0:-2]
bN = b[0:-2, 1:-1]
bS = b[2:,   1:-1]


q0 = ( (b0>bN) & (b0>bS) & (b0>bE) & (b0>bW) & (b0>bThershold) )

q = np.full((5, 5), False)

q[1:-1,  1:-1] = q0

ind = np.where(q)
bmax = b[ind]



a = np.random.rand(5,5)
q = (a>0.5)


matplotlib.pyplot.savefig('filename')

plt.ion()
plt.plot([1.6, 2.7])
ax = plt.gca()
ax.plot([3.1, 2.2])

#------------------------------------


from copy import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# compute some interesting data
x0, x1 = -5, 5
y0, y1 = -3, 3
x = np.linspace(x0, x1, 500)
y = np.linspace(y0, y1, 500)
X, Y = np.meshgrid(x, y)
Z1 = np.exp(-X**2 - Y**2)
Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
Z = (Z1 - Z2) * 2

# Set up a colormap:
# use copy so that we do not mutate the global colormap instance
palette = copy(plt.cm.gray)
palette.set_over('r', 1.0)
palette.set_under('g', 1.0)
palette.set_bad('b', 1.0)
# Alternatively, we could use
# palette.set_bad(alpha = 0.0)
# to make the bad region transparent.  This is the default.
# If you comment out all the palette.set* lines, you will see
# all the defaults; under and over will be colored with the
# first and last colors in the palette, respectively.
Zm = np.ma.masked_where(Z > 1.2, Z)

# By setting vmin and vmax in the norm, we establish the
# range to which the regular palette color scale is applied.
# Anything above that range is colored based on palette.set_over, etc.

# set up the Axes objects
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(6, 5.4))

# plot using 'continuous' color map
im = ax1.imshow(Zm, interpolation='bilinear',
                cmap=palette,
                norm=colors.Normalize(vmin=-1.0, vmax=1.0),
                aspect='auto',
                origin='lower',
                extent=[x0, x1, y0, y1])
ax1.set_title('Green=low, Red=high, Blue=masked')
cbar = fig.colorbar(im, extend='both', shrink=0.9, ax=ax1)
cbar.set_label('uniform')
for ticklabel in ax1.xaxis.get_ticklabels():
    ticklabel.set_visible(False)

# Plot using a small number of colors, with unevenly spaced boundaries.
im = ax2.imshow(Zm, interpolation='nearest',
                cmap=palette,
                norm=colors.BoundaryNorm([-1, -0.5, -0.2, 0, 0.2, 0.5, 1],
                                         ncolors=palette.N),
                aspect='auto',
                origin='lower',
                extent=[x0, x1, y0, y1])
ax2.set_title('With BoundaryNorm')
cbar = fig.colorbar(im, extend='both', spacing='proportional',
                    shrink=0.9, ax=ax2)
cbar.set_label('proportional')

fig.suptitle('imshow, with out-of-range and masked data')
plt.show()

#############################################################################

from copy import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

plt.ion()

x0, x1 = -5, 5 ;  y0, y1 = -3, 3 ;  x = np.linspace(x0, x1, 500); y = np.linspace(y0, y1, 500)
X, Y = np.meshgrid(x, y) ;  Z1 = np.exp(-X**2 - Y**2); Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2);  Z = (Z1 - Z2) * 2

palette = copy(plt.cm.gray); palette.set_over('r', 1.0); palette.set_under('g', 1.0); palette.set_bad('b', 1.0)

Zm = np.ma.masked_where(Z > 1.2, Z)


fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(6, 5.4))


im = ax1.imshow(Zm, interpolation='bilinear', cmap=palette,
                norm=colors.Normalize(vmin=-1.0, vmax=1.0),
                aspect='auto', origin='lower', extent=[x0, x1, y0, y1])


ax1.set_title('Green=low, Red=high, Blue=masked')

cbar = fig.colorbar(im, extend='both', shrink=0.9, ax=ax1)

cbar.set_label('uniform')
for ticklabel in ax1.xaxis.get_ticklabels():
    ticklabel.set_visible(False)


im = ax2.imshow(Zm, interpolation='nearest',  cmap=palette,
                norm=colors.BoundaryNorm([-1, -0.5, -0.2, 0, 0.2, 0.5, 1], ncolors=palette.N),
                aspect='auto', origin='lower',  extent=[x0, x1, y0, y1])

ax2.set_title('With BoundaryNorm')

cbar = fig.colorbar(im, extend='both', spacing='proportional', shrink=0.9, ax=ax2)
cbar.set_label('proportional')

fig.suptitle('imshow, with out-of-range and masked data')
plt.show()


#############################################################################


from copy import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

plt.ion()

x0, x1 = -5, 5 ;  y0, y1 = -3, 3 ;  x = np.linspace(x0, x1, 500); y = np.linspace(y0, y1, 500)
X, Y = np.meshgrid(x, y) ;  Z1 = np.exp(-X**2 - Y**2); Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2);  Z = (Z1 - Z2) * 2

palette = copy(plt.cm.gray); palette.set_over('r', 1.0); palette.set_under('g', 1.0); palette.set_bad('b', 1.0)

Zm = np.ma.masked_where(Z > 1.2, Z)

fig, ax = plt.subplots()

im =  ax.imshow(Zm, interpolation='bilinear', cmap=palette,
                norm=colors.Normalize(vmin=-1.0, vmax=1.0),
                aspect='auto', origin='lower', extent=[x0, x1, y0, y1])


im =  ax.imshow(Zm, cmap=palette, norm=colors.Normalize(vmin=-1.0, vmax=1.0),  origin='lower')


im = ax.imshow(Zm, interpolation='nearest',  cmap=palette,
                norm=colors.BoundaryNorm([-1, -0.5, -0.2, 0, 0.2, 0.5, 1], ncolors=palette.N),
                aspect='auto', origin='lower',  extent=[x0, x1, y0, y1])


plt.show()


#############################################################################



import numpy as np
import matplotlib.pyplot as plt

dat = np.loadtxt("chi4a_t10150-10160.txt")


n1 = dat[:,1]*dat[:,1] + dat[:,2]*dat[:,2] ;  n2 = dat[:,3]*dat[:,3] + dat[:,4]*dat[:,4] ; t = dat[:,0]
R = 0.001; gamma = 0.1; nu1 = 2*R*R/gamma

plt.ion()

fig, ax = plt.subplots(); ax.plot(t,n1/nu1); ax.plot(t,n2/nu1);













