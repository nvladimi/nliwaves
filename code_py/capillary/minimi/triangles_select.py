import numpy as np
import resonant

kmin=4
kmax=10

#-- all triangles with max vector shorter than kmax --

allT = resonant.Triangles(kmax, debug=False)


#-- remove triangles with any vector shorter than kman --

T = np.array(list(allT)).squeeze()
kk = (T[:,:,0]**2 + T[:,:,1]**2).squeeze()
kkmin = np.min(kk,1)

ind = np.argwhere(kkmin>kmin**2).squeeze()

T = T[ind]
kk = kk[ind]


#-- sort triangles by largest vector --

kkmax = np.max(kk,1).squeeze()
ind = np.argsort(kkmax)

T = T[ind]
kk = kk[ind]


#-- in each trianble, sort vectors by amplitude ---

for i in np.arange(len(kk)):
    t = T[i]
    ind = np.argsort(kk[i])
    T[i] = t[ind]

#-- output --
print("above kmin={}:  {}".format(kmin, len(T)) )  
for i in np.arange(len(kk)):
    t = T[i]
    kk = (t[:,0]**2 + t[:,1]**2)
    print("[{:3d} {:3d}]  [{:3d} {:3d}]  [{:3d} {:3d}]    [{:3d} {:3d}  {:3d}] ".
          format(t[0][0], t[0][1],  t[1][0], t[1][1],  t[2][0], t[2][1],  kk[0], kk[1], kk[2]) )
    
