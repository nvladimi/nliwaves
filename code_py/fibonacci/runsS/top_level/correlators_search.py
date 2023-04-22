import numpy as np
import itertools

F = np.array([1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 932 ])
#F = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192])


def partition(number):
    answer = set()
    answer.add((number, ))
    for x in range(1, number):
        for y in partition(number - x):
            answer.add(tuple(sorted((x, ) + y)))
    return answer


c3 = ( [1,1,1,], [2,1,0],) # [2,0,0])
c4 = ( (1,1,1,1), (2,1,1,0), (2,2,0,0), (3,1,0,0),) # (4,0,0,0) )
c5 = ( (1,1,1,1,1), (2,1,1,1,0), (2,2,1,0,0), (3,1,1,0,0), (3,2,0,0,0), (4,1,0,0,0), ) # (5,0,0,0,0) )
c6 = ( (1,1,1,1,1,1), (2,1,1,1,1,0), (3,1,1,1,0,0), (4,1,1,0,0,0), (5,1,0,0,0,0), #(6,0,0,0,0,0), 
       (2,2,1,1,0,0), (2,2,2,0,0,0), (3,2,1,0,0,0), (3,3,0,0,0,0), (4,2,0,0,0,0) )

c7=  ( (1,1,1,4,0,0,0), (1,1,1,1,1,2,0), (2,2,3,0,0,0,0), (1,1,1,1,3,0,0), (1,1,1,1,1,1,1),
       (2,5,0,0,0,0,0), (1,1,2,3,0,0,0), (1,6,0,0,0,0,0), (1,2,2,2,0,0,0), (1,1,1,2,2,0,0), 
       (1,1,5,0,0,0,0), (1,2,4,0,0,0,0), (3,4,0,0,0,0,0), (1,3,3,0,0,0,0) ) #(7,0,0,0,0,0,0)

#m = 5; kmax=6; cc = c5
#m = 6; kmax=8; cc = c6
m = 7; kmax=10; cc = c7


extc=list(itertools.repeat(0, kmax-m))

ss = []
for i in range(m):
    s1 = list(itertools.repeat(1, m-i))
    s2 = list(itertools.repeat(-1, i))
    s1.extend(s2)
    perm = set(itertools.permutations(s1))
    for i in perm:
        ss.append(list(i))

#print(ss)

allC = []
for c in cc:
    for s in ss:
        ci = list(np.array(c)*np.array(s))
        ci.extend(extc)
        perm = set(itertools.permutations(ci))
        for C in perm:
            FC = np.sum(np.array(C) * F[3:kmax+3])
            if FC == 0:
                ac = np.array(C)
                ind = np.nonzero(ac)[0];
                ac= ac[np.min(ind):np.max(ind)+1]
                if(ac[-1]>0):
                    ac = -ac
                #print(ac)
                allC.append(list(ac))

allC.sort()
allC = list(allC for allC, _ in itertools.groupby(allC))

for c in allC:
    print(c)

print( "Total: ", len(allC) ) 





#--------------------------------------------------------------
quit()
#--------------------------------------------------------------




k=3;

print({{1, 2}, {2, 3}})
print([{4, 2}, {5, 3}])

allsets = [{}, {1}, {{1,1},{2}} ]

def make_sumsets(m, allsets):
    print("   m = ", m)
    msets = {}
    for k in range(m):
        print("   k = ", k)
        print("   allsets[k] = ",  allsets[k])
        for s in allsets[k]:
            print("       s = ", s)
            s0 = s.copy()
            print("       s0 = ", s0)
            s0.add(m-k)
            print("       s0 = ", s0)
            msets.add(s0)
            print("       msets = ", msets)
    return(msets)

for i in range(1,m):
    print("i = ", i)
    msets = make_sumsets(i, allsets)
    print("msets = ", msets)
    allsets.append(msets)
    






