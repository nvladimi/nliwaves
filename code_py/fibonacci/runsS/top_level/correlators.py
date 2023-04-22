import numpy as np

mmax=8; cmax=1000
#mmax=7; cmax=200
#mmax=6; cmax=50
#mmax=5; cmax=15

#-- find all --

corr = np.zeros((cmax, mmax), dtype=int)
count = np.zeros(cmax, dtype=int)
currcorr = np.zeros(mmax, dtype=int)

corr[0][0]=1
corr[0][1]=1
count[0]  = 1

c1=0
c2=1

while c1 < cmax:
    #print("\n c1=", c1, corr[c1])
    for i in range(mmax-2):
        if ((corr[c1][i] > 0) and (c2 < cmax)) :
            np.copyto(currcorr, corr[c1])
            currcorr[i]   -= 1
            currcorr[i+1] += 1
            currcorr[i+2] += 1
            #print("currcorr=", currcorr)
            for j in range(c2+1):
                if np.all(corr[j]-currcorr == 0):
                    count[j] += count[c1]
                    break
            #print("j=", j, "count=", count[j])
            if (j == c2):
                corr[c2] = currcorr
                count[c2] += count[c1]
                #print("c2=", c2, corr[c2], count[c2])
                c2 += 1
    c1 += 1

imax = np.max(np.where(count > 0))

if (imax+1 == cmax):
    print("Increase cmax.")
    quit()
    
for i in range(imax+1):
    c = corr[i]
    s = ''
    m = np.sum(c)+1
    #if (m>9):
    #    break
    for k in range(len(c),0,-1):
        for j in range(c[k-1]):
            s += str(k)
    n=0
    for ss in s:
        n += int(ss)
            
    print(c, '   ', s,  '    ', m, '    ', n, '   ', count[i])

#--------------------------------------------------------------
quit()
#--------------------------------------------------------------
#-- slow, but nore explicit way --

#-- find all --

corr = np.zeros((cmax, mmax), dtype=int)
corr[0][0]=1
corr[0][1]=1

c1=0
c2=1

while c1 < cmax:
    for i in range(mmax-2):
        if ((corr[c1][i] > 0) and (c2 < cmax)) :
            corr[c2] = corr[c1]
            corr[c2][i]   -= 1
            corr[c2][i+1] += 1
            corr[c2][i+2] += 1
            c2 += 1
    c1 += 1

print(corr)


    
if (np.max(corr[-1]) > 0):
    print("Increase cmax.")
    quit()

#-- find multiples --
    
summary = np.zeros((cmax, mmax), dtype=int)
count = np.zeros(cmax, dtype=int)
currcorr = np.zeros(mmax, dtype=int)

c2=0
for c1 in range(cmax):
    if np.any(corr[c1] > 0):
        np.copyto(currcorr, corr[c1])
        for i in range(c1, cmax):
            if np.all(corr[i]-currcorr == 0):
                count[c2] += 1
                corr[i][:] = 0
        summary[c2] = currcorr
        c2 += 1

imax = np.max(np.where(count > 0))

for i in range(imax+1):
    print(summary[i], count[i])

