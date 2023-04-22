import numpy as np
import itertools

F = np.array([1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 932 ])

filename="correlators_m7_found.txt"

with open(filename) as file:
    lines = file.readlines()
    #lines = [line.rstrip() for line in lines]

for line in lines:
    cc=line.rstrip()
    if cc[0] == "T":
        continue
    cc=cc[1:-1].split(",")
    jmax = len(cc)
    cp=""
    cm=""
    for j in range(0,jmax):
        d = int(cc[j])
        if d<0:
            cm = cm + str(jmax-j-1)*(-d)
        if d>0:
            cp = cp + str(jmax-j-1)*d
        
    #print(type(c))
    print(cp+"-"+cm)
    
