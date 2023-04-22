import numpy as np
import itertools

F = np.array([1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 932 ])

filename="correlators_m7_sorted.txt"

#   fiboPost.modeCorr(f, [3, 3, 1], [5,0],    fbaseout = outdir + '/' + f + '_331-50');

with open(filename) as file:
    lines = file.readlines()
    #lines = [line.rstrip() for line in lines]

for line in lines:
    cc=line.rstrip()    
    if cc == "":
        continue
    if cc[0] == "#":
        continue
    cc = cc.split()
    cc = cc[0]
    cs = cc.split("-")
    cp="["
    cm="["
    for c in cs[0]:
        cp = cp + c + ","
    for c in cs[1]:
        cm = cm + c + ","
    s1 = "fiboPost.modeCorr(f, "
    s2 = cp + "], " + cm  + "]"
    s3 = ",   fbaseout = outdir + '/' + f + '_" + cc + "')"
    print(s1+s2+s3)
