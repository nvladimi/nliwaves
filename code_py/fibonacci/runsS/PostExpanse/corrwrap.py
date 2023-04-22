
import sys
import fiboPost

outdir="Post06"
f="s5m10"
cc = sys.argv[1]

fbaseout = outdir + '/' + f + "_" + cc

cs = cc.split("-")
cp = []
cm = []
for c in cs[0]:
    cp.append(int(c))
for c in cs[1]:
    cm.append(int(c))
    
print(f, cp, cm, fbaseout)

fiboPost.modeCorr(f, cp, cm,  fbaseout = fbaseout)


