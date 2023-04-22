import sys
import fiboPost

fbase = sys.argv[1]
fset  = sys.argv[2]

iset  = int(fset)


fbaseout = "POSTmts/" + fbase + "_set" + fset

istart = iset*100 + 1
iend   = istart + 100

fiboPost.Moments(fbase, nmom=12, istart=istart, iend=iend, fbaseout=fbaseout)




