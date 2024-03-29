import numpy as np
import nlwTools as nlw
import sys
run=sys.argv[1]
outdir=sys.argv[2]

#run="d1"
#ourdir="PostA"

#=========================================================

if run=="c1":
    N=60
    No=128
    seeds = (10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
    istart = 21
    iend = -1
    
if run=="c5":
    N=60
    No=128
    seeds = (10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29)
    istart = 11
    iend = -1


if run=="d1":
    N=40
    No=64
    seeds = (10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
    istart = 21
    iend = -1


if run=="d5":
    N=40
    No=64
    seeds = (10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
    istart = 21
    iend = -1

    
#=========================================================
#triplets:

T = [
[[  0,   7],  [  7,   0]],
[[ -7,   0],  [  0,   7]],
[[  0,  -7],  [  7,   0]],
[[ -7,   0],  [  0,  -7]],
[[ -1,   7],  [  7,   1]],
[[ -7,  -1],  [  1,  -7]],
[[ -7,   1],  [ -1,  -7]],
[[ -7,  -1],  [ -1,   7]],
[[  1,   7],  [  7,  -1]],
[[  1,  -7],  [  7,   1]],
[[ -7,   1],  [  1,   7]],
[[ -1,  -7],  [  7,  -1]],
[[  5,  -5],  [  5,   5]],
[[ -5,  -5],  [ -5,   5]],
[[ -5,  -5],  [  5,  -5]],
[[ -5,   5],  [  5,   5]],
[[  2,  -4],  [ -8,  -4]],
[[  2,   4],  [ -8,   4]],
[[ -2,   4],  [  8,   4]],
[[ -4,   2],  [ -4,  -8]],
[[  4,  -2],  [  4,   8]],
[[  4,   2],  [  4,  -8]],
[[ -4,  -2],  [ -4,   8]],
[[ -2,  -4],  [  8,  -4]],
[[  3,   4],  [  4,  -3]],
[[  0,   5],  [  5,   0]],
[[  3,   3],  [ -4,   4]],
[[  4,  -4],  [  4,   4]],
[[  0,   6],  [  6,   0]],
[[  0,   6],  [  7,   0]],
[[  1,   4],  [  8,  -2]],
[[  1,   4],  [ -8,   2]],
[[  0,   7],  [  7,   0]],
[[  1,   7],  [  7,  -1]],
[[  5,   5],  [  5,  -5]],
[[  0,   6],  [  8,   0]],
[[  2,   4],  [  8,  -4]],
[[  2,   4],  [ -8,   4]],
[[  3,   6],  [ -8,   4]],
[[  0,   5],  [ 10,   0]],
[[  3,   4],  [ -8,   6]],
[[  0,   8],  [  8,   0]],
[[  1,   8],  [  8,  -1]],
[[  0,   7],  [  9,   0]],
[[  2,   6],  [ -9,   3]],
[[  1,   5],  [-10,   2]],
[[  0,   8],  [  9,   0]],
[[  2,   5],  [-10,   4]],
[[  2,   4],  [-10,   5]],
[[  0,   9],  [  9,   0]],
]


#quadruplets:

Q = [
[[  0,   5],  [  5,   0],  [  3,  -3]],
[[  0,   5],  [  5,   0],  [  4,  -4]],
[[  0,   5],  [  5,   0],  [  5,  -5]],
[[  3,   4],  [  4,  -3],  [ -1,   7]],
[[  3,   3],  [ -4,   4],  [  7,   1]],
[[  4,  -4],  [  4,   4],  [  0, -5]],
[[  4,  -4],  [  4,   4],  [  0,   6]],
[[  4,  -4],  [  4,   4],  [  0,   7]],
[[  4,  -4],  [  4,   4],  [  0,   8]],
[[  4,  -4],  [  4,   4],  [  0,   9]],
[[  4,  -4],  [  4,   4],  [  0,  10]],
[[  0,   6],  [  6,   0],  [  3,  -3]],
[[  0,   6],  [  6,   0],  [  4,  -4]],
[[  0,   7],  [  7,   0],  [  3,  -3]],
[[  1,   7],  [  7,  -1],  [ -3,  4]], 
[[  0,   6],  [  8,   0],  [ -3,  4]], 
[[  2,   4],  [ -8,   4],  [ -4, -3]], 
[[  5,   5],  [  5,  -5],  [  0,   5]], 
[[  5,   5],  [  5,  -5],  [  0,   6]], 
[[  5,   5],  [  5,  -5],  [  0,   7]], 
[[  5,   5],  [  5,  -5],  [  0,   8]], 
[[  5,   5],  [  5,  -5],  [  0,   9]], 
[[  5,   5],  [  5,  -5],  [  0,  10]], 
[[  2,   4],  [  8,  -4],  [  0,   5]], 
[[  2,   4],  [  8,  -4],  [  0,   6]], 
[[  2,   4],  [  8,  -4],  [  0,   7]], 
[[  2,   4],  [  8,  -4],  [  0,   8]], 
[[  2,   4],  [  8,  -4],  [  0,   9]], 
[[  2,   4],  [  8,  -4],  [  0,  10]], 
[[  3,   6],  [ -8,   4],  [ -4,  -2]], 
[[  0,   5],  [ 10,   0],  [ -2,   4]], 
[[  3,   4],  [ -8,   6],  [ -4,  -2]], 
]

#-- additional data:

T1= [
[[  1,   4],  [  4,  -1]],
[[  1,  -4],  [  4,   1]],
[[ -1,   4],  [  4,   1]],
[[ -4,  -1],  [  1,  -4]],
[[ -4,   1],  [  1,   4]],
[[ -4,  -1],  [ -1,   4]],
[[ -1,  -4],  [  4,  -1]],
[[ -4,   1],  [ -1,  -4]],
[[ -3,  -3],  [ -3,   3]],
[[ -3,  -3],  [  3,  -3]],
[[ -3,   3],  [  3,   3]],
[[  3,  -3],  [  3,   3]],
[[  1,   4],  [  4,  -1]],
[[  3,   3],  [  3,  -3]],
[[  2,   4],  [  4,  -2]],
[[  4,   4],  [  4,  -4]],
[[  0,   6],  [  6,   0]],
[[  6,   0],  [  0,   7]],
[[  1,   4],  [ -8,   2]],
[[  1,   4],  [  8,  -2]],
]

Q1 = [
[[ -1,  4],  [  4,   1],  [  5,  -3]],
[[  3, -3],  [  3,   3],  [  0,   5]],
[[ -3,  3],  [  3,   3],  [  6,   0]],
[[  3, -3],  [  3,   3],  [  0,   7]],
[[ -3,  3],  [  3,   3],  [  8,   0]],
[[  3, -3],  [  3,   3],  [  0,   9]],
[[ -3,  3],  [  3,   3],  [ 10,   0]],
[[ -2,  4],  [  4,   2],  [  6,  -2]],
[[ -2,  4],  [  4,   2],  [ -9,   3]],
]



#==========================================================


for C in T1:
    name = "{:+03d}{:+03d}{:+03d}{:+03d}".format(
        C[0][0], C[0][1], C[1][0], C[1][1])
    name = outdir + '/' + run + name
    prefix = run + "s"
    print(C, name)
    corr = nlw.CorrelatorsKlo(name, prefix, C, N, No, seeds, istart = istart,  iend = iend)


for C in Q1:
    name = "{:+03d}{:+03d}{:+03d}{:+03d}{:+03d}{:+03d}".format(
        C[0][0], C[0][1], C[1][0], C[1][1], C[2][0], C[2][1] )
    name = outdir + '/' + run + name
    prefix = run + "s"
    print(C, name)
    corr = nlw.CorrelatorsKlo(name, prefix, C, N, No, seeds, istart = istart,  iend = iend)

