
import fiboPost


outdir="Post05"



#runs = ['s4v14',  's4v38',  's4v58',  's4v34']
#runs = ['s5m10', 's5m90']


f='s5m90'
fiboPost.modeCorr(f, [10,9,7,5,3,1,], [0,],   fbaseout = outdir + '/' + f + '_A97531-0')


"""

#01
fiboPost.modeCorr(f, [5,4,4,4,3,3,], [0,],   fbaseout = outdir + '/' + f + '_544433-0')
fiboPost.modeCorr(f, [5,5,4,4,4,2,], [0,],   fbaseout = outdir + '/' + f + '_554442-0')
fiboPost.modeCorr(f, [6,5,4,3,3,3,], [0,],   fbaseout = outdir + '/' + f + '_654333-0')
fiboPost.modeCorr(f, [6,5,5,4,3,2,], [0,],   fbaseout = outdir + '/' + f + '_655432-0')
fiboPost.modeCorr(f, [6,6,5,5,5,1,], [0,],   fbaseout = outdir + '/' + f + '_665551-0')
fiboPost.modeCorr(f, [7,4,4,4,4,2,], [0,],   fbaseout = outdir + '/' + f + '_744442-0')
fiboPost.modeCorr(f, [7,6,4,4,3,2,], [0,],   fbaseout = outdir + '/' + f + '_764432-0')
fiboPost.modeCorr(f, [7,6,6,3,3,2,], [0,],   fbaseout = outdir + '/' + f + '_766332-0')
fiboPost.modeCorr(f, [7,6,6,5,2,2,], [0,],   fbaseout = outdir + '/' + f + '_766522-0')
fiboPost.modeCorr(f, [7,6,6,5,4,1,], [0,],   fbaseout = outdir + '/' + f + '_766541-0')

#02

fiboPost.modeCorr(f, [8,5,5,5,5,1,], [0,],   fbaseout = outdir + '/' + f + '_855551-0')
fiboPost.modeCorr(f, [8,7,3,3,3,3,], [0,],   fbaseout = outdir + '/' + f + '_873333-0')
fiboPost.modeCorr(f, [8,7,5,3,3,2,], [0,],   fbaseout = outdir + '/' + f + '_875332-0')
fiboPost.modeCorr(f, [8,7,5,5,2,2,], [0,],   fbaseout = outdir + '/' + f + '_875522-0')
fiboPost.modeCorr(f, [8,7,5,5,4,1,], [0,],   fbaseout = outdir + '/' + f + '_875541-0')
fiboPost.modeCorr(f, [8,7,7,4,2,2,], [0,],   fbaseout = outdir + '/' + f + '_877422-0')
fiboPost.modeCorr(f, [8,7,7,4,4,1,], [0,],   fbaseout = outdir + '/' + f + '_877441-0')
fiboPost.modeCorr(f, [8,7,7,6,3,1,], [0,],   fbaseout = outdir + '/' + f + '_877631-0')

#03
fiboPost.modeCorr(f, [9,8,6,4,2,2,], [0,],   fbaseout = outdir + '/' + f + '_986422-0')
fiboPost.modeCorr(f, [9,8,6,4,4,1,], [0,],   fbaseout = outdir + '/' + f + '_986441-0')
fiboPost.modeCorr(f, [9,8,6,6,3,1,], [0,],   fbaseout = outdir + '/' + f + '_986631-0')
fiboPost.modeCorr(f, [9,8,8,5,3,1,], [0,],   fbaseout = outdir + '/' + f + '_988531-0')
fiboPost.modeCorr(f, [10,9,7,5,3,1,], [0,],   fbaseout = outdir + '/' + f + '_A97531-0')

#04
fiboPost.modeCorr(f, [4,3,3,3,3,], [5,0,],   fbaseout = outdir + '/' + f + '_43333-50')
fiboPost.modeCorr(f, [4,4,4,3,2,], [6,0,],   fbaseout = outdir + '/' + f + '_44432-60')
fiboPost.modeCorr(f, [5,2,2,2,2,], [1,0,],   fbaseout = outdir + '/' + f + '_52222-10')
fiboPost.modeCorr(f, [5,3,3,3,3,], [7,0,],   fbaseout = outdir + '/' + f + '_53333-70')
fiboPost.modeCorr(f, [5,5,3,3,2,], [7,0,],   fbaseout = outdir + '/' + f + '_55332-70')
fiboPost.modeCorr(f, [5,5,5,2,2,], [7,0,],   fbaseout = outdir + '/' + f + '_55522-70')
fiboPost.modeCorr(f, [5,5,5,4,1,], [7,0,],   fbaseout = outdir + '/' + f + '_55541-70')
fiboPost.modeCorr(f, [6,3,3,3,1,], [2,0,],   fbaseout = outdir + '/' + f + '_63331-20')
fiboPost.modeCorr(f, [6,3,3,3,2,], [4,0,],   fbaseout = outdir + '/' + f + '_63332-40')
fiboPost.modeCorr(f, [6,5,1,1,1,], [0,0,],   fbaseout = outdir + '/' + f + '_65111-00')

#05
fiboPost.modeCorr(f, [6,6,4,2,2,], [8,0,],   fbaseout = outdir + '/' + f + '_66422-80')
fiboPost.modeCorr(f, [6,6,4,4,1,], [8,0,],   fbaseout = outdir + '/' + f + '_66441-80')
fiboPost.modeCorr(f, [6,6,6,3,1,], [8,0,],   fbaseout = outdir + '/' + f + '_66631-80')
fiboPost.modeCorr(f, [7,4,4,2,2,], [5,0,],   fbaseout = outdir + '/' + f + '_74422-50')
fiboPost.modeCorr(f, [7,4,4,4,1,], [5,0,],   fbaseout = outdir + '/' + f + '_74441-50')
fiboPost.modeCorr(f, [7,6,2,2,2,], [3,0,],   fbaseout = outdir + '/' + f + '_76222-30')
fiboPost.modeCorr(f, [7,6,6,1,1,], [2,0,],   fbaseout = outdir + '/' + f + '_76611-20')
fiboPost.modeCorr(f, [7,7,5,3,1,], [9,0,],   fbaseout = outdir + '/' + f + '_77531-90')
fiboPost.modeCorr(f, [8,7,3,3,1,], [4,0,],   fbaseout = outdir + '/' + f + '_87331-40')
fiboPost.modeCorr(f, [8,7,5,1,1,], [2,0,],   fbaseout = outdir + '/' + f + '_87511-20')
fiboPost.modeCorr(f, [8,5,5,3,1,], [6,0,],   fbaseout = outdir + '/' + f + '_85531-60')


#06
fiboPost.modeCorr(f, [3,1,1,1,], [5,0,0,],   fbaseout = outdir + '/' + f + '_3111-500')
fiboPost.modeCorr(f, [3,3,3,2,], [7,6,0,],   fbaseout = outdir + '/' + f + '_3332-760')
fiboPost.modeCorr(f, [4,2,2,2,], [6,3,0,],   fbaseout = outdir + '/' + f + '_4222-630')
fiboPost.modeCorr(f, [4,4,2,2,], [8,7,0,],   fbaseout = outdir + '/' + f + '_4422-870')
fiboPost.modeCorr(f, [4,4,4,1,], [8,7,0,],   fbaseout = outdir + '/' + f + '_4441-870')
fiboPost.modeCorr(f, [5,3,3,1,], [7,4,0,],   fbaseout = outdir + '/' + f + '_5331-740')
fiboPost.modeCorr(f, [5,5,1,1,], [7,2,0,],   fbaseout = outdir + '/' + f + '_5511-720')
fiboPost.modeCorr(f, [5,5,3,1,], [9,8,0,],   fbaseout = outdir + '/' + f + '_5531-980')

#07
fiboPost.modeCorr(f, [6,3,1,1,], [4,2,0,],   fbaseout = outdir + '/' + f + '_6311-420')
fiboPost.modeCorr(f, [7,2,2,2,], [5,5,0,],   fbaseout = outdir + '/' + f + '_7222-550')
fiboPost.modeCorr(f, [8,3,3,1,], [6,6,0,],   fbaseout = outdir + '/' + f + '_8331-660')
fiboPost.modeCorr(f, [8,7,1,1,], [4,4,0,],   fbaseout = outdir + '/' + f + '_8711-440')
fiboPost.modeCorr(f, [9,6,2,1,], [7,7,0,],   fbaseout = outdir + '/' + f + '_9621-770')

#08
fiboPost.modeCorr(f, [1,1,1,], [5,2,2,0,],   fbaseout = outdir + '/' + f + '_111-5220')
fiboPost.modeCorr(f, [2,2,2,], [7,6,6,0,],   fbaseout = outdir + '/' + f + '_222-7660')
fiboPost.modeCorr(f, [2,2,2,], [8,7,5,0,],   fbaseout = outdir + '/' + f + '_222-8750')
fiboPost.modeCorr(f, [3,1,1,], [7,6,2,0,],   fbaseout = outdir + '/' + f + '_311-7620')
fiboPost.modeCorr(f, [3,3,1,], [8,7,7,0,],   fbaseout = outdir + '/' + f + '_331-8770')
fiboPost.modeCorr(f, [3,3,1,], [9,8,6,0,],   fbaseout = outdir + '/' + f + '_331-9860')
fiboPost.modeCorr(f, [5,1,1,], [7,4,4,0,],   fbaseout = outdir + '/' + f + '_511-7440')
fiboPost.modeCorr(f, [7,1,1,], [5,5,5,0,],   fbaseout = outdir + '/' + f + '_711-5550')
fiboPost.modeCorr(f, [8,1,1,], [6,6,4,0,],   fbaseout = outdir + '/' + f + '_811-6640')

#09
fiboPost.modeCorr(f, [1,1,], [7,6,6,5,0,],   fbaseout = outdir + '/' + f + '_11-76650')
fiboPost.modeCorr(f, [1,1,], [8,7,5,5,0,],   fbaseout = outdir + '/' + f + '_11-87550')
fiboPost.modeCorr(f, [1,1,], [8,7,7,4,0,],   fbaseout = outdir + '/' + f + '_11-87740')
fiboPost.modeCorr(f, [1,1,], [9,8,6,4,0,],   fbaseout = outdir + '/' + f + '_11-98640')

"""





