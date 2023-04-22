
import fiboPost



#-- Post.6 --

#runs = ['s3m10', 's3m50']

#runs = ['s4m10', 's4m90', 's5m10', 's5m90',
#        's3v14',  's3v38',  's3v58',  's3v34',
#        's4v14',  's4v38',  's4v58',  's4v34']

#runs = ['s3v58',  's3v34',  's4v58',  's4v34']


runs = ['s3v58']

print(runs)

"""
for f in runs:
   fiboPost.modeCorr(f, [-2, -1], convert=False, fbaseout='Post/'+f+'_21');
   fiboPost.modeCorr(f, [-3, -2, -2], convert=False, fbaseout='Post/'+f+'_322');
   fiboPost.modeCorr(f, [-4, -3, -1], convert=False, fbaseout='Post/'+f+'_431');
   fiboPost.modeCorr(f, [-4, -3, -3, -2], convert=False, fbaseout='Post/'+f+'_4332');
   fiboPost.modeCorr(f, [-5, -4, -2, -2], convert=False, fbaseout='Post/'+f+'_5422');
   fiboPost.modeCorr(f, [-5, -4, -4, -1], convert=False, fbaseout='Post/'+f+'_5441');
   fiboPost.modeCorr(f, [-6, -5, -3, -1], convert=False, fbaseout='Post/'+f+'_6531');

for f in runs:
   fiboPost.modeCorr(f, [-4, -4, -3, -3, -3], convert=False, fbaseout='Post/'+f+'_44333');
   fiboPost.modeCorr(f, [-5, -4, -4, -3, -2], convert=False, fbaseout='Post/'+f+'_54432');
   fiboPost.modeCorr(f, [-6, -3, -3, -3, -3], convert=False, fbaseout='Post/'+f+'_63333');
   fiboPost.modeCorr(f, [-6, -5, -3, -3, -2], convert=False, fbaseout='Post/'+f+'_65332');
   fiboPost.modeCorr(f, [-6, -5, -5, -2, -2], convert=False, fbaseout='Post/'+f+'_65522');
   fiboPost.modeCorr(f, [-6, -5, -5, -4, -1], convert=False, fbaseout='Post/'+f+'_65541');
   fiboPost.modeCorr(f, [-7, -6, -4, -2, -2], convert=False, fbaseout='Post/'+f+'_76422');
   fiboPost.modeCorr(f, [-7, -6, -4, -4, -1], convert=False, fbaseout='Post/'+f+'_76441');
   fiboPost.modeCorr(f, [-7, -6, -6, -3, -1], convert=False, fbaseout='Post/'+f+'_76631');
   fiboPost.modeCorr(f, [-8, -7, -5, -3, -1], convert=False, fbaseout='Post/'+f+'_87531');

for f in runs:
   fiboPost.modeCorr(f, [-5, -4, -4, -4, -3, -3], convert=False, fbaseout='Post/'+f+'_544433');
   fiboPost.modeCorr(f, [-5, -5, -4, -4, -4, -2], convert=False, fbaseout='Post/'+f+'_554442');
   fiboPost.modeCorr(f, [-6, -5, -4, -3, -3, -3], convert=False, fbaseout='Post/'+f+'_654333');
   fiboPost.modeCorr(f, [-6, -5, -5, -4, -3, -2], convert=False, fbaseout='Post/'+f+'_655432');
   fiboPost.modeCorr(f, [-6, -6, -5, -5, -5, -1], convert=False, fbaseout='Post/'+f+'_665551');

for f in runs: 
   fiboPost.modeCorr(f, [-7, -6, -4, -4, -3, -2], convert=False, fbaseout='Post/'+f+'_764432');
   fiboPost.modeCorr(f, [-7, -6, -6, -3, -3, -2], convert=False, fbaseout='Post/'+f+'_766332');
   fiboPost.modeCorr(f, [-7, -6, -6, -5, -2, -2], convert=False, fbaseout='Post/'+f+'_766522');
   fiboPost.modeCorr(f, [-7, -6, -6, -5, -4, -1], convert=False, fbaseout='Post/'+f+'_766541');
   fiboPost.modeCorr(f, [-8, -7, -5, -3, -3, -2], convert=False, fbaseout='Post/'+f+'_875332');

for f in runs:
   fiboPost.modeCorr(f, [-8, -7, -5, -5, -2, -2], convert=False, fbaseout='Post/'+f+'_875522');
   fiboPost.modeCorr(f, [-8, -7, -7, -4, -2, -2], convert=False, fbaseout='Post/'+f+'_877422');
   fiboPost.modeCorr(f, [-8, -7, -5, -5, -4, -1], convert=False, fbaseout='Post/'+f+'_875541');
   fiboPost.modeCorr(f, [-8, -7, -7, -4, -4, -1], convert=False, fbaseout='Post/'+f+'_877441');
   fiboPost.modeCorr(f, [-8, -7, -7, -6, -3, -1], convert=False, fbaseout='Post/'+f+'_877631');
   fiboPost.modeCorr(f, [-5, -5, -4, -4, -4, -4, -3], convert=False, fbaseout='Post/'+f+'_5544443');
   fiboPost.modeCorr(f, [-5, -5, -5, -4, -4, -4, -4, -4], convert=False, fbaseout='Post/'+f+'_55544444');
"""

for f in runs:
   fiboPost.Moments(f, nmom=12, convert=False, fbaseout='Post/'+f);
   






