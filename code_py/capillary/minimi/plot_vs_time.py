#import importlib ; import plot_vs_time
#importlib.reload(plot_vs_time) ; plot_vs_time.All()

def All():
    Plot(ampname="a1e3", varname = "eta", fnameout='plot/a1e3eta.pdf', every = 1)
    Plot(ampname="a1e3", varname = "H",   fnameout='plot/a1e3ham.pdf', every = 1)
    Plot(ampname="a1e3", varname = "NL",  fnameout='plot/a1e3nl.pdf',  every = 1)

    Plot(ampname="a2e3", varname = "eta", fnameout='plot/a2e3eta.pdf', every = 1)
    Plot(ampname="a2e3", varname = "H",   fnameout='plot/a2e3ham.pdf', every = 1)
    Plot(ampname="a2e3", varname = "NL",  fnameout='plot/a2e3nl.pdf',  every = 1)

    
#=================================================================

def Plot(ampname="a1e3", varname = "eta", fnameout='test.pdf', every = 1):
    import matplotlib.pyplot as plt
    import numpy as np

    prefix = "runsA/dat/sub_"

    modenames = (("m05", "m07"), ("m09", "m13"), ("m19", "m21"))
    
        
    fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(9,10))

    setcanvas(ax, varname)

    marker_etamin = dict(linestyle='-', linewidth=0.5, color='b', markersize=6, mfc="none", mec="y")
    marker_etamax = dict(linestyle='-', linewidth=0.5, color='r', markersize=6, mfc="none", mec="y")
    marker_etarms = dict(linestyle='-', linewidth=0.5, color='g', markersize=6, mfc="none", mec="y")

    if varname == "eta":

        for i in (0,1,2):
            for j in (0,1):
                try:
                    runname = prefix + modenames[i][j] + ampname + ".dat"
                    print(runname)
                    plot_data(ax[i,j], runname, 'eta_min', '-', marker_etamin, every=every)
                    plot_data(ax[i,j], runname, 'eta_max', '-', marker_etamax, every=every)
                    plot_data(ax[i,j], runname, 'eta_rms', '-', marker_etarms, every=every)
                except:
                    pass

    else:
                
        for i in (0,1,2):
            for j in (0,1):
                try:
                    runname = prefix + modenames[i][j] + ampname + ".dat"
                    print(runname)
                    plot_data(ax[i,j], runname, varname, '-', marker_etamin, every=every)
                except:
                    pass
   
            
    #add_legends(ax, loc)

    plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    #plt.show()
    plt.close()


#=================================================================

def setcanvas(ax, varname):
    import matplotlib.pyplot as plt

    modes = (("5x5", "7x7"), ("9x9", "13x13"),  ("19x19", "21x21"))

    if varname == 'eta':
        ylabel = "$\eta$"
        y1=-1
        y2= 1
    elif varname == "NL":
        ylabel = "$H_{NL} / H$"
        y1 = -0.04
        y2 = 0.01
    elif varname == "H":
        ylabel = "$H$"
        y1 = 0
        y2 = 2
    else:
        ylabel = varname
    
    for i in (0,1,2):
        for j in (0,1):

            ax[i,j].set_title(modes[i][j])
            ax[i,j].set_xlabel('$t$')
            ax[i,j].set_ylabel(ylabel)
            ax[i,j].set_ylim(y1,y2)
            ax[i,j].grid()
            ax[i,j].set_xlim(0,1e6)
            ax[i,j].ticklabel_format(axis='x', style='sci', scilimits =(-3,0)) 
            if i == 2 and varname == "NL":
                ax[i,j].set_ylim(-0.12,0.02)
               
            
            
#=================================================================
   
def plot_data(iax, fname, var, marker, marker_style, every=1):
    import numpy as np
    import nlwTools
    
    v = nlwTools.data_index()
        
    dat = np.loadtxt(fname, comments="#")   
    ind = np.arange(0,len(dat), every)
    dat = dat[ind,:]
 
    t = dat[:,v["time"]]
        
    if (var == "eta_min" or var == "eta_max" or var == "eta_rms"):
        f = dat[:,v[var]]
    if (var == "H"):
        f = dat[:,v["E_pot"]] + dat[:,v["E_kin"]] + dat[:,v["E_nl"]]
    if (var == "NL"):
        h = dat[:,v["E_pot"]] + dat[:,v["E_kin"]] + dat[:,v["E_nl"]]
        f = dat[:,v["E_nl"]]/h

    #print(dat.shape, t.shape, f.shape)
    iax.plot(t, f, marker, **marker_style)
   
#=================================================================



#=================================================================

    
#=================================================================

def add_legends(ax, loc):
    import matplotlib.pyplot as plt
    
    for i in (0,1,2,3):
        for j in (0,1):
 
            ax[i,j].legend(frameon=False, loc=loc)

   
#=================================================================


#=================================================================

#=================================================================

