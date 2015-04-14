import numpy as np
import scipy
import scipy.stats

"""
Ordering used for the generation of result tables.
"""
ordering = ['20NG-m5','AP','DOE','FR','SJMN','WSJ','ZF','classic','fbis','hitech','reviews','sports','la12','new3','tr31','tr41','ohscal','re0','re1','k1a','k1b','wap']


"""
Result files grouped by hash family (minwise or rhyperplanes), estimation method(mode 1 or mode 2) and signature block size (500 or 1000).
"""
results = {'minwise_mod1_500': ["20ng-results_minwise-results_mod1-r1b500", "20ng-results_minwise-results_mod1-r2b500", "20ng-results_minwise-results_mod1-r3b500",   "tipster-results_minwise-results_mod1-r1b500", "tipster-results_minwise-results_mod1-r2b500", "tipster-results_minwise-results_mod1-r3b500",  "cluto-results_minwise-results_mod1-r1b500", "cluto-results_minwise-results_mod1-r2b500","cluto-results_minwise-results_mod1-r3b500"], 

    'minwise_mod1_1000': ["20ng-results_minwise-results_mod1-r1b1000", "20ng-results_minwise-results_mod1-r2b1000", "20ng-results_minwise-results_mod1-r3b1000",   "tipster-results_minwise-results_mod1-r1b1000", "tipster-results_minwise-results_mod1-r2b1000", "tipster-results_minwise-results_mod1-r3b1000",  "cluto-results_minwise-results_mod1-r1b1000", "cluto-results_minwise-results_mod1-r2b1000","cluto-results_minwise-results_mod1-r3b1000"], 
    
    'minwise_mod2_500': ["20ng-results_minwise-results_mod2-r1b500", "20ng-results_minwise-results_mod2-r2b500", "20ng-results_minwise-results_mod2-r3b500",   "tipster-results_minwise-results_mod2-r1b500", "tipster-results_minwise-results_mod2-r2b500", "tipster-results_minwise-results_mod2-r3b500",  "cluto-results_minwise-results_mod2-r1b500", "cluto-results_minwise-results_mod2-r2b500","cluto-results_minwise-results_mod2-r3b500"], 

        'minwise_mod2_1000': ["20ng-results_minwise-results_mod2-r1b1000", "20ng-results_minwise-results_mod2-r2b1000", "20ng-results_minwise-results_mod2-r3b1000",   "tipster-results_minwise-results_mod2-r1b1000", "tipster-results_minwise-results_mod2-r2b1000", "tipster-results_minwise-results_mod2-r3b1000",  "cluto-results_minwise-results_mod2-r1b1000", "cluto-results_minwise-results_mod2-r2b1000","cluto-results_minwise-results_mod2-r3b1000"],
    
    'rhyperplanes_mod1_500': ["20ng-results_rhyperplanes-results_mod1-r1b500", "20ng-results_rhyperplanes-results_mod1-r2b500", "20ng-results_rhyperplanes-results_mod1-r3b500",   "tipster-results_rhyperplanes-results_mod1-r1b500", "tipster-results_rhyperplanes-results_mod1-r2b500", "tipster-results_rhyperplanes-results_mod1-r3b500",  "cluto-results_rhyperplanes-results_mod1-r1b500", "cluto-results_rhyperplanes-results_mod1-r2b500","cluto-results_rhyperplanes-results_mod1-r3b500"], 

        'rhyperplanes_mod1_1000': ["20ng-results_rhyperplanes-results_mod1-r1b1000", "20ng-results_rhyperplanes-results_mod1-r2b1000", "20ng-results_rhyperplanes-results_mod1-r3b1000",   "tipster-results_rhyperplanes-results_mod1-r1b1000", "tipster-results_rhyperplanes-results_mod1-r2b1000", "tipster-results_rhyperplanes-results_mod1-r3b1000",  "cluto-results_rhyperplanes-results_mod1-r1b1000", "cluto-results_rhyperplanes-results_mod1-r2b1000","cluto-results_rhyperplanes-results_mod1-r3b1000"], 
    
        'rhyperplanes_mod2_500': ["20ng-results_rhyperplanes-results_mod2-r1b500", "20ng-results_rhyperplanes-results_mod2-r2b500", "20ng-results_rhyperplanes-results_mod2-r3b500",   "tipster-results_rhyperplanes-results_mod2-r1b500", "tipster-results_rhyperplanes-results_mod2-r2b500", "tipster-results_rhyperplanes-results_mod2-r3b500",  "cluto-results_rhyperplanes-results_mod2-r1b500", "cluto-results_rhyperplanes-results_mod2-r2b500","cluto-results_rhyperplanes-results_mod2-r3b500"], 

            'rhyperplanes_mod2_1000': ["20ng-results_rhyperplanes-results_mod2-r1b1000", "20ng-results_rhyperplanes-results_mod2-r2b1000", "20ng-results_rhyperplanes-results_mod2-r3b1000",   "tipster-results_rhyperplanes-results_mod2-r1b1000", "tipster-results_rhyperplanes-results_mod2-r2b1000", "tipster-results_rhyperplanes-results_mod2-r3b1000",  "cluto-results_rhyperplanes-results_mod2-r1b1000", "cluto-results_rhyperplanes-results_mod2-r2b1000","cluto-results_rhyperplanes-results_mod2-r3b1000"]}
            
            
"""
Performance values attained by the clustering algorithm over the exact similarities.
Jaccard and Cosine.
"""            
jacc = {'20NG-m5':{'E':0.365,'P':0.8283},
'AP': {'E':0.3066, 'P':0.5853 }, 
'DOE': {'E':0.2153, 'P':0.7825},  
'FR': {'E':0.254, 'P':0.7701 }, 
'SJMN': {'E':0.2264, 'P':0.773 }, 
'WSJ': {'E':0.3294, 'P':0.5712 }, 
'ZF': {'E':0.3576, 'P':0.6443}, 
'classic': {'E':0.2056, 'P':0.9307},  
'fbis': {'E':0.3921, 'P':0.6217 }, 
'hitech': {'E':0.7115, 'P':0.4846 }, 
'reviews': {'E':0.5678, 'P':0.6079},  
'sports': {'E':0.2547, 'P':0.8164},  
'la12': {'E':0.5869, 'P':0.6327 }, 
'new3': {'E':0.4453, 'P':0.5033 }, 
'tr31': {'E':0.316, 'P':0.7374 }, 
'tr41': {'E':0.2406, 'P':0.795 }, 
'ohscal': {'E':0.7476, 'P':0.3942},  
're0': {'E':0.4479, 'P':0.5855 }, 
're1': {'E':0.3386, 'P':0.6421},  
'k1a': {'E':0.3728, 'P':0.641},  
'k1b': {'E':0.1724, 'P':0.8872},  
'wap': {'E':0.3568, 'P':0.6628}}
            
            
cosine = {"20NG-m5":{"E":0.3574,"P":0.8312},
"AP":{"E":0.3042,"P":0.5837},
"DOE":{"E":0.2261,"P":0.7596},
"FR":{"E":0.2489,"P":0.7714},
"SJMN":{"E":0.2263,"P":0.7719},
"WSJ":{"E":0.3327,"P":0.5723},
"ZF":{"E":0.3601,"P":0.6466},
"classic":{"E":0.2297,"P":0.9189},
"fbis":{"E":0.3802,"P":0.6337},
"hitech":{"E":0.7084,"P":0.508},
"reviews":{"E":0.5305,"P":0.6822},
"sports":{"E":0.2791,"P":0.799},
"la12":{"E":0.5946,"P":0.6027},
"new3":{"E":0.4324,"P":0.5223},
"tr31":{"E":0.327,"P":0.7434},
"tr41":{"E":0.2302,"P":0.8018},
"ohscal":{"E":0.7404,"P":0.3859},
"re0":{"E":0.4711,"P":0.5682},
"re1":{"E":0.307,"P":0.6874},
"k1a":{"E":0.3666,"P":0.6504},
"k1b":{"E":0.1745,"P":0.891},
"wap":{"E":0.3386,"P":0.6865}}
            
def _conf_interval(xbar, s, N):
    std_error = s/np.sqrt(N)
    interv = scipy.stats.t.interval(0.95, N-1)[1]
    lower = xbar - interv*std_error
    upper = xbar + interv*std_error
    return lower, upper          

            
def conf_interval(X):#sample data given
    """
    Computes the confidence interval with 95% confidence of the
    sample X, whose population deviation is Unknown.
    Note: Good tutorial in http://onlinestatbook.com/2/estimation/mean.html
    """
    xbar, s2 = np.mean(X), np.var(X,ddof=1)
    return _conf_interval(xbar, np.sqrt(s2), len(X))
            
            
def remove_extension(L):
    return L.replace("_out.dat","").replace(".dat","").replace(".mat","")

def datasets_mean_std_differences(fname):
    """
    Given a result file, reads it and then outpus a data structure
    containg the means a std for the measures appearing in the file for each dataset.
    Note: The dataset name is filtered by the function remove_extension and the resulting string
    must be contained into the ordering list.
    """
    nruns = 0
    dsets = {}
    fieldnames = []
    
    fieldslist = "E,P"

    f = open(fname)

    for l in f:
        if len(l.strip().split(";")) == 1:
            if l.startswith("RUN"):
                nruns += 1
            else:# A filename appears.
                curFile = remove_extension(l.strip())
        else: # the measures names appear.
            if l.startswith("E"):
                nfields = len(l.strip().split(";"))
                
                if len(fieldnames) == 0: # initialize the field names.
                    fieldnames = l.strip().replace(" ",'').split(";")
                    
                if not curFile in dsets:
                    dsets[curFile] = [[] for i in range(nfields)]
            else:# the measures values appear.
                fields = map(float, l.strip().split(";"))
                nfields = len(fields)

                for i in range(nfields):
                    a = dsets[curFile][i]
                    a.append(fields[i])
                
    f.close()
    #print "Number of runs:",nruns

    fieldsToReport = []
    for ff in fieldslist.replace(" ",'').split(","):
        fieldsToReport.append( fieldnames.index(ff) )

    out = {}
    for f in dsets:
        #print "%FILE:",f,"num. measures:",len(dsets[f]),"reported measures:",sys.argv[2]
        #out[f] = {'mean':[], 'std':[]}        
        out[f] = {}
        i = fieldnames.index('E')
        if 'minwise' in fname:
            REF_VALUE = jacc[f]['E']
        else:
            REF_VALUE = cosine[f]['E']
        out[f]['E'] = {}
        out[f]['E']['mean'] = np.mean(np.array(dsets[f][i]) -  REF_VALUE )
        out[f]['E']['std'] = np.std(np.array(dsets[f][i]) - REF_VALUE , ddof=1)

        i = fieldnames.index('P')
        if 'minwise' in fname:
            REF_VALUE = jacc[f]['P']
        else:
            REF_VALUE = cosine[f]['P']
        out[f]['P'] = {}
        out[f]['P']['mean'] = np.mean(REF_VALUE - np.array(dsets[f][i]) )
        out[f]['P']['std'] = np.std(REF_VALUE - np.array(dsets[f][i]) , ddof=1 )

        if out[f]['P']['mean'] < 0:
            print "[",fname,"]NEgative Purity mean difference for",f

        if out[f]['E']['mean'] < 0:
            print "[",fname,"]NEgative Entropy mean difference for",f
        
    return out


def datasets_mean_std(fname, fieldslist):
    """
    Given a result file, reads it and then outpus a data structure
    containg the means a std for the measures appearing in the file for each dataset.
    Note: The dataset name is filtered by the function remove_extension and the resulting string
    must be contained into the ordering list.
    """
    nruns = 0
    dsets = {}
    fieldnames = []

    f = open(fname)

    for l in f:
        if len(l.strip().split(";")) == 1:
            if l.startswith("RUN"):
                nruns += 1
            else:# A filename appears.
                curFile = remove_extension(l.strip())
        else: # the measures names appear.
            if l.startswith("E"):
                nfields = len(l.strip().split(";"))
                
                if len(fieldnames) == 0: # initialize the field names.
                    fieldnames = l.strip().replace(" ",'').split(";")
                    
                if not curFile in dsets:
                    dsets[curFile] = [[] for i in range(nfields)]
            else:# the measures values appear.
                fields = map(float, l.strip().split(";"))
                nfields = len(fields)

                for i in range(nfields):
                    a = dsets[curFile][i]
                    a.append(fields[i])
                
    f.close()
    #print "Number of runs:",nruns

    fieldsToReport = []
    for ff in fieldslist.replace(" ",'').split(","):
        fieldsToReport.append( fieldnames.index(ff) )

    out = {}
    for f in dsets:
        #print "%FILE:",f,"num. measures:",len(dsets[f]),"reported measures:",sys.argv[2]
        #out[f] = {'mean':[], 'std':[]}
        out[f] = {}
        for ff in fieldslist.replace(" ",'').split(","):
            i = fieldnames.index(ff)
            out[f][ff] = {}
            out[f][ff]['mean'] = np.mean(dsets[f][i])
            out[f][ff]['std'] = np.std(dsets[f][i] , ddof=1)
    return out

def _plot_entropy_purity(R, outpath=None):
    """
    Given a datastructure generated by the method 'datasets_mean_std' generates a barplot
    for the results for each dataset.
    The ordering of the datasets is given by the list of the first line of this file.
    """
    
    meanEs = []
    meanPs = []
    stdEs = []
    stdPs = []
    labels = []
    for dset in ordering: #  ensures the fullfilment of the desired order for the bars
        if dset in R:
            labels.append(dset)
            meanEs.append( R[dset]['E']['mean'] )
            stdEs.append(  R[dset]['E']['std'] )
            meanPs.append( R[dset]['P']['mean'] )
            stdPs.append(  R[dset]['P']['std'] )
        else:
            print "ERROR",dset," must be present in the ordering list."
            


    import numpy as np
    import matplotlib.pyplot as plt
    N = len(R)
    ind = np.arange(0,5*N,5)    # the x locations for the groups
    width = 1.8       # the width of the bars: can also be len(x) sequence
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, meanEs, width, color='y', yerr=stdEs, error_kw={'ecolor': '0.1'}, hatch=".",alpha=0.4,log=True)
    
    rects2 = ax.bar(ind+width, meanPs, width, color='brown', yerr=stdPs, error_kw={'ecolor': '0.1'}, hatch="//",alpha=0.7,log=True)
    
    # add some text for labels, title and axes ticks
    ax.set_ylabel('Measure value')
    #ax.set_title('Entropy and Purity value for different datasets')
    ax.set_xticks(ind+width-0.2)
    
    ax.set_xticklabels( labels , rotation="vertical")

    plt.xlim((-3.0, 5*N + 3))
    plt.ylim((-0.4, 0.6)) # for plotting raw levels

    
    ax.legend( (rects1[0], rects2[0]), ('Entropy', 'Purity') )
    plt.grid()
    if not outpath:
        plt.show()
    else:
        plt.savefig(outpath,bbox_inches='tight')
        plt.close()

def plot_conf_intervals(R, outpath=None):
    x = np.arange(0, len(R)*3, 3)
    Elower = []
    Eupper = []
    Plower = []
    Pupper = []

    labels = []
    for dset in ordering: #  ensures the fullfilment of the desired order for the bars
        if dset in R:
            labels.append(dset)
            Elow, Eup = _conf_interval(R[dset]['E']['mean'], R[dset]['E']['std'], 10)
            Plow, Pup = _conf_interval(R[dset]['P']['mean'], R[dset]['P']['std'], 10)
            Elower.append(Elow)
            Eupper.append(Eup)
            Plower.append(Plow)
            Pupper.append(Pup)
        else:
            print "ERROR",dset," must be present in the ordering list."
               
    N = len(R)
    ind = np.arange(0,5*N,5)    # the x locations for the groups
    width = 1.8       # the width of the bars: can also be len(x) sequence
    
    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)    
    ax1.fill_between(ind+width, Elower, Eupper, color='b', alpha = 0.4)
    ax1.plot(ind+width, np.zeros(N),'black')
    ax1.set_ylabel('$\Delta=\mathbf{E}_{est}-\mathbf{E}_{exact}$')
    ax1.grid()
    
    ax2.fill_between(ind+width, Plower, Pupper, color='r', alpha = 0.4)
    ax2.plot(ind+width, np.zeros(N),'black')
    ax2.set_ylabel('$\Delta=\mathbf{P}_{exact}-\mathbf{P}_{est}$')
    ax2.grid()
    ax2.set_xticks(ind+width-0.2)    
    ax2.set_xticklabels( labels , rotation="vertical")
    ax2.set_xlim(right=5*N-2)
    #ax2.legend((rects1[0], rects2[0]),('Entropy', 'Purity') )  
    if not outpath:
        plt.show()
    else:
        plt.savefig(outpath,format='eps',bbox_inches='tight')
        plt.close()



def plot_entropy_purity(fname, outpath=None):
    R = datasets_mean_std(fname, "E,P")
    _plot_entropy_purity(R, outpath)
    

########### END OF METHOD DEFINITIONS ##########

if __name__ == '2__main__':
    import sys
    import numpy as np


    if len(sys.argv) < 2:
        print "Usage:",sys.argv[0],"result-file [plot-output-file]"
        sys.exit(-1)

    outpath = None
    if len(sys.argv) == 3:
        outpath = sys.argv[2]
        
    plot_entropy_purity(sys.argv[1], outpath)
    

if __name__ == '2__main__':
    import sys
    import numpy as np
    
    # Environment variables to indicate the desired output.
    GEN_LATEX_CODE = True
    GEN_AND_SAVE_PLOTS = False

    s_texrow = '{DATA} &{E_r1:.4f} ({E_r1_std:.4f}) &{E_r2:.4f} ({E_r2_std:.4f}) &{E_r3:.4f} ({E_r3_std:.4f}) & \multicolumn{{1}}{{l}}{{}}  & \multicolumn{{1}}{{l}}{{{P_r1:.4f} ({P_r1_std:.4f}) }} & {P_r2:.4f} ({P_r2_std:.4f}) & {P_r3:.4f} ({P_r3_std:.4f}) \\\\'

    #R = datasets_mean_std(sys.argv[1], "E,P")
    for texp in results:# E.g. minwise_mod1_1000
        r1b = {}
        r2b = {}
        r3b = {}
        #print "[Processing] Task", texp
        for fname in results[texp]:
            # r1b
            if 'r1b' in fname:
                #print "[r1b] adding",fname
                r1b.update( datasets_mean_std(fname, "E,P") )
            elif 'r2b' in fname:
                #print "[r2b] adding",fname                
                r2b.update( datasets_mean_std(fname, "E,P") )
            else:
                #print "[r3b] adding",fname                
                r3b.update( datasets_mean_std(fname, "E,P") )
        # latex code printing
        if GEN_LATEX_CODE:
            print "\n%LaTeX code for task", texp
            for dset in ordering:
                if dset in r1b:
                    print s_texrow.format(DATA=dset,E_r1=r1b[dset]['E']['mean'],E_r2=r2b[dset]['E']['mean'],E_r3=r3b[dset]['E']['mean'],P_r1=r1b[dset]['P']['mean'],P_r2=r2b[dset]['P']['mean'],P_r3=r3b[dset]['P']['mean'],E_r1_std=r1b[dset]['E']['std'],E_r2_std=r2b[dset]['E']['std'],E_r3_std=r3b[dset]['E']['std'],P_r1_std=r1b[dset]['P']['std'],P_r2_std=r2b[dset]['P']['std'],P_r3_std=r3b[dset]['P']['std'])
            
        # Generation of plots 
        if GEN_AND_SAVE_PLOTS:
            _plot_entropy_purity(r1b, '/Users/jz/Desktop/plots/{0}_{1}.png'.format(texp,'r1b') )
            _plot_entropy_purity(r2b, '/Users/jz/Desktop/plots/{0}_{1}.png'.format(texp,'r2b') )
            _plot_entropy_purity(r3b, '/Users/jz/Desktop/plots/{0}_{1}.png'.format(texp,'r3b') )

if __name__ == '__main__':
    import sys
    import numpy as np
    
    # Environment variables to indicate the desired output.
    GEN_LATEX_CODE = False
    GEN_AND_SAVE_PLOTS = True

    s_texrow = '{DATA} &{E_r1:.4f} ({E_r1_std:.4f}) &{E_r2:.4f} ({E_r2_std:.4f}) &{E_r3:.4f} ({E_r3_std:.4f}) & \multicolumn{{1}}{{l}}{{}}  & \multicolumn{{1}}{{l}}{{{P_r1:.4f} ({P_r1_std:.4f}) }} & {P_r2:.4f} ({P_r2_std:.4f}) & {P_r3:.4f} ({P_r3_std:.4f}) \\\\'

    #R = datasets_mean_std(sys.argv[1], "E,P")
    for texp in results:# E.g. minwise_mod1_1000
        r1b = {}
        r2b = {}
        r3b = {}
        #print "[Processing] Task", texp
        for fname in results[texp]:
            # r1b
            if 'r1b' in fname:
                #print "[r1b] adding",fname
                r1b.update( datasets_mean_std_differences("/Users/juan/Downloads/experimental_results/"+fname) )
            elif 'r2b' in fname:
                #print "[r2b] adding",fname                
                r2b.update( datasets_mean_std_differences("/Users/juan/Downloads/experimental_results/"+fname) )
            else:
                #print "[r3b] adding",fname                
                r3b.update( datasets_mean_std_differences("/Users/juan/Downloads/experimental_results/"+fname) )

        # latex code printing
        if GEN_LATEX_CODE:
            print "\n%LaTeX code for task", texp
            for dset in ordering:
                if dset in r1b:
                    print s_texrow.format(DATA=dset,E_r1=r1b[dset]['E']['mean'],E_r2=r2b[dset]['E']['mean'],E_r3=r3b[dset]['E']['mean'],P_r1=r1b[dset]['P']['mean'],P_r2=r2b[dset]['P']['mean'],P_r3=r3b[dset]['P']['mean'],E_r1_std=r1b[dset]['E']['std'],E_r2_std=r2b[dset]['E']['std'],E_r3_std=r3b[dset]['E']['std'],P_r1_std=r1b[dset]['P']['std'],P_r2_std=r2b[dset]['P']['std'],P_r3_std=r3b[dset]['P']['std'])
        # Generation of plots 
        if GEN_AND_SAVE_PLOTS:
            plot_conf_intervals(r1b,'/Users/juan/Desktop/plots/diff_{0}_{1}.eps'.format(texp,'r1b') )
            plot_conf_intervals(r2b,'/Users/juan/Desktop/plots/diff_{0}_{1}.eps'.format(texp,'r2b') )
            plot_conf_intervals(r3b,'/Users/juan/Desktop/plots/diff_{0}_{1}.eps'.format(texp,'r3b') )
        
           