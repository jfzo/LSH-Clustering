import sys
import numpy as np

nruns = 0
dsets = {}
fieldnames = []

if len(sys.argv) < 3:
    print "Usage:",sys.argv[0],"result-file","field-list-separated-by-comma"
    sys.exit(-1)
    

f = open(sys.argv[1])

for l in f:
    if len(l.strip().split(";")) == 1:
        if l.startswith("RUN"):
            nruns += 1
        else:# A filename appears.
            curFile = l.strip()
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
print "Number of runs:",nruns

fieldsToReport = []
for ff in sys.argv[2].replace(" ",'').split(","):
    fieldsToReport.append( fieldnames.index(ff) )

for f in dsets:
    print "%FILE:",f,"num. measures:",len(dsets[f]),"reported measures:",sys.argv[2]
    for i in fieldsToReport:
        print "{0:.4f}/{1:.4f} ".format(np.mean(dsets[f][i]), np.std(dsets[f][i]) ),
    print ""
