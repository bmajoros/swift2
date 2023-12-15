#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import math
import ProgramName
import statsmodels
from statsmodels import stats
from statsmodels.stats import multitest

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <swift2.out> <#nulls> <#alts> <alpha> <out:q-values.txt>\n")
(infile,numNulls,numAlts,alpha,outFile)=sys.argv[1:]
numNulls=int(numNulls)
numAlts=int(numAlts)
alpha=float(alpha)

# Read p-values from file
pvalues=[]
IN=open(infile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=7): continue
    (id,theta,Palt,CI_left,CI_right,median_p,area_p)=fields
    pvalues.append(float(median_p))
    #pvalues.append(float(area_p))
IN.close()

# Do FDR correction to produce q-values
#out=statsmodels.stats.multitest.fdrcorrection(pvalues,alpha=0.05,
#                                               method='indep',
#                                               is_sorted=False)
out=statsmodels.stats.multitest.multipletests(pvalues, alpha=0.5,
                                              method='fdr_bh',
                                              #method="bonferroni",
                                              maxiter=1,
                                              is_sorted=False,
                                              returnsorted=False)
qValues=out[1]
numValues=len(qValues)
if(numValues!=numNulls+numAlts):
    raise Exception("#values does not equal #nulls + #alts")
with open(outFile,"wt") as OUT:
    for x in qValues: print(x,file=OUT)

# Compute FDR and power
#qValues=[float(x) for x in q]
qNulls=qValues[:numNulls]
qAlts=qValues[numNulls:]
FP=0; TP=0
for q in qNulls:
    if(q<=alpha):
        FP+=1
        #print("Another FP:",q,"<=",alpha)
for q in qAlts:
    if(q<=alpha): TP+=1
FDR=float(FP)/float(FP+TP)
power=float(TP)/numAlts

# Compute binomial confidence interval for estimate of FDR
n=TP+FP
lower=FDR-1.96*math.sqrt(FDR*(1-FDR)/n)
upper=FDR+1.96*math.sqrt(FDR*(1-FDR)/n)
digits=3
print("FDR=",round(FDR,digits),", 95% CI=(",round(lower,digits),",",
      round(upper,digits),"), power=",round(power,digits),
      ", FP=",FP,", TP=",TP,sep="")

