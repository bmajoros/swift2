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
import ProgramName
from scipy.stats import binom

def computePvalue(alpha,beta,k,m):
    n=k+m
    testMean=(alpha+k)/(alpha+k+beta+m)
    moreExtreme=0
    for x in range(n+1):
        weight=binom.pmf(x,n,0.5)
        nullMean=(alpha+x)/(alpha+x+beta+n-x)
        if(nullMean<=testMean): moreExtreme+=weight
    return (testMean,moreExtreme)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <k> <m> <alpha> <beta>\n")
(k,m,alpha,beta)=sys.argv[1:]
k=int(k)
m=int(m)
alpha=float(alpha)
beta=float(beta)

(testMean,pValue)=computePvalue(alpha,beta,k,m)
binomTest=binom.cdf(k,k+m,0.5)
print(round(testMean,3),round(pValue,6),round(binomTest,6))






