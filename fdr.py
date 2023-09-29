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
import statsmodels
from statsmodels import stats
from statsmodels.stats import multitest

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <swift2.out>\n")
(infile,)=sys.argv[1:]

pvalues=[]
IN=open(infile,"rt")
for line in IN:
    fields=line.rstrip().split()
    if(len(fields)!=7): continue
    (id,theta,Palt,CI_left,CI_right,median_p,area_p)=fields
    pvalues.append(float(median_p))
IN.close()

out=statsmodels.stats.multitest.fdrcorrection(pvalues,alpha=0.05,
                                               method='indep',
                                               is_sorted=False)
q=out[1]
for x in q:
    print(x)
    



