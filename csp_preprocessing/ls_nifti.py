#!/usr/bin/env python

import glob
import sys
import os

if len(sys.argv) > 1:
    fn = sys.argv[1]
else:
    fn = '../E*/*_DWI_*_MB2*/*.nii.gz'

def argsort(seq):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__)

lst = glob.glob(fn)
lst_num = [ int( fn.split(os.sep)[2].split('_')[0] ) for fn in lst ]
lst_sorted = [ lst[ ind ] for ind in argsort(lst_num) ]

for fn in lst_sorted:
    print(fn)

