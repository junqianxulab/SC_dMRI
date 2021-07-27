#!/usr/bin/env python

import sys
import os

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s basename [tol=0.15]\n' % os.path.basename(sys.argv[0]))
    sys.exit(-1)

if len(sys.argv) > 2:
    tol = float(sys.argv[2])
else:
    tol = 0.15

with open(sys.argv[1] + '.bval') as f:
    bvals = [val for val in f.readline().strip().split()]

with open(sys.argv[1] + '.bvec') as f:
    bvecs = [[float(val) for val in line.strip().split()] for line in f.readlines()]

for i in range(len(bvals)):
    if -tol < bvecs[2][i] < tol:
        print i, bvals[i], bvecs[0][i], bvecs[1][i], bvecs[2][i]

