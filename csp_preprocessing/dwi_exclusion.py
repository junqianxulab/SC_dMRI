#!/usr/bin/env python

# create / modify exclusion file

# input csv file:
# slice,frame
#
# e.g.
# slice,frame
# 0, 3
# 2, 1
#
#
# exclusion file:
# 1st row: n_slice n_frame (the number of slice and the number of frame)
# 2nd to 1+n_slice row: n_frame 0 or 1's (0: normal frame, 1: bad frame)
#
# e.g.
# 3 5
# 0 0 1 0 0
# 0 0 0 0 0
# 0 1 0 0 0
# 

import argparse
import os
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='input csv filename', dest='fn_csv', required=True)
    parser.add_argument('-e', help='exclusion slice/frame filename', dest='fn_exc', default='badend.dat')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-o', help='overwrite(creating new) existing file', dest='overwrite', action='store_true')
    group.add_argument('-m', help='modify(append) existing file', dest='modify', action='store_true')
    parser.add_argument('-d', help='dwi filename for #slice/#frame extraction', dest='fn_dwi')

    ag = parser.parse_args()

    print ag
    # if exclusion file exists,
    if os.path.isfile(ag.fn_exc) and ag.overwrite is False and ag.modify is False:
        s = raw_input('%s exists. [o]verwrite/[m]odify/e[x]it? ' % ag.fn_exc)
        if s[0].lower() == 'o':
            ag.overwrite = True
        elif s[0].lower() == 'm':
            ag.modify = True
        else:
            exit(0)

    if not os.path.isfile(ag.fn_csv):
            sys.stderr.write('fn_csv %s not exist.\n' % ag.fn_csv)
            sys.exit(-1)
    # 
    if ag.modify is True and os.path.isfile(ag.fn_exc):
        with open(ag.fn_exc) as fin:
            n_slice, n_frame = [int(value) for value in fin.readline().strip().split()]
            bad_frames = [ fin.readline().strip().split() for s in range(n_slice) ]
    else:
        if ag.fn_dwi is None:
            sys.stderr.write('fn_dwi should be set unless modifying an existing exclustion file.\n')
            sys.exit(-1)
        if not os.path.isfile(ag.fn_dwi):
            sys.stderr.write('fn_dwi %s not exist.\n' % ag.fn_dwi)
            sys.exit(-1)
        import nibabel as nib
        img = nib.load(ag.fn_dwi)
        n_slice, n_frame = img.shape[2:]
        bad_frames = [ [ '0' for f in range(n_frame) ] for s in range(n_slice) ]

    print '# slices = %s, # frames = %s' % (n_slice, n_frame)

    with open(ag.fn_csv) as fin:
        line = fin.readline() # header
        line = fin.readline()
        while line:
            s, f = [int(value) for value in line.strip().split(',')]
            bad_frames[s][f] = '1'
            line = fin.readline()

    with open(ag.fn_exc, 'w') as fout:
        fout.write('%s %s\n' % (n_slice, n_frame))
        fout.write('\n'.join( [ ' '.join(s) for s in bad_frames ] ) )

