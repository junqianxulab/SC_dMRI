#!/usr/bin/env python

# exclude bad slice/frame from DWI
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
import tempfile
import subprocess



def run_command(cmd, print_cmd = True):
    if print_cmd:
        print('>> %s' % cmd)

    p = subprocess.call(cmd, shell=True)
    return p

def run_command_with_prefix(cmd, print_cmd = True, prefix="export FSLOUTPUTTYPE='NIFTI' ; "):
    if print_cmd is True:
        print('>> %s' % cmd)

    p = subprocess.call(prefix + cmd, shell=True)
    return p

def get_filename_wo_ext(filename):
    rtn = filename
    if rtn[-3:] == '.gz':
        rtn = rtn[:-3]
    if rtn[-4:] == '.nii':
        rtn = rtn[:-4]
    return rtn

def read_excludion(filename):
    with open(filename) as fin:
        nz, nf = [int(value) for value in fin.readline().strip().split()]
        bad_frames = [ fin.readline().strip().split() for z in range(nz) ]
    return nz, nf, bad_frames

def split_dwi(filename, nz, tmpdir=None):
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix='slice', dir='.')
    run_command_with_prefix('fslslice %s %s/dwi' % (filename, tmpdir))
    for i in range(nz):
        fn_in = '%s/dwi_slice_%04d' % (tmpdir, i)
        run_command_with_prefix('fslsplit %s %s_' % (fn_in, fn_in))
    return tmpdir

def merge_frame(fn_out, z, lst_f, tmpdir):
    lst_fn_f = ['%s/dwi_slice_%04d_%04d.nii' % (tmpdir, z, f) for f in lst_f]
    run_command('fslmerge -t %s %s' % (fn_out, ' '.join(lst_fn_f)))

def read_b_files(bn_b):
    with open(bn_b + '.bval') as fin:
        bval = fin.readline().strip().split()
    with open(bn_b + '.bvec') as fin:
        bvec = [fin.readline().strip().split() for v in range(3)]
    return bval, bvec

def modify_b_files(bn_out, bval, bvec, lst_f):
    with open(bn_out + '.bval', 'w') as fout:
        fout.write(' '.join([bval[ind] for ind in lst_f]))
    with open(bn_out + '.bvec', 'w') as fout:
        fout.write(' '.join([bvec[0][ind] for ind in lst_f]))
        fout.write('\n')
        fout.write(' '.join([bvec[1][ind] for ind in lst_f]))
        fout.write('\n')
        fout.write(' '.join([bvec[2][ind] for ind in lst_f]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='DWI filename', dest='fn_dwi', required=True)
    parser.add_argument('-b', help='bval/bvec basename', dest='bn_b', required=True)
    parser.add_argument('-e', help='exclusion slice/frame filename', dest='fn_exc', required=True)
    parser.add_argument('-o', help='output basename', dest='bn_out', default=None)

    ag = parser.parse_args()

    if ag.bn_out is None:
        ag.bn_out = os.path.basename(get_filename_wo_ext(ag.fn_dwi))

    bval, bvec = read_b_files(ag.bn_b)
    nz, nf, bad_frames = read_excludion(ag.fn_exc)

    tmpdir = split_dwi(ag.fn_dwi, nz)
    if nz > 1:
        for z in range(nz):
            lst_f = [f for f in range(nf) if bad_frames[z][f] == '0']
            modify_b_files(ag.bn_out + '_%04d_sub' % z, bval, bvec, lst_f)
            merge_frame(ag.bn_out + '_%04d_sub' % z, z, lst_f, tmpdir)
    else:
        lst_f = [f for f in range(nf) if bad_frames[0][f] == '0']
        modify_b_files(ag.bn_out + '_sub', bval, bvec, lst_f)
        merge_frame(ag.bn_out + '_sub', 0, lst_f, tmpdir)

    
    run_command('rm %s -rf' % tmpdir)

