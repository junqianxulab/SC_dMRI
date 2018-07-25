#!/usr/bin/env python

# exclude high b-value frame from DWI

import argparse
import os
import sys
import tempfile
import subprocess


def run_command(cmd, print_cmd = True):
    if print_cmd:
        print '>> %s' % cmd

    p = subprocess.call(cmd, shell=True)
    return p

def run_command_with_prefix(cmd, print_cmd = True, prefix="export FSLOUTPUTTYPE='NIFTI' ; "):
    if print_cmd is True:
        print '>> %s' % cmd

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

def split_dwi(filename, tmpdir=None):
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix='slice', dir='.')
    run_command_with_prefix('fslsplit %s %s/dwi_' % (filename, tmpdir))
    return tmpdir

def merge_frame(fn_out, lst_f, tmpdir):
    lst_fn_f = ['%s/dwi_%04d.nii' % (tmpdir, f) for f in lst_f]
    run_command('fslmerge -t %s %s' % (fn_out, ' '.join(lst_fn_f)))

def read_b_files(bn_b):
    with open(bn_b + '.bval') as fin:
        bval = [int(value) for value in fin.readline().strip().split()]
    with open(bn_b + '.bvec') as fin:
        bvec = [fin.readline().strip().split() for v in range(3)]
    return bval, bvec


def modify_b_files(bn_out, bval, bvec, lst_f):
    with open(bn_out + '.bval', 'w') as fout:
        fout.write(' '.join([str(bval[ind]) for ind in lst_f]))
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
    parser.add_argument('-o', help='output basename', dest='bn_out', default=None)
    parser.add_argument('-b0', help='b0 threshold', dest='b0', default=90)
    parser.add_argument('-high', help='high b value threshold', dest='high', default=400)
    ag = parser.parse_args()

    if ag.bn_out is None:
        ag.bn_out = os.path.basename(get_filename_wo_ext(ag.fn_dwi))

    thr_b0 = int(ag.b0)
    thr_high = int(ag.high)

    bval, bvec = read_b_files(ag.bn_b)

    tmpdir = split_dwi(ag.fn_dwi)
    lst_f = [f for f in range(len(bval)) if bval[f] < thr_b0 or bval[f] > thr_high ]
    modify_b_files(ag.bn_out + '_sub', bval, bvec, lst_f)
    merge_frame(ag.bn_out + '_sub', lst_f, tmpdir)

    run_command('rm %s -rf' % tmpdir)

