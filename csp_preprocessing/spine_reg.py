#!/usr/bin/env python

import os
import sys
import re
import subprocess
from tempfile import mkstemp
import time
#from dwi_utils import run_command
import dwi_utils
import shutil

run_command = dwi_utils.run_command

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from parameter import Parameter_reg2d

def split_slice(param, filename, basename):
    run_command('mkdir -p %s' % param.tempdir)
    run_command('fslslice %s %s' % (filename, os.path.join(param.tempdir, basename)))

    for i in range(param.ngroup):
        dir_slice = param.get_dir_slice(i)
        run_command('mkdir -p %s' % dir_slice)
        if param.multiband > 1:
            slclst = [ '%s_%s' % ( os.path.join(param.tempdir, basename), param.get_merge_slice_number(i+param.ngroup*j, merged=False)) for j in range(param.multiband) ]
            run_command('fslmerge -z %s_%s %s' % (os.path.join(dir_slice, basename), param.get_merge_slice_number(i), ' '.join(slclst)))
        else:
            run_command('mv %s_%s.nii.gz %s' % (os.path.join(param.tempdir, basename), param.get_merge_slice_number(i, merged=False), dir_slice))

def extract_frame(param, slc):
    '''
    extract_frame(param, slc)
    '''
    dir_slice = param.get_dir_slice(slc)
    for k in param.allframes:
        run_command('fslroi %s %s %s 1' % (param.get_fn('all', slc=slc, xenc=False), param.get_fn(slc=slc, frame=k, xenc=False), k))

def merge_frame(param, mode, slc=None, seg=None):
    '''
    merge_frame(param, mode, seg=0)
    mode: 'b0', 'segdwi', 'seggeom', 'all'

    '''
    alist = param.get_frame_list(mode, seg)
    lst = [ param.get_fn('frame', slc=slc, frame=k, xenc=True) for k in alist ]
    run_command('fslmerge -t %s %s' % (param.get_fn_slc_concate(mode, slc), ' '.join(lst)))

def merge_slice(param, mode, seg=None):
    '''
    merge_slice(param, mode, seg)
    mode: 'b0', 'segdwi', 'seggeom', 'all'
    '''

    if param.multiband > 1:
        slice_groups = [ [] for tmp in range(param.multiband) ]
        for i in range(param.ngroup):
            filename = param.get_fn_slc_concate(mode, i)
            basename = '%s_mbgrp_%04d' % (os.path.join(param.tempdir, param.bn_img), i)
            run_command('fslslice %s %s' % (filename, basename))
            for j in range(param.multiband):
                slice_groups[j].append('%s_slice_%04d.nii.gz' % (basename, j))
        slice_list = []
        for a_list in slice_groups:
            slice_list += a_list
    else:
        slice_list = [ param.get_fn_slc_concate(mode, slc) for slc in range(param.ngroup) ]
    if len(slice_list) > 1:
        run_command('fslmerge -z %s %s' % (param.get_fn(mode), ' '.join(slice_list)))
    else:
        run_command('cp %s %s' % (slice_list[0], param.get_fn(mode)))

def init(param, mode, slc, seg=None):
    '''
    init(param, mode, slc, seg=None)
    mode: 'b0', 'seggeom', 'segdwi'
    If mode is 'segdwi', set seg.
    '''
    alist = param.get_frame_list(mode, seg)

    for k in alist:
        mat = param.get_fn_mat(mode, frame=k, seg=seg, slc=slc) + '.mat'
        cmd = 'cp -f %s.nii.gz %s.nii.gz' % (param.get_fn(mode, k, slc=slc, xenc=False), param.get_fn(mode, k, slc=slc, xenc=True))
        run_command(cmd)
        run_command('cp -f %s %s' % (param.identmat, mat))

def calc(param, mode, slc, seg=None):
    '''
    calc(param, mode, slc, seg=None)
    mode: 'b0', 'seggeom, 'segdwi'
    If mode is 'seggdwi', set seg.
    '''
    alist = param.get_frame_list(mode, seg)
    mean_name = param.get_fn_mean(mode, seg=seg, slc=slc)

    lst = [ param.get_fn(mode, frame=k, slc=slc) for k in alist ]

    run_command('fslmerge -t %s_concat %s' % (mean_name, ' '.join(lst)))
    run_command('fslmaths %s_concat -Tmean %s' % (mean_name, mean_name))

def check_marker_filenames(filenames):
    block_filenames = filenames[:]

    i = 0
    while i < len(block_filenames):
        if os.path.isfile(block_filenames[i]):
            i += 1
        else:
            block_filenames.pop(i)

    return block_filenames

def run_command_bg(cmd, block_filenames=[], max_cpu=1, sleep_time=10, tempdir='.'):
    if '"' in cmd:
        print '" in command, %s' % cmd
        print 'operation runs as foreground not background'
        run_command(cmd)
        return block_filenames

    print max_cpu

    block_filenames = check_marker_filenames(block_filenames)
    while len(block_filenames) >= max_cpu:
        time.sleep(sleep_time)
        block_filenames = check_marker_filenames(block_filenames)

    # create temp file
    fd, temp_path = mkstemp(prefix='temp_marker_', dir=tempdir)
    os.close(fd)
    block_filenames.append(temp_path)
    run_command('run_with_marker_file.sh %s "%s" &' % (temp_path, cmd))
    return block_filenames

def xalign(param, mode, slc, seg=None):
    '''
    xalign(param, mode, slc, seg=None)
    mode: 'b0', 'segmean, 'segdwi'
    '''
    bg = False
    if bg:
        block_filenames = []

    alist = param.get_frame_list(mode, seg)
    mean_name = param.get_fn_mean(mode, seg=seg, slc=slc)

    for k in alist:
        # FIXME
        #cmd = 'flirt -in %s -inweight %s -ref %s -refweight %s -out %s -omat %s %s %s %s -interp sinc -cost mutualinfo' %\
        cmd = 'flirt -in %s -inweight %s -ref %s -refweight %s -out %s -omat %s %s %s %s -interp sinc' %\
                (param.get_fn(mode, frame=k, slc=slc), param.get_fn_mask(slc), mean_name, param.get_fn_mask(slc), param.get_fn(mode, frame=k, slc=slc), param.get_fn_tempmat(slc=slc, frame=k, seg=seg), param.schedule, param.noresample, param.nosearch)
        if bg:
            block_filenames = run_command_bg(cmd, block_filenames, max_cpu=param.max_cpu, sleep_time=param.sleep_time, tempdir=param.tempdir)
        else:
            run_command(cmd)
    if bg:
        while len(block_filenames) > 0:
            time.sleep(param.sleep_time)
            block_filenames = check_marker_filenames(block_filenames)

    for k in alist:
        mat = param.get_fn_mat(mode, frame=k, seg=seg, slc=slc) + '.mat'
        run_command('convert_xfm -omat %s -concat %s %s' % (mat, param.get_fn_tempmat(slc=slc, frame=k, seg=seg), mat))

        # skip?
        cmd = 'flirt -in %s -ref %s -out %s -applyxfm -init %s -interp sinc' %\
                (param.get_fn(mode, frame=k, slc=slc, xenc=False), mean_name, param.get_fn(mode, frame=k, slc=slc), mat)
        run_command(cmd)

def seggeommean_b0mean(param, slc):
    b0mean = param.get_fn_mean('b0', slc=slc)
    seggeommean = param.get_fn_mean('segmean', slc=slc)
    mask = param.get_fn_mask(slc=slc)
    mat = param.get_fn_mat('seggeom', slc=slc)
    #invmat = ...

    run_command('flirt -ref %s -refweight %s -in %s -inweight %s -out %s -omat %s %s %s -nosearch -interp sinc' %
            (b0mean, mask, seggeommean, mask, mat, param.get_fn_tempmat(slc=slc, seg='seggeommean_b0mean'), param.schedule, param.noresample) )
    mat += '.mat'
    run_command('cp -f %s %s' % (param.get_fn_tempmat(slc=slc, seg='seggeommean_b0mean'), mat))

    #invmat = '%s%s%s%s%s.mat' % (b0mean, slc, '_to_', os.path.basename(seggeommean), slc)
    #run_command('convert_xfm -omat %s -inverse %s' % (invmat, mat))
        
def compose_xfm(param, mode, slc=None, seg=None):
    '''
    compose_xfm(param, mode, slc=None, seg=None)
    mode: 'b0', 'segdwi'
    '''
    alist = param.get_frame_list(mode, seg)
    #mean_name = param.get_fn_mean(mode, seg, slc=slc)
    seggeommean = param.get_fn_mean('segmean', slc=slc)


    for k in alist:
        matAtoD = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
        #matAtoC = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
        matCtoD = param.get_fn_mat('seggeom', slc=slc) + '.mat'
        matBtoC = param.get_fn_mat('segmean', frame=seg, slc=slc) + '.mat'
        matAtoB = param.get_fn_mat(mode, frame=k, seg=seg, slc=slc) + '.mat'
        tempmat = param.get_fn_tempmat(slc=slc, frame=k, seg=seg)
        if mode == 'b0':
            run_command('cp -f %s %s' % (matAtoB, matAtoD))
        else:
            run_command('convert_xfm -omat %s -concat %s %s' %
                    (tempmat, matBtoC, matAtoB) )
            run_command('convert_xfm -omat %s -concat %s %s' %
                    (matAtoD, matCtoD, tempmat) )
        #run_command('flirt -in %s -ref %s -applyxfm -init %s -out %s -interp sinc' %
        #        (param.get_fn(mode, frame=k, slc=slc, xenc=False), seggeommean, matAtoD, param.get_fn(mode, frame=k, slc=slc)))

def apply_xfm(param, mode, slc=None, seg=None):
    '''
    apply_xfm(param, mode, slc=None, seg=None)
    mode: 'b0', 'segdwi'
    '''
    alist = param.get_frame_list(mode, seg)
    #mean_name = param.get_fn_mean(mode, seg, slc=slc)
    seggeommean = param.get_fn_mean('segmean', slc=slc)

    for k in alist:
        matAtoD = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
        run_command('flirt -in %s -ref %s -applyxfm -init %s -out %s -interp sinc' %
                (param.get_fn(mode, frame=k, slc=slc, xenc=False), seggeommean, matAtoD, param.get_fn(mode, frame=k, slc=slc)))

def get_xy_trans(param, mode, slc=None, seg=None):
    '''
    apply_xfm(param, mode, slc=None, seg=None)
    mode: 'b0', 'segdwi'
    '''

    alist = param.get_frame_list(mode, seg)
    xy_trans = [0.0 for i in range(len(alist))]

    # b0
    mode = 'b0'
    alist = param.get_frame_list(mode, seg)
    for k in alist:
        matAtoD = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
        with open(matAtoD) as fin:
            x = float(fin.readline().strip().split()[-1])
            y = float(fin.readline().strip().split()[-1])
            xy_trans[k] = x*x+y*y

    mode = 'segdwi'
    for seg in param.seglist:
        alist = param.get_frame_list(mode, seg)
        for k in alist:
            matAtoD = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
            with open(matAtoD) as fin:
                x = float(fin.readline().strip().split()[-1])
                y = float(fin.readline().strip().split()[-1])
                xy_trans[k] = x*x+y*y

    return xy_trans

def align(param, mode, slc, seg=None):
    '''
    align(param, mode, slc, seg=None)
    mode: 'b0', 'seggeom, 'segdwi'
    If mode is 'segdwi', set seg.
    '''
    init(param, mode, slc, seg)
    calc(param, mode, slc, seg)
    
    for t in range(param.nitr):
        xalign(param, mode, slc, seg)
        calc(param, mode, slc, seg)

def run_applywarp(param):
    for slc in range(param.ngroup):
        apply_xfm(param, 'b0', slc=slc)
        for j in param.seglist:
            apply_xfm(param, 'segdwi', slc=slc, seg=j)
        merge_frame(param, 'all', slc=slc)

    print '## Merge frames'
    merge_slice(param, 'all')
    if param.is_param_object:
        fn_out = param.get_fn('all')
        shutil.copy(fn_out + '.nii.gz', '../%s' % fn_out)
        return os.path.basename(fn_out) + '.nii.gz'

def run_registration(param):
    #
    print '## Extract slices'
    split_slice(param, param.fn_mas, param.bn_mas)
    split_slice(param, param.fn_img, param.bn_img)

    for slc in range(param.ngroup):
        # extract frames
        extract_frame(param, slc)

        # register individual dwi segment
        align(param, 'b0', slc)

        for j in param.seglist:
            align(param, 'segdwi', slc, seg=j)

        if len(param.seglist) > 0:
            align(param, 'segmean', slc)
            seggeommean_b0mean(param, slc)

        #compose transformations
        compose_xfm(param, 'b0', slc=slc)
        for j in param.seglist:
            compose_xfm(param, 'segdwi', slc=slc, seg=j)

        #apply xfm
        apply_xfm(param, 'b0', slc=slc)
        for j in param.seglist:
            apply_xfm(param, 'segdwi', slc=slc, seg=j)
        #apply_xfm(param, 'seggeom')

        merge_frame(param, 'all', slc=slc)
  
    print '## Merge frames'
    merge_slice(param, 'all')

    # remove intermediate files
    #run_command('rm -rf %s' % param.tempdir)

    # check the quality of registration by playing a movie in fsl
    #run_command('fslview %s_xenc.nii.gz' % param.p)
    if param.is_param_object:
        fn_out = param.get_fn('all')
        shutil.copy(fn_out + '.nii.gz', '../')
        return os.path.basename(fn_out) + '.nii.gz'
    else:
        print 'fslview %s_xenc.nii.gz' % param.bn_img

def generate_xy_trans(param, obj_return_value=None):
    xy_trans_group = []
    fn_out = param.get_fn_xy_trans()
    fout = open(fn_out, 'w')
    fout.write('%s,%s,%s\n' % ('mbgroup', 'frame', 'distance'))
    print 'Translation > 2.0 mm:'
    print '%s\t%s\t%s' % ('mbgroup', 'frame', 'distance')
    for slc in range(param.ngroup):
        xy_trans = get_xy_trans(param, 'all', slc=slc)
        xy_trans_group.append(xy_trans)
        for frame in range(len(xy_trans)):
            fout.write('%s,%s,%s\n' % (slc, frame, xy_trans[frame]))
            if xy_trans[frame] > 2.0:
                print '%s\t%s\t%s' % (slc, frame, xy_trans[frame])
    fout.close()

    if obj_return_value is not None:
        obj_return_value.delete(0, len(obj_return_value.get()))
        obj_return_value.insert(0, fn_out)
        return

    return xy_trans_group

def set_print_cmd(print_cmd):
    global run_command
    run_command = lambda x:dwi_utils.run_command(x, print_cmd=print_cmd)

if __name__ == '__main__':
    # start
    print '## Load a param file'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        print 'Usage: %s filename.params' % sys.argv[0]
        sys.exit(-1)

    param = Parameter_reg2d(paramfile=filename)
    #param = Param('CSPV032_dwi1_qa.params')
    #param.set_multiband(2)
    run_registration(param)


