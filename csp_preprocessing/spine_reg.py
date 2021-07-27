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
import nibabel as nib
import numpy as np

#run_command = dwi_utils.run_command
def run_command_with_prefix(cmd, print_cmd = True, prefix="export FSLOUTPUTTYPE='NIFTI' ; "):
    if print_cmd is True:
        print '>> %s' % cmd
    elif type(print_cmd) == type('') and print_cmd != '':
        dwi_utils.append_log(print_cmd, cmd)

    p = subprocess.call(prefix + cmd, shell=True)
    return p

run_command = run_command_with_prefix

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from parameter import Parameter_reg2d

def split_slice_wo_fsl(param, filename, basename):
    run_command('mkdir -p %s' % param.tempdir)
    img = nib.load(filename)
    shape = list(img.shape)
    dat = img.get_data()
    if len(shape) > 3:
        is_4d = True
    else:
        is_4d = False

    #run_command('fslslice %s %s' % (filename, os.path.join(param.tempdir, basename)))
    for z in range(shape[2]):
        if is_4d:
            img_out = nib.Nifti1Image(dat[:,:,z,:].reshape(shape[:2]+[1]+[shape[3]]), img.affine, img.header)
        else:
            img_out = nib.Nifti1Image(dat[:,:,z], img.affine, img.header)
        nib.save(img_out, os.path.join(param.tempdir, '%s_%04d.nii' % (basename, z)))


    if param.multiband > 1:
        if is_4d:
            dat_sub = np.empty(shape[:2] + [param.multiband] + [shape[3]])
        else:
            dat_sub = np.empty(shape[:2] + [param.multiband])
        for i in range(param.ngroup):
            dir_slice = param.get_dir_slice(i)
            run_command('mkdir -p %s' % dir_slice)
            slclst = [ '%s_%s' % ( os.path.join(param.tempdir, basename), param.get_merge_slice_number(i+param.ngroup*j, merged=False)) for j in range(param.multiband) ]
            #run_command('fslmerge -z %s_%s %s' % (os.path.join(dir_slice, basename), param.get_merge_slice_number(i), ' '.join(slclst)))
            for j in range(param.multiband):
                dat_sub[:,:,j] = dat[:,:,i+param.ngroup*j]
            img_out = nib.Nifti1Image(dat_sub, img.affine, img.header)
            nib.save(img_out, '%s_%s.nii' % (os.path.join(dir_slice, basename), param.get_merge_slice_number(i)))

    else:
        for i in range(param.ngroup):
            dir_slice = param.get_dir_slice(i)
            run_command('mkdir -p %s' % dir_slice)
            run_command('mv %s_%s.nii %s' % (os.path.join(param.tempdir, basename), param.get_merge_slice_number(i, merged=False), dir_slice))

    return img

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
            run_command('mv %s_%s.nii %s' % (os.path.join(param.tempdir, basename), param.get_merge_slice_number(i, merged=False), dir_slice))
            #run_command('mv %s_%s.nii.gz %s' % (os.path.join(param.tempdir, basename), param.get_merge_slice_number(i, merged=False), dir_slice))

def extract_frame_wo_fsl(param, slc):
    '''
    extract_frame(param, slc)
    '''
    dir_slice = param.get_dir_slice(slc)
    img = nib.load('%s.nii' % param.get_fn('all', slc=slc, xenc=False))
    shape = list(img.shape)
    dat = img.get_data().reshape(shape[:2]+shape[-1:])
    for k in param.allframes:
        #run_command('fslroi %s %s %s 1' % (param.get_fn('all', slc=slc, xenc=False), param.get_fn(slc=slc, frame=k, xenc=False), k))
        img_out = nib.Nifti1Image(dat[:,:,k], img.affine, img.header)
        nib.save(img_out, '%s.nii' % param.get_fn(slc=slc, frame=k, xenc=False))

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

def merge_slice(param, mode, seg=None, existing=[]):
    '''
    merge_slice(param, mode, seg)
    mode: 'b0', 'segdwi', 'seggeom', 'all'
    '''

    fn_out = param.get_fn(mode)
    if param.multiband > 1:
        slice_groups = [ [] for tmp in range(param.multiband) ]
        for i in range(param.ngroup):
            # merge with existing registration
            if i in existing:
                for j in range(param.multiband):
                    slice_groups[j].append('%s/%s_slice_%04d.nii' % (param.tempdir, os.path.basename(fn_out), i+j*param.ngroup))
                continue
            filename = param.get_fn_slc_concate(mode, i)
            basename = '%s_mbgrp_%04d' % (os.path.join(param.tempdir, param.bn_img), i)
            run_command('fslslice %s %s' % (filename, basename))
            for j in range(param.multiband):
                slice_groups[j].append('%s_slice_%04d.nii' % (basename, j))
                #slice_groups[j].append('%s_slice_%04d.nii.gz' % (basename, j))
        slice_list = []
        for a_list in slice_groups:
            slice_list += a_list
    else:
        slice_list = [ param.get_fn_slc_concate(mode, slc) for slc in range(param.ngroup) ]
        # merge with existing registration
        for i in existing:
            slice_list[i] = '%s/%s_slice_%04d.nii' % (param.tempdir, os.path.basename(fn_out), i)
    if len(slice_list) > 1:
        run_command('fslmerge -z %s %s' % (param.get_fn(mode), ' '.join(slice_list)))
    else:
        run_command('cp %s %s' % (slice_list[0], fn_out))

def init(param, mode, slc, seg=None):
    '''
    init(param, mode, slc, seg=None)
    mode: 'b0', 'segdwi', 'lowsegdwi', 'segmean'
    If mode is 'segdwi', set seg.
    '''
    alist = param.get_frame_list(mode, seg)

    for k in alist:
        mat = param.get_fn_mat(mode, frame=k, seg=seg, slc=slc) + '.mat'
        cmd = 'cp -f %s.nii %s.nii' % (param.get_fn(mode, k, slc=slc, xenc=False), param.get_fn(mode, k, slc=slc, xenc=True))
        #cmd = 'cp -f %s.nii.gz %s.nii.gz' % (param.get_fn(mode, k, slc=slc, xenc=False), param.get_fn(mode, k, slc=slc, xenc=True))
        run_command(cmd)
        run_command('cp -f %s %s' % (param.identmat, mat))

def calc(param, mode, slc, seg=None, itr=None):
    '''
    calc(param, mode, slc, seg=None)
    mode: 'b0', 'segdwi', 'lowsegdwi', 'b0mean', 'segmean', 'lowsegmean'
    If mode is 'b0', 'seggdwi' or 'lowsegdwi', set seg.
    '''
    alist = param.get_frame_list(mode, seg)
    mean_name = param.get_fn_mean(mode, seg=seg, slc=slc)

    lst = [ param.get_fn(mode, frame=k, slc=slc) for k in alist ]

    run_command('fslmerge -t %s_concat %s' % (mean_name, ' '.join(lst)))
    run_command('fslmaths %s_concat -Tmean %s' % (mean_name, mean_name))
    if itr is not None:
        run_command('cp %s_concat.nii %s_concat_itr%s.nii' % (mean_name, mean_name, itr))
        run_command('cp %s.nii %s_itr%s.nii' % (mean_name, mean_name, itr))

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

def xalign(param, mode, slc, seg=None, itr=None):
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
                (param.get_fn(mode, frame=k, slc=slc), param.get_fn_mask(slc), mean_name, param.get_fn_mask(slc), param.get_fn(mode, frame=k, slc=slc), param.get_fn_tempmat(slc=slc, frame=k, seg=mode+str(seg)), param.schedule, param.noresample, param.nosearch)
        if bg:
            block_filenames = run_command_bg(cmd, block_filenames, max_cpu=param.max_cpu, sleep_time=param.sleep_time, tempdir=param.tempdir)
        else:
            run_command(cmd)
            if itr is not None:
                run_command('cp %s %s' % (param.get_fn_tempmat(slc=slc, frame=k, seg=mode+str(seg)), param.get_fn_tempmat(slc=slc, frame=k, seg=mode+str(seg)) + '_itr%s' % itr) )
    if bg:
        while len(block_filenames) > 0:
            time.sleep(param.sleep_time)
            block_filenames = check_marker_filenames(block_filenames)

    for k in alist:
        mat = param.get_fn_mat(mode, frame=k, seg=seg, slc=slc) + '.mat'
        run_command('convert_xfm -omat %s -concat %s %s' % (mat, param.get_fn_tempmat(slc=slc, frame=k, seg=mode+str(seg)), mat))

        # skip?
        cmd = 'flirt -in %s -ref %s -out %s -applyxfm -init %s -interp sinc' %\
                (param.get_fn(mode, frame=k, slc=slc, xenc=False), mean_name, param.get_fn(mode, frame=k, slc=slc), mat)
        run_command(cmd)
        if itr is not None:
            run_command('cp %s %s' % (param.get_fn(mode, frame=k, slc=slc)+'.nii', param.get_fn(mode, frame=k, slc=slc)+'_itr%s'%itr+'.nii'))

#def seggeommean_b0mean(param, slc):
def seggeommean_b0geom(param, mode, slc):
    '''
    mode: high or low
    '''
    b0geom = param.get_fn_mean('b0mean', slc=slc)
    geommean = param.get_fn(mode, slc=slc, xenc=False)
    mask = param.get_fn_mask(slc=slc)
    mat = param.get_fn_mat(mode, slc=slc)
    #invmat = ...

    run_command('flirt -ref %s -refweight %s -in %s -inweight %s -out %s -omat %s %s %s -nosearch -interp sinc' %
            (b0geom, mask, geommean, mask, mat, param.get_fn_tempmat(slc=slc, seg=mode+'_b0geom'), param.schedule, param.noresample) )
    mat += '.mat'
    run_command('cp -f %s %s' % (param.get_fn_tempmat(slc=slc, seg=mode+'_b0geom'), mat))

    #invmat = '%s%s%s%s%s.mat' % (b0mean, slc, '_to_', os.path.basename(geommean), slc)
    #run_command('convert_xfm -omat %s -inverse %s' % (invmat, mat))
        
def compose_xfm(param, mode, slc=None, seg=None):
    '''
    compose_xfm(param, mode, slc=None, seg=None)
    mode: 'b0', 'segdwi', 'lowsegdwi'
    '''
    alist = param.get_frame_list(mode, seg)
    #mean_name = param.get_fn_mean(mode, seg, slc=slc)
    #seggeommean = param.get_fn_mean('segmean', slc=slc)


    for k in alist:
        matAtoD = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
        #matAtoC = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
        if mode == 'b0':
            matBtoD = param.get_fn_mat('b0mean', frame=seg, slc=slc) + '.mat'
            matAtoB = param.get_fn_mat(mode, frame=k, seg=seg, slc=slc) + '.mat'
            tempmat = param.get_fn_tempmat(slc=slc, frame=k, seg=mode+str(seg))

        elif mode == 'segdwi':
            matCtoD = param.get_fn_mat('seggeom', slc=slc) + '.mat'
            matBtoC = param.get_fn_mat('segmean', frame=seg, slc=slc) + '.mat'
            matAtoB = param.get_fn_mat(mode, frame=k, seg=seg, slc=slc) + '.mat'
            tempmat = param.get_fn_tempmat(slc=slc, frame=k, seg=mode+str(seg))
        elif mode == 'lowsegdwi':
            matCtoD = param.get_fn_mat('lowseggeom', slc=slc) + '.mat'
            matBtoC = param.get_fn_mat('lowsegmean', frame=seg, slc=slc) + '.mat'
            matAtoB = param.get_fn_mat(mode, frame=k, seg=seg, slc=slc) + '.mat'
            tempmat = param.get_fn_tempmat(slc=slc, frame=k, seg=mode+str(seg))

        if mode == 'b0':
            #run_command('cp -f %s %s' % (matAtoB, matAtoD))
            run_command('convert_xfm -omat %s -concat %s %s' %
                    (matAtoD, matBtoD, matAtoB) )
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
    #seggeommean = param.get_fn_mean('segmean', slc=slc)
    b0geom = param.get_fn_mean('b0mean', slc=slc)

    for k in alist:
        matAtoD = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
        run_command('flirt -in %s -ref %s -applyxfm -init %s -out %s -interp sinc' %
                (param.get_fn(mode, frame=k, slc=slc, xenc=False), b0geom, matAtoD, param.get_fn(mode, frame=k, slc=slc)))
        #       (param.get_fn(mode, frame=k, slc=slc, xenc=False), seggeommean, matAtoD, param.get_fn(mode, frame=k, slc=slc)))

def get_xy_trans(param, mode, slc=None, seg=None):
    '''
    apply_xfm(param, mode, slc=None, seg=None)
    mode: 'b0', 'segdwi'
    '''

    alist = param.get_frame_list(mode, seg)
    xy_trans = [0.0 for i in range(len(alist))]
    x_trans = [0.0 for i in range(len(alist))]
    y_trans = [0.0 for i in range(len(alist))]

    # b0
    mode = 'b0'
    alist = param.get_frame_list(mode, seg)
    for k in alist:
        matAtoD = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
        with open(matAtoD) as fin:
            x = float(fin.readline().strip().split()[-1])
            y = float(fin.readline().strip().split()[-1])
            x_trans[k] = x
            y_trans[k] = y
            xy_trans[k] = x*x+y*y

    mode = 'segdwi'
    for seg in param.seglist:
        alist = param.get_frame_list(mode, seg)
        for k in alist:
            matAtoD = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
            with open(matAtoD) as fin:
                x = float(fin.readline().strip().split()[-1])
                y = float(fin.readline().strip().split()[-1])
                x_trans[k] = x
                y_trans[k] = y
                xy_trans[k] = np.sqrt( x*x+y*y )

    mode = 'lowsegdwi'
    for seg in param.lowseglist:
        alist = param.get_frame_list(mode, seg)
        for k in alist:
            matAtoD = param.get_fn_mat('%s_concate' % mode, frame=k, slc=slc) + '.mat'
            with open(matAtoD) as fin:
                x = float(fin.readline().strip().split()[-1])
                y = float(fin.readline().strip().split()[-1])
                x_trans[k] = x
                y_trans[k] = y
                xy_trans[k] = np.sqrt( x*x+y*y )

    return xy_trans, x_trans, y_trans

def align(param, mode, slc, seg=None):
    '''
    align(param, mode, slc, seg=None)
    mode: 'b0', 'segdwi', 'lowsegdwi', 'b0mean', 'segmean', 'lowsegmean'
    If mode is 'segdwi', set seg.
    '''
    init(param, mode, slc, seg)
    calc(param, mode, slc, seg)
    
    for t in range(param.nitr):
        xalign(param, mode, slc, seg, t)
        calc(param, mode, slc, seg, t)

def run_applywarp(param):
    for slc in range(param.ngroup):
        apply_xfm(param, 'b0', slc=slc)
        for j in param.seglist:
            apply_xfm(param, 'segdwi', slc=slc, seg=j)
        merge_frame(param, 'all', slc=slc)

    print '## Merge frames'
    merge_slice(param, 'all')
    dwi_utils.gzip(param.get_fn('all') + '.nii')

    if param.is_param_object:
        fn_out = param.get_fn('all')
        shutil.copy(fn_out + '.nii.gz', '../%s' % fn_out)
        return os.path.basename(fn_out) + '.nii.gz'

def run_registration(param, sub_slices=None):
    #
    if sub_slices is None:
        lst_slices = range(param.ngroup)
    else:
        # sub slices
        # there should be exisiting registration result
        fn_reg_out = param.get_fn('all')
        if not (os.path.isfile(fn_reg_out) or os.path.isfile(fn_reg_out+'.nii') or os.path.isfile(fn_reg_out+'.nii.gz')):
            sys.stderr.write('%s not exist\n' % fn_reg_out)
            sys.exit(-1)

        if param.multiband > 1:
            lst_slices = []
            for s in sub_slices:
                s = int(s)
                if s > param.ngroup and s not in lst_slices:
                    lst_slices.append(s)
                elif s not in lst_slices:
                    lst_slices.append(s)
        else:
            lst_slices = [int(s) for s in sub_slices]

    print '## Extract slices'
    split_slice(param, param.fn_mas, param.bn_mas)
    split_slice(param, param.fn_img, param.bn_img)

    #for slc in range(param.ngroup):
    for slc in lst_slices:
        # extract frames
        extract_frame(param, slc)

        # register individual dwi segment
        #align(param, 'b0', slc)

        for j in param.b0list:
            align(param, 'b0', slc, seg=j)

        if len(param.b0list) > 0:
            align(param, 'b0mean', slc)

        for j in param.seglist:
            align(param, 'segdwi', slc, seg=j)

        if len(param.seglist) > 0:
            align(param, 'segmean', slc)
            seggeommean_b0geom(param, 'seggeom', slc)

        for j in param.lowseglist:
            align(param, 'lowsegdwi', slc, seg=j)

        if len(param.lowseglist) > 0:
            align(param, 'lowsegmean', slc)
            seggeommean_b0geom(param, 'lowseggeom', slc)

        #compose transformations
        #compose_xfm(param, 'b0', slc=slc)

        for j in param.b0list:
            compose_xfm(param, 'b0', slc=slc, seg=j)

        for j in param.seglist:
            compose_xfm(param, 'segdwi', slc=slc, seg=j)

        for j in param.lowseglist:
            compose_xfm(param, 'lowsegdwi', slc=slc, seg=j)

        #apply xfm
        for j in param.b0list:
            apply_xfm(param, 'b0', slc=slc, seg=j)
        for j in param.seglist:
            apply_xfm(param, 'segdwi', slc=slc, seg=j)
        for j in param.lowseglist:
            apply_xfm(param, 'lowsegdwi', slc=slc, seg=j)
        #apply_xfm(param, 'seggeom')

        merge_frame(param, 'all', slc=slc)
  
    print '## Merge frames'
    if sub_slices is not None:
        existing = [slc for slc in range(param.ngroup) if slc not in lst_slices]
        basename = os.path.basename(fn_reg_out)
        run_command('fslslice %s %s/%s' % (fn_reg_out, param.tempdir, basename))
    else:
        existing = []

    merge_slice(param, 'all', existing=existing)
    dwi_utils.gzip(param.get_fn('all') + '.nii')


    # remove intermediate files
    run_command('rm -rf %s' % param.tempdir)

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
    fout.write('%s,%s,%s,%s,%s\n' % ('mbgroup', 'frame', 'x_trans', 'y_trans', 'distance'))
    print 'Translation > 2.0 mm:'
    print '%s\t%s\t%s\t%s\t%s' % ('mbgroup', 'frame', 'x_trans', 'y_trans', 'distance')
    for slc in range(param.ngroup):
        xy_trans, x_trans, y_trans = get_xy_trans(param, 'all', slc=slc)
        xy_trans_group.append(xy_trans)
        for frame in range(len(xy_trans)):
            fout.write('%s,%s,%s,%s,%s\n' % (slc, frame, x_trans[frame], y_trans[frame], xy_trans[frame]))
            if xy_trans[frame] > 2.0:
                print '%s\t%s\t%s\t%s\t%s' % (slc, frame, x_trans[frame], y_trans[frame], xy_trans[frame])
    fout.close()

    if obj_return_value is not None:
        obj_return_value.delete(0, len(obj_return_value.get()))
        obj_return_value.insert(0, fn_out)
        return

    return xy_trans_group

def set_print_cmd(print_cmd):
    global run_command
    run_command = lambda x:run_command_with_prefix(x, print_cmd=print_cmd)

def generate_xy_trans_lazy(base_directory, fn_out='xy_reg_distance.csv'):
    import glob
    lst = glob.glob(os.path.join(base_directory, '*_*', '*seggeom_xenc.mat'))

    def extract_mbgroup_slice(filename):
        group_dn, fn = filename.split(os.sep)[-2:]
        group = group_dn.split('_')[1]

        ind_to = fn.rfind('_to_')
        ind_frame = fn[:ind_to].rfind('_frame_')
        frame = fn[ind_frame+7:ind_to]
        try:
            group = int(group)
            frame = int(frame)
        except:
            return None
        return (group, frame, filename),
    
    lst_g_f_fn = [ ext for filename in lst for ext in extract_mbgroup_slice(filename) if ext is not None ]
    lst_g_f_fn.sort()

    with open(fn_out, 'w') as fout:
        fout.write('%s,%s,%s,%s,%s\n' % ('mbgroup', 'frame', 'x_trans', 'y_trans', 'distance'))
        print 'Translation > 2.0 mm:'
        print '%s\t%s\t%s\t%s\t%s' % ('mbgroup', 'frame', 'x_trans', 'y_trans', 'distance')
        for (group, frame, filename) in lst_g_f_fn:
            with open(filename) as fin:
                x = float(fin.readline().strip().split()[-1])
                y = float(fin.readline().strip().split()[-1])
                x_trans = x
                y_trans = y
                xy_trans = np.sqrt(x*x+y*y)

            fout.write('%s,%s,%s,%s,%s\n' % (group, frame, x_trans, y_trans, xy_trans))
            if xy_trans > 2.0:
                print '%s\t%s\t%s\t%s\t%s' % (group, frame, x_trans, y_trans, xy_trans)


def usage():
    sys.stderr.write('Usage: %s filename.params\n' % os.path.basename(sys.argv[0]))
    sys.stderr.write('Usage: %s --xytrans base_directory [output_filename=xy_reg_distance.csv]\n' % os.path.basename(sys.argv[0]))

if __name__ == '__main__':
    # start
    if len(sys.argv) < 2:
        usage()
        sys.exit(-1)
    
    if sys.argv[1][:2] == '--':
        if sys.argv[1] == '--xytrans':
            if len(sys.argv) > 3:
                generate_xy_trans_lazy(sys.argv[2], sys.argv[3])
            elif len(sys.argv) > 2:
                generate_xy_trans_lazy(sys.argv[2])
            else:
                usage()
                sys.exit(-1)

    else:
        filename = sys.argv[1]
        print '## Load a param file'
        param = Parameter_reg2d(paramfile=filename)
        #param = Param('CSPV032_dwi1_qa.params')
        #param.set_multiband(2)

        sub_slices = None if len(sys.argv) == 2 else [int(s) for s in sys.argv[2:]]
        print 'sub_slices:', sub_slices
        run_registration(param, sub_slices=sub_slices)


