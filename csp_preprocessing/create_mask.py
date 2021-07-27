#!/usr/bin/env python

import sys
import re
import subprocess
import os
import numpy as np
import nibabel as nib
from dwi_utils import run_command, filename_wo_ext

def nbd_4(ij, shape):
    nbds = [ (-1,0), (1,0), (0,-1), (0,1) ]
    return [ (ij[0] + nbd[0], ij[1] + nbd[1]) for nbd in nbds if 0 <= ij[0] + nbd[0] < shape[0] if 0 <= ij[1] + nbd[1] < shape[1] ]

def nbd_8(ij, shape):
    nbds = [ (-1,0), (1,0), (0,-1), (0,1),
            (-1,-1), (1,1), (1,-1), (-1,1) ]
    return [ (ij[0] + nbd[0], ij[1] + nbd[1]) for nbd in nbds if 0 <= ij[0] + nbd[0] < shape[0] if 0 <= ij[1] + nbd[1] < shape[1] ]

def center_of_mass_2d(matrix_2d):
    center_x = 0.0
    center_y = 0.0
    sum_intensity = 0.0
    for i in range(matrix_2d.shape[0]):
        for j in range(matrix_2d.shape[1]):
            if matrix_2d[i,j] > 0:
                #value = np.log(matrix_2d[i,j])
                value = matrix_2d[i,j]
                center_x += (i*value)
                center_y += (j*value)
                sum_intensity += value

    center_x /= sum_intensity
    center_y /= sum_intensity

    return int(np.ceil(center_x)), int(np.ceil(center_y))

def dilate_threshold(base, threshold, dat, k):
    to_do = base
    done = []
    while to_do:
        done += to_do
        ing = to_do[:]
        to_do = []

        for ij in ing:
            for nbd in nbd_8(ij, dat.shape):
                if nbd not in to_do and nbd not in ing and nbd not in done and \
                        dat[nbd[0],nbd[1],k] > threshold:
                    to_do.append(nbd)
                    #dat_th[nbd[0],nbd[1],k] = 1
    return done
    
def dilate_nth(base, nth, shape):
    done = base[:]
    to_do = base[:]
    for i in range(nth):
        ing = to_do[:]
        to_do = []

        for ij in ing:
            for nbd in nbd_4(ij, shape):
                if nbd not in to_do and nbd not in ing and nbd not in done:
                    to_do.append(nbd)
        done += to_do
    return done

def simple_hist(data, bins, normed=False, histtype='step', rangex=None):
    if type(bins) == type(0):
        if rangex:
            bin = np.linspace( rangex[0], rangex[1], bins+1 )
        else:
            bin = np.linspace( min(data), max(data), bins+1 )
    else:
        bin = np.array(bins)

    if normed:
        y = np.zeros( (len(bin)-1), dtype=float )
    else:
        y = np.zeros( (len(bin)-1), dtype=int )


    for value in data:
        if value < bin[0] or value > bin[-1]:
            continue
        done = False
        for i in range(1, len(bin)):
            if value < bin[i]:
                y[i-1] += 1
                done = True
                break
        if not done:    # maybe value == bin[-1]
            y[-1] += 1

    if normed:
        y /= sum(y)
        y /= np.array([bins[i+1]-bins[i] for i in range(len(bins)-1)])

    return y, bin, None

def make_base_image(filename, filename_fltd, b0frames):
    if b0frames:
        img = nib.load(filename)
        dat = img.get_data()
        mean_b0 = np.zeros( (dat.shape[0], dat.shape[1], dat.shape[2]), dtype=float)
        for b0 in b0frames:
            mean_b0 += dat[:,:,:,b0]

        mean_b0 /= len(b0frames)
        hdr = img.get_header()
        hdr.set_data_dtype(float)

        img_b0 = nib.Nifti1Image(mean_b0, img.get_affine(), hdr)
        nib.save(img_b0, filename_fltd)
        
    else:
        run_command('fslroi %s %s 0 1' % (filename, filename_fltd))
    run_command('fslmaths %s -kernel 2D -fmean %s' % (filename_fltd, filename_fltd))

def set_filenames(filename):
    filename_woext = filename_wo_ext(filename)
    filename_wodir = os.path.basename(filename_woext)

    filename_fltd = '%s_mask_base_filtered' % filename_wodir
    filename_mask = '%s_mask'        % filename_wodir

    return filename, filename_fltd, filename_mask

def calculate_com(dat, zooms):
    shape = dat.shape
    dxy = int(min([shape[0], shape[1]]) / 2)
    range_x = int(shape[0]/2 - dxy, shape[0]/2 + dxy)
    range_y = int(shape[1]/2 - dxy, shape[1]/2 + dxy)

    dx_mm = 15
    dy_mm = 15
    dx = int(dx_mm / zooms[0])
    dy = int(dy_mm / zooms[1])

    roi_x = [[]] * shape[2]
    roi_y = [[]] * shape[2]
    com_x = [[]] * shape[2]
    com_y = [[]] * shape[2]


    for k in range(shape[2]):
        com_x[k], com_y[k] = center_of_mass_2d( dat[range_x[0]:range_x[1], range_y[0]:range_y[1], k] )
        com_x[k] += range_x[0]
        com_y[k] += range_y[0]

        roi_x[k] = [max([0, com_x[k] - dx]), min(shape[0], com_x[k] + dx)]
        roi_y[k] = [max([0, com_y[k] - dy]), min(shape[1], com_y[k] + dy)]

    return com_x, com_y, roi_x, roi_y

def fill_in(dat_th, done, k):
    expanded_done = []
    for row in range(dat_th.shape[1]):
        index = [ i for (i,j) in done if j == row ]
        if not index:
            continue

        for col in range(min(index), max(index)+1):
            dat_th[col,row,k] = 1
            expanded_done.append((col,row))

    for col in range(dat_th.shape[0]):
        index = [ j for (i,j) in done if i == col ]
        if not index:
            continue

        for row in range(min(index), max(index)+1):
            dat_th[col,row,k] = 1
            expanded_done.append((col,row))
    return expanded_done

def show_stat(dat_th):
    for k in range(dat_th.shape[2]):
        marked = (dat_th[:,:,k] == 1).nonzero()
        done = zip(marked[0], marked[1])
        done_expand = dilate_nth(done, 3, dat_th.shape)

        inside  = np.array([dat[ij[0], ij[1], k] for ij in done])
        outside = np.array([dat[ij[0], ij[1], k] for ij in done_expand if ij not in done])

        print(' #mask voxel: %s, #outside voxel: %s' % (len(inside), len(outside)))
        print(' Inside mask:  mean=%s, std=%s' % (inside.mean(), inside.std()))
        print(' Outside mask: mean=%s, std=%s' % (outside.mean(), outside.std()))

def mask_from_threshold(com_x, com_y, threshold, dat, dat_th, k):
    done = dilate_threshold([(com_x, com_y)], threshold, dat, k)
    inside  = np.array([dat[ij[0], ij[1], k] for ij in done])
    done = dilate_nth(done, 2, dat.shape)
    done = fill_in(dat_th, done, k)

# start
if __name__ == '__main__':
    if len(sys.argv) > 1:
        filename, filename_fltd, filename_mask = set_filenames(sys.argv[1])
    else:
        print('Usage: %s filename b0frames(1st frame is 0)' % sys.argv[0])
        sys.exit(-1)

    print('Input: %s' % filename)
    print('Output: %s' % filename_mask)

    make_base_image(filename, filename_fltd+'.nii.gz', sys.argv[2:])

    img = nib.load(filename_fltd+'.nii.gz')
    dat = img.get_data()
    hdr = img.get_header()
    zooms = hdr.get_zooms()
    shape = dat.shape
    dat_th = np.zeros(shape, dtype='int16')

    com_x, com_y, roi_x, roi_y = calculate_com(dat, zooms)

    plt = None
    for k in range(shape[2]):
        dat_subflt = dat[roi_x[k][0]:roi_x[k][1], roi_y[k][0]:roi_y[k][1], k].flatten()
        if plt:
            if k < maxsub:
                plt.subplot('%s%s' % (subs,k+1))
                plt.hist(dat_subflt, bins=40)
            if not k < maxsub:
                plt.subplot('%s%s' % (subs,k-5))
                plt.hist(dat_subflt, bins=40)

        #threshold = np.mean(dat_subflt)
        threshold = np.mean(dat_subflt) + np.std(dat_subflt)
        mask_from_threshold(com_x[k], com_y[k], threshold, dat, dat_th, k)

        print('slice: %s, center = (%s, %s) voxel, threshold = %s' % (k, com_x[k], com_y[k], threshold))

    img_out = nib.Nifti1Image(dat_th, img.get_affine(), hdr)
    nib.save(img_out, filename_mask)

    if plt:
        import matplotlib.pyplot as plt
        subs = '23'
        maxsub = 6
        plt.figure( figsize = (10,5) )

        for k in range(shape[2]):
            dat_subflt = dat[roi_x[k][0]:roi_x[k][1], roi_y[k][0]:roi_y[k][1], k].flatten()
            if plt:
                if k < maxsub:
                    plt.subplot('%s%s' % (subs,k+1))
                    plt.hist(dat_subflt, bins=40)
                if not k < maxsub:
                    plt.subplot('%s%s' % (subs,k-5))
                    plt.hist(dat_subflt, bins=40)

        plt.figure()
        for k in range(shape[2]):
            if k < maxsub:
                plt.subplot('%s%s' % (subs, k+1))
                plt.imshow(dat[roi_x[k][0]:roi_x[k][1], roi_y[k][0]:roi_y[k][1], k].T, origin='lower', cmap='Greys_r')

        plt.figure( figsize = (10,5) )
        for k in range(shape[2]):

            if k < maxsub:
                plt.subplot('%s%s' % (subs,k+1))
                plt.imshow(dat[:,:,k].T, origin='lower', cmap='Greys_r')

        plt.figure( figsize = (10,5) )
        for k in range(shape[2]):

            if k < maxsub:
                plt.subplot('%s%s' % (subs,k+1))
                plt.imshow(dat_th[:,:,k].T, origin='lower', cmap='Greys_r')

        plt.show()

