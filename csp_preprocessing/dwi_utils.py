import shutil
import subprocess
import nibabel as nib
import sys
import os
import numpy as np

total_readout_time = 0.86 * (32-1) * 0.001
#gzip = 'gzip'
#gunzip = 'gunzip'

def gzip(filename):
    run_command('gzip -f %s' % filename)

def gunzip(filename):
    run_command('gunzip -c %s > %s' % (filename, filename[:-3]))

def append_log(filename, log):
    with open(filename, 'a') as fout:
        fout.write('%s\n' % log)

def run_command(cmd, print_cmd = True):
    if print_cmd is True:
        print('>> %s' % cmd)
    elif type(print_cmd) == type('') and print_cmd != '':
        append_log(print_cmd, cmd)

    p = subprocess.call(cmd, shell=True)
    return p

def filename_wo_ext(filename):
    if filename[-7:] == '.nii.gz':
        return filename[:-7]
    elif filename[-4:] == '.nii':
        return filename[:-4]
    elif filename[-3:] == '.gz':
        return filename[:-3]
    elif filename[-5:] == '.bval':
        return filename[:-5]
    elif filename[-5:] == '.bvec':
        return filename[:-5]
    elif filename[-9:] == '.4dfp.img':
        return filename[:-9]
    elif filename[-5:] == '.4dfp':
        return filename[:-5]

    return filename

def extname(filename):
    if filename[-7:] == '.nii.gz':
        return '.nii.gz'
    elif filename[-4:] == '.nii':
        return '.nii'
    elif filename[-3:] == '.gz':
        return '.gz'
    elif filename[-5:] == '.bval':
        return '.bval'
    elif filename[-5:] == '.bvec':
        return '.bvec'
    elif filename[-9:] == '.4dfp.img':
        return filename[-9:]
    elif filename[-5:] == '.4dfp':
        return filename[-5:]
    return ''

def create_merge(fn_output, lst_fn, is_dwi=False, verbose=False):
    run_command('fslmerge -t %s %s' % (fn_output, ' '.join(lst_fn)), verbose)
    if is_dwi:
        fn_out_wo_ext = filename_wo_ext(fn_output)
        lst_fn_wo_ext = [filename_wo_ext(filename) for filename in lst_fn]
        run_command('paste %s.bval > %s.bval' % ('.bval '.join(lst_fn_wo_ext), fn_out_wo_ext), verbose)
        run_command('paste %s.bvec > %s.bvec' % ('.bvec '.join(lst_fn_wo_ext), fn_out_wo_ext), verbose)


def create_acqparams(fn_output, lst_fn, direction=None, total_readout_time=total_readout_time):
    if direction is None:
        appa = [ 'ap' in filename.lower() for filename in lst_fn ]
    else:
        appa = [ 'AP' == tmp for tmp in direction ]

    with open(fn_output, 'w') as f:
        for i in range(len(lst_fn)):
            filename = lst_fn[i]

            if appa[i]:
                f.write('0 -1 0 %s\n' % total_readout_time)
                print('%s: AP' % filename)
            else:
                f.write('0 1 0 %s\n' % total_readout_time)
                print('%s: PA' % filename)

def create_index(fn_output, lst_fn_b0, lst_fn_dw, direction_b0, direction_dw):
    index = []
    i_dw = 0
    i_b0 = 0
    if direction_b0 is None:
        appa_b0 = [ 'ap' in filename.lower() for filename in lst_fn_b0 ]
    else:
        appa_b0 = [ 'AP' == tmp for tmp in direction_b0 ]
    if direction_dw is None:
        appa_dw = [ 'ap' in filename.lower() for filename in lst_fn_dw ]
    else:
        appa_dw = [ 'AP' == tmp for tmp in direction_dw ]

    while i_b0 < len(lst_fn_b0) and i_dw < len(lst_fn_dw):
        # update direction for new i_dwi
        if appa_dw[i_dw] == appa_b0[i_b0]:
            if appa_dw[i_dw]:
                direction = 'AP'
            else:
                direction = 'PA'
            img = nib.load(lst_fn_dw[i_dw])
            i_b0 += 1
            i_dw += 1
            index.append( ('%s ' % i_b0) * img.shape[3] )

            print('%s %s: %s' % (i_b0, direction, img.shape[3]))

        else:
            i_b0 += 1

    with open(fn_output, 'w') as f:
        f.write('%s\n' % ''.join(index))

def bval_bvec_to_4dfp(bval, bvec, outputname='bval_bvec'):
    with open(bval) as f:
        bval = f.readline().strip().split()
    with open(bvec) as f:
        bvecx = f.readline().strip().split()
        bvecy = f.readline().strip().split()
        bvecz = f.readline().strip().split()
    with open(outputname, 'w') as f:
        f.write('%s\n' % len(bval))
        for i in range(len(bval)):
            f.write('%s\t%s\t%s\t%s\n' % (bval[i], bvecx[i], bvecy[i], bvecz[i]))

def niftigz_4dfp(filename, direction='n'):
    if direction == 'n':
        fn = filename_wo_ext(filename)
        run_command('nifti_4dfp -n %s %s' % (fn, fn))
        gzip(fn + '.nii')
        return fn
    else:
        if filename[-3:] == '.gz':
            is_gz = True
            gunzip(filename)
        else:
            is_gz = False

        fn = filename_wo_ext(filename)
        run_command('nifti_4dfp -4 %s %s' % (fn, fn))
        #run_command('%s %s' % (gzip, filename))
        run_command('rm %s -f' % (filename[:-3]))
        return fn

def name_inc_find_last(base, search_dir=False):
    i = 0
    if search_dir:
        f = os.path.isdir
    else:
        f = os.path.isfile

    fn = ''
    while True:
        fn = '%s_%02d' % (base, i)
        if f(fn) is True:
            i += 1
            continue
        break

    return fn

def generate_dti_maps(filename, bval_bvec=None, bval=None, bvec=None, prefix='', outlier=None):
    if bval_bvec is None and (bval is None or bvec is None):
        sys.stderr.write('Either of bval_bvec or (bval and bvec) should be given\n')
        return

    basedir = os.path.dirname(filename)
    if bval_bvec is None:
        bval_bvec = os.path.join(basedir, prefix+'bval_bvec')
        bval_bvec_to_4dfp(bval, bvec, bval_bvec)

    if extname(filename) == '.nii' or extname(filename) == '.nii.gz':
        fn_img = niftigz_4dfp(filename, '4')
    else:
        fn_img = filename_wo_ext(filename)

    fn_out = fn_img + '_dti'
    fn_res = fn_img + '_res'

    if outlier is not None and outlier is not False:
        name_inc = name_inc_find_last(fn_img + '_itr', search_dir=True)
        os.mkdir(name_inc)
        run_command('mv %s_dti* %s' % (fn_img, name_inc))
        run_command('mv %s_res* %s' % (fn_img, name_inc))
        run_command('mv %s/*mask.nii.gz ./' % (name_inc))

    if outlier is not None and outlier is not False:
        res_command = 'G'
        res_post = outlier if type(outlier) == type('') else 'badenc.dat'
    else:
        res_command = ''
        res_post = ''

    cmd = 'diff_4dfp_js -r%sNEFV3' % res_command
    run_command('%s %s %s %s' % (cmd, bval_bvec, fn_img, res_post))

    niftigz_4dfp(fn_out, 'n')
    niftigz_4dfp(fn_res, 'n')

    run_command('fslroi %s %s_MD     0 1' % (fn_out, fn_out))
    run_command('fslroi %s %s_RA     1 1' % (fn_out, fn_out))
    run_command('fslroi %s %s_FA     2 1' % (fn_out, fn_out))
    run_command('fslroi %s %s_Eig3   3 1' % (fn_out, fn_out))
    run_command('fslroi %s %s_Eig2   4 1' % (fn_out, fn_out))
    run_command('fslroi %s %s_Eig1   5 1' % (fn_out, fn_out))
    run_command('fslroi %s %s_RD     6 1' % (fn_out, fn_out))
    run_command('fslroi %s %s_Evec1  7 3' % (fn_out, fn_out))
    run_command('fslroi %s %s_Evec2 10 3' % (fn_out, fn_out))
    run_command('fslroi %s %s_Evec3 13 3' % (fn_out, fn_out))

    return fn_out

def create_md_mask(fn_md, fn_out=None, thr_min=0.1, thr_max=1.1):
    import scipy.ndimage

    #constants
    minimum_number_of_voxel = 3

    if fn_out is None:
        fn_out = fn_md
        if fn_out[-7:] == '.nii.gz':
            fn_out = fn_out[:-7]
        elif fn_out[-4:] == '.nii':
            fn_out = fn_out[:-4]
        fn_out = fn_out + '_mask.nii.gz'

    img_md = nib.load(fn_md)
    dat_md = img_md.get_data()

    # Generate mask
    s = [[0,1,0], [1,1,1], [0,1,0]]
    c_x = int(dat_md.shape[0] / 2)
    c_y = int(dat_md.shape[1] / 2)
    max_dist = c_x + c_y
    dat_mask_volume = np.zeros(dat_md.shape, dtype=np.int8)
    dat_roi_volume = dat_md.copy()
    for z in range(dat_md.shape[-1]):
        dat_slice = dat_roi_volume[:,:,z]
        dat_mask = np.zeros(dat_slice.shape, dtype=dat_mask_volume.dtype)

        # apply min thr
        dat_mask[dat_slice > 0] = 1
        dat_mask[dat_slice < thr_min] = 0

        # erode
        dat_mask_base = dat_mask.copy()
        dat_mask = scipy.ndimage.binary_erosion(dat_mask, iterations=1)

        # remove non-major segment
        lbl = scipy.ndimage.label(dat_mask, s)
        lbl_count = [ len((lbl[0]==tmp).nonzero()[0]) for tmp in range(1, lbl[1]+1) ]
        lbl_dist = []
        for i_lbl in range(1, lbl[1]+1):
            if lbl_count[i_lbl-1] < minimum_number_of_voxel:
                lbl_dist.append(max_dist)
            else:
                lbl_dist.append(min([ abs(c_x-tmp[0]) + abs(c_y-tmp[1]) for tmp in zip(*(lbl[0]==i_lbl).nonzero()) ]))

        #arg_max = np.argmax(lbl_count) + 1
        #dat_mask[ lbl[0] != arg_max ] = 0
        arg_min = np.argmin(lbl_dist) + 1
        dat_mask[ lbl[0] != arg_min ] = 0

        # dilate
        if False:
            dat_mask = scipy.ndimage.binary_propagation(dat_mask, mask=dat_mask_base)
            dat_mask = scipy.ndimage.binary_propagation(dat_mask, mask=dat_mask_base)

        # apply max thr
        dat_mask[dat_slice > thr_max] = 0

        if True:
            # erode
            dat_mask = scipy.ndimage.binary_erosion(dat_mask, iterations=1)

            # remove non-major segment
            lbl = scipy.ndimage.label(dat_mask, s)
            lbl_count = [ len((lbl[0]==tmp).nonzero()[0]) for tmp in range(1, lbl[1]+1) ]
            lbl_dist = []
            for i_lbl in range(1, lbl[1]+1):
                if lbl_count[i_lbl-1] < minimum_number_of_voxel:
                    lbl_dist.append(max_dist)
                else:
                    lbl_dist.append(min([ abs(c_x-tmp[0]) + abs(c_y-tmp[1]) for tmp in zip(*(lbl[0]==i_lbl).nonzero()) ]))

            #arg_max = np.argmax(lbl_count) + 1
            #dat_mask[ lbl[0] != arg_max ] = 0
            if len(lbl_dist) > 0:
                arg_min = np.argmin(lbl_dist) + 1
                dat_mask[ lbl[0] != arg_min ] = 0
            else:
                dat_mask[:] = 0

        dat_slice[dat_mask == 0]= 0
        dat_mask_volume[:,:,z] = dat_mask

    # save
    hdr = img_md.get_header()
    hdr.set_data_dtype(dat_mask_volume.dtype)
    img_out = nib.Nifti1Image(dat_mask_volume, img_md.get_affine(), hdr)
    nib.save(img_out, fn_out)

    return fn_out

def calculate_voxel_values(filenames, filename_roi, filename_csv=None):
    img_roi = nib.load(filename_roi)
    roi = img_roi.get_data()

    if type(filenames) == type(''):
        filenames = [filenames]

    nz = roi.shape[2]
    nf = len(filenames)

    stats = [ ('mean', np.mean), ('std', np.std), ('min', np.min), ('max', np.max), ('median', np.median) ]
    ns = len(stats)

    stat_matrix = np.zeros((nf*ns, nz+1), dtype=float)

    if filename_csv is not None:
        fout = open(filename_csv, 'w')
        fout.write('filename,%s,total\n' % ','.join([str(z) for z in range(nz)]))

    for f in range(nf):
        dat = nib.load(filenames[f]).get()
        for z in range(nz):
            dat_in = dat[:,:,z][roi[:,:,z]]
            for s in range(len(stats)):
                title, stat = stats[s]
                stat_matrix[f*ns+s,z] = stat(dat_in)

        dat_in = dat[roi]
        for s in range(len(stats)):
            title, stat = stats[s]
            stat_matrix[f*ns+s,-1] = stat(dat_in)

            if filename_csv is not None:
                fout.write('%s_%s,%s\n' % (filenames[f], title, ','.join( [str(tmp) for tmp in stat_matrix[f*ns+s,:] ] )))


    if filename_csv is not None:
        fout.close()
    return stat_matrix

def abspath_to_relpath(src, dst):
    sep = os.path.sep
    src_l = src.strip(sep).split(sep)
    dst_l = dst.strip(sep).split(sep)

    i = 0
    while i < len(src_l) and i < len(dst_l) and src_l[i] == dst_l[i]:
        i += 1

    rel_path = sep.join( ['..' for tmp in range(len(src_l)-i)] + dst_l[i:] )
    return rel_path


