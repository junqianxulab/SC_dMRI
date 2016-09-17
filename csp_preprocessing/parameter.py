import re
import nibabel as nib
import os.path
import socket
import sys
import argparse

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-fn_img', dest='fn_img')
    parser.add_argument('-fn_mas', dest='fn_mas')
    parser.add_argument('-nitr', dest='nitr')
    parser.add_argument('-multiband', dest='multiband')
    parser.add_argument('-b0frames', dest='b0frames')
    parser.add_argument('-b0_threshold', dest='b0_threshold')
    parser.add_argument('-schedule', dest='schedule')
    parser.add_argument('-noresample', dest='noresample')
    parser.add_argument('-nosearch', dest='nosearch')
    parser.add_argument('-dir_slice', dest='dir_slice')
    parser.add_argument('-max_cpu', dest='max_cpu')
    parser.add_argument('-sleep_time', dest='sleep_time')
    parser.add_argument('-verbose', dest='verbose')
    return parser

###
# subject
# working_dir
# fn_bval
# fn_bvec
# b0_threshold

# topup_app
# topup_config
# eddy_app

# fn_b0
# fn_b0_mask
# fn_dwi
# fn_dwi_mask
# fn_dwi_eddy_out

# multiband
# niter
# nosearch
# noresample

# fn_dwi_eddy
# fn_reg_mask
# fn_reg_outlier
# fn_reg_out

# fn_dwi_dti
# fn_dwi_dti_mask
# fn_dwi_dti_outlier

# fn_dti_value
# fn_dti_value_roi
# fn_dti_value_out

class Parameter:
    def __init__(self):
        self.subject = ''
        self.working_dir = ''
        self.fn_bval = ''
        self.fn_bvec = ''
        self.b0_threshold = ''

        self.topup_app = ''
        self.topup_config =''
        self.eddy_app = ''

        self.fn_b0 = ''
        self.fn_b0_mask = ''
        self.fn_dwi = ''
        self.fn_dwi_mask = ''
        self.fn_dwi_eddy_out = ''

        self.multiband = ''
        self.niter = ''
        self.nosearch = ''
        self.noresample = ''
        self.schedule = ''

        self.fn_dwi_eddy = ''
        self.fn_reg_mask = ''
        self.fn_reg_outlier = ''
        self.fn_reg_out = ''
        self.fn_dwi_dti = ''
        self.fn_dti = ''
        self.fn_dwi_dti_mask = ''
        self.fn_dwi_dti_outlier = ''
        self.fn_dti_value = ''
        self.fn_dti_value_roi = ''
        self.fn_dti_value_out = ''

        self.set_variable_types()

    def set_variable_types(self):
        try:
            self.b0_threshold = int(self.b0_threshold)
        except:
            self.b0_threshold = 0

        try:
            self.multiband = int(self.multiband)
        except:
            self.multiband = 0

        try:
            self.niter = int(self.niter)
        except:
            self.niter = 0

        self.nosearch = (self.nosearch == True or self.nosearch == 'True')
        self.noresample = (self.noresample == True or self.noresample == 'True')
        

    def save(self, filename):
        variables = vars(self)
        with open(filename, 'w') as fout:
            for a_key in variables:
                fout.write("%s = %s\n" % (str(a_key), str(variables[a_key])))

    def read(self, filename):
        variables = vars(self)
        with open(filename) as fin:
            for line in fin.readlines():
                a_pair = re.sub(r'( )*=( )*', '\n', line.strip()).split('\n')
                variables[a_pair[0]] = a_pair[1]
        self.set_variable_types()


class Parameter_reg2d:
    def __init__(self, fn_img=None, fn_mas=None, nitr=3, multiband=None, b0frames=None, b0_threshold=90,
                 schedule=None, noresample=False, nosearch=True, dir_slice=True, max_cpu=1, sleep_time=30, verbose=True,
                 paramfile=None, args=None, param_object=None):
        '''
        Parameter(
            fn_image=None,
            fn_mask=None,
            niter=3,
            multiband=None,
            b0frames=None,
            b0_threshold=90,
            schedule=None,
            noresample=True,
            dir_slice = True,
            max_cpu = 1,
            sleep_time = 30,
            verbose=True)
        '''
        self.fn_img = fn_img
        self.fn_mas = fn_mas
        self.nitr = nitr
        self.multiband = multiband
        self.b0frames = b0frames
        self.b0_threshold = b0_threshold
        self.schedule = schedule
        self.noresample = noresample
        self.nosearch = nosearch
        self.dir_slice = dir_slice
        self.max_cpu = max_cpu
        self.sleep_time = sleep_time
        self.verbose = verbose
      
        if paramfile is not None:
            self.read_param(paramfile)
        if args is not None:
            self.parse_args(args)

        if param_object is not None:
            self.read_parameter_object(param_object)
            self.is_param_object = True
        else:
            self.is_param_object = False

        if self.fn_img is not None and self.fn_mas is not None and self.b0frames is not None:
            self.set_variables()
        print vars(self)

    def read_parameter_object(self, param_object):
        dirname = 'xy_reg'
        os.chdir(param_object.working_dir)

        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        os.chdir(dirname)

        self.fn_img = os.path.join('..', param_object.fn_dwi_eddy)
        self.fn_mas = os.path.join('..', param_object.fn_reg_mask)
        self.nitr = param_object.niter
        self.multiband = param_object.multiband
        self.b0frames = os.path.join('..', param_object.fn_bval)
        self.b0_threshold = param_object.b0_threshold
        self.schedule = param_object.schedule
        self.noresample = param_object.noresample
        self.nosearch = param_object.nosearch

    def read_param(self, filename):
        fin = open(filename)
        args = [ '-%s' % (tmp.strip()) for tmp in fin.readlines() if tmp != '\n' and tmp[0] != '#' ]
        self.parse_args(re.sub(r'( )*=( )*', '=', '\n'.join(args)).split('\n'))
        fin.close()

    def save_param(self, filename):
        with open(filename, 'w') as fout:
            variables = [
                    ('fn_img', self.fn_img),
                    ('fn_mas', self.fn_mas),
                    ('nitr', self.nitr),
                    ('multiband', self.multiband),
                    ('b0frames', self.b0frames),
                    ('b0_threshold', self.b0_threshold),
                    ('schedule', self.schedule),
                    ('noresample', self.noresample),
                    ('nosearch', self.nosearch),
                    ('dir_slice', self.dir_slice),
                    ('max_cpu', self.max_cpu),
                    ('sleep_time', self.sleep_time),
                    ('verbose', self.verbose),
                    ]
            for var in variables:
                if var[1] is not None:
                    fout.write('%s = %s\n' % (var[0], var[1]))

    def parse_args(self, args):
        parser = create_parser()
        ag = parser.parse_args(args)
        parse_dict = vars(ag)
        self_dict = vars(self)
        for key in parse_dict:
            if parse_dict[key] is not None:
                self_dict[key] = parse_dict[key]

    def set_variables(self):
        self.nitr = int(self.nitr)
        #self.b0frames = b0frames
        self.b0_threshold = int(self.b0_threshold)
        if type(self.verbose) == type(''):
            if self.verbose[0].lower() == 't':
                self.verbose = True
            else:
                self.verbose = False

        test_miss = ['fn_img', 'fn_mas', 'b0frames']
        test_miss = [tmp for tmp in test_miss if tmp not in vars(self)]
        if len(test_miss) > 0:
            sys.stderr.write('%s variable(s) should be set first\n' % (', '.join(test_miss)))
            return

        self.identmat = os.path.join('$FSLDIR','etc', 'flirtsch', 'ident.mat')
        self.tempmat = 'temp.mat'
        self.temp1mat = 'temp1.mat'
        self.tempdir= 'temp.dir'
        self.tag_frame = 'frame'
        self.tag_xenc = 'xenc'

        self.read_fn_img(self.fn_img)
        self.set_multiband(self.multiband)
        if self.multiband > 1:
            self.ngroup = self.nslice / self.multiband
        else:
            self.ngroup = self.nslice

        if self.fn_img[-7:] == '.nii.gz':
            self.bn_img = os.path.basename(self.fn_img)[:-7]
        elif self.fn_img[-4:] == '.nii':
            self.bn_img = os.path.basename(self.fn_img)[:-4]
        else:
            self.bn_img = os.path.basename(self.fn_img)

        if self.fn_mas[-7:] == '.nii.gz':
            self.bn_mas = os.path.basename(self.fn_mas)[:-7]
        elif self.fn_mas[-4:] == '.nii':
            self.bn_mas = os.path.basename(self.fn_mas)[:-4]
        else:
            self.bn_mas = os.path.basename(self.fn_mas)


        if type(self.b0frames) == type('str'):
            if os.path.isfile(self.b0frames):
                self.b0frames = get_b0_from_bval(self.b0frames, threshold=self.b0_threshold, return_verbose=False, verbose=self.verbose)
            elif len(self.b0frames.split(' ')) > 1:
                self.b0frames = [int(tmp) for tmp in self.b0frames.split(' ')]
            elif len(self.b0frames.split(',')) > 1:
                self.b0frames = [int(tmp) for tmp in self.b0frames.split(',')]
            else:
                sys.stderr.write('invalid b0frames or file not exist.: %s\n' % self.b0frames)
                sys.exit(-1)

        self.segframes = set_seg_frames(self.b0frames, self.nvolume)
        self.seglist = range(len(self.segframes))
        self.allframes = range(self.nvolume)

        if type(self.noresample) == type(''):
            self.noresample = (self.noresample[0].lower() == 't')
        if self.noresample:
            self.noresample = '-noresample'
        else:
            self.noresample = ''

        if type(self.nosearch) == type(''):
            self.nosearch = (self.nosearch[0].lower() == 't')
        if self.nosearch:
            self.nosearch = '-nosearch'
        else:
            self.nosearch = ''

        if type(self.max_cpu) == type(''):
            self.max_cpu = int(self.max_cpu)

        if type(self.sleep_time) == type(''):
            self.sleep_time = int(self.sleep_time)

        if self.schedule is None:
            dirname = os.path.dirname(os.path.realpath(__file__))
            fn_schedule = os.path.join(dirname, 'xytrans.sch')
            if not os.path.isfile(fn_schedule):
                if False:
                    hostname = socket.gethostname()
                    if 'some' == hostname[:4].lower():
                        fn_schedule = '/home/xugroup/bin/xytrans.sch'
                    elif 'mssm-jw' == hostname[:7].lower():
                        fn_schedule = '/Users/joowon/research/code/spineDTI/xytrans.sch'
                    else:
                        sys.stderr.write('schedule file not found\n')
                        sys.exit(-1)
                sys.stderr.write('schedule file not found\n')
                sys.exit(-1)
            self.schedule = '-schedule %s' % fn_schedule
        elif self.schedule[:9] != '-schedule':
            self.schedule = '-schedule ' + self.schedule

    def set_multiband(self, multiband):
        if type(multiband) == type(''):
            multiband = int(multiband)
        if type(multiband) == type(1) and multiband > 1:
            self.multiband = multiband
        else:
            self.multiband = 0

    def read_fn_img(self, filename):
        img = nib.load(filename)
        shape = img.shape
        if len(shape) > 3:
            self.nvolume = shape[3]
            self.nslice = shape[2]
        else:
            self.nvolume = shape[2]
            self.nslice = 1
        
        del img

    def get_fn_tempmat(self, slc=None, frame=None, seg=None):
        fn_lst = [ os.path.join(self.tempdir, self.bn_img) ]
        if slc is not None:
            fn_lst.append('_slice_%s' % slc)
        if seg is not None:
            fn_lst.append('_seg_%s' % seg)
        if frame is not None:
            fn_lst.append('_frame_%s' % frame)
        fn_lst.append('.mat')
        return ''.join(fn_lst)

    def get_fn(self, mode='frame', frame=None, slc=None, xenc=True, merged=True):

        if slc is not None:
            fn_lst = [ os.path.join(self.get_dir_slice(slc), self.bn_img) ]
            fn_lst.append(self.get_merge_slice_number(slc, merged))
        else:
            fn_lst = [self.bn_img]

        if mode == 'b0mean':            # n = 1
            fn_lst.append('b0mean')
        elif mode == 'segmean':         # n = n_seg
            fn_lst.append('segmean')
        elif mode == 'seggeom':         # n = 1
            fn_lst.append('seggeom')
        elif mode == 'all':
            pass
        elif frame is not None:                           # n = n_frame
            fn_lst.append(self.tag_frame)

        if frame is not None:
            fn_lst.append(str(frame))

        if xenc:
            fn_lst.append(self.tag_xenc)

        return '_'.join(fn_lst)

    def get_fn_slc_concate(self, mode='frame', slc=None, merged=True, xenc=True):
        fn_lst = [ os.path.join(self.get_dir_slice(slc), self.bn_img) ]
        fn_lst.append(self.get_merge_slice_number(slc, merged))
        fn_lst.append(self.tag_xenc)
        return '_'.join(fn_lst)

    def get_fn_mean(self, mode, seg=None, slc=None, xenc=False):
        if mode == 'b0':
            return self.get_fn('b0mean', xenc=xenc, slc=slc)
        elif mode == 'segdwi':
            return self.get_fn('segmean', xenc=xenc, frame=seg, slc=slc)
        elif mode == 'segmean':
            return self.get_fn('seggeom', xenc=xenc, slc=slc)
        elif mode == 'seggeom':
            return self.get_fn('b0mean', xenc=xenc, slc=slc)

    def get_fn_mask(self, slc=None, merged=True):
        return '%s_%s' % (os.path.join(self.get_dir_slice(slc), self.bn_mas), self.get_merge_slice_number(slc, merged))

    def get_merge_slice_number(self, slc, merged=True):
        if slc is None:
            return ''
        if self.multiband > 1 and merged is True:
            return 'mbgrp_%04d' % slc
        else:
            return 'slice_%04d' % slc

    def get_frame_list(self, mode, seg=None):
        if mode == 'b0':
            return self.b0frames
        elif mode == 'segdwi':
            return self.segframes[seg]
        elif mode == 'segmean':
            return self.seglist
        elif mode == 'all':
            return self.allframes

    def get_dir_slice(self, slc):
        if slc is None:
            return os.path.curdir
        elif self.dir_slice:
            return '%s' % (self.get_merge_slice_number(slc))
        return os.path.curdir

    def get_fn_mat(self, mode, frame=None, seg=None, slc=None):
        if mode == 'b0' or mode == 'segdwi' or mode == 'segmean':
            return '%s_to_%s' % (self.get_fn(mode=mode, frame=frame, slc=slc, xenc=False, merged=True), os.path.basename(self.get_fn_mean(mode=mode, seg=seg, slc=slc)))
        elif mode == 'b0_concate':
            return '%s_to_%s' % (self.get_fn(mode='b0', frame=frame, slc=slc, xenc=False, merged=True), os.path.basename(self.get_fn(mode='seggeom', slc=slc)))
        elif mode == 'segdwi_concate':
            return '%s_to_%s' % (self.get_fn(mode='segdwi', frame=frame, slc=slc, xenc=False, merged=True), os.path.basename(self.get_fn(mode='seggeom', slc=slc)))
        elif mode == 'seggeom':
            return '%s_to_%s' % (self.get_fn(mode='seggeom', slc=slc, xenc=False, merged=True), os.path.basename(self.get_fn_mean(mode='seggeom', slc=slc)))

    def get_fn_xy_trans(self):
        return os.path.basename(self.bn_img) + '_xy_trans.csv'

def get_b0_from_bval(filename, threshold=90, return_verbose=False, verbose=True):
    fin = open(filename)
    frames = [int(value) for value in fin.readline().strip().split()]
    b0frames = []
    rtn_str = '# '
    for i in range(len(frames)):
        if frames[i] < threshold:
            b0frames.append(i)
            rtn_str += str('(%s):%s ' % (i,frames[i]))
        else:
            rtn_str += str('%s:%s ' % (i,frames[i]))
        if (i+1)%5 == 0:
            rtn_str += '\n# '
    rtn_str += '\n'

    fin.close()
    if verbose:
        print rtn_str
    if return_verbose:
        return b0frames, rtn_str
    return b0frames

def set_seg_frames(b0frames, nvolume):
    if len(b0frames) < 1:
        return range(nvolume)
    segframes = [ range(1+int(b0frames[j-1]), int(b0frames[j]))
            for j in range(1, len(b0frames)) if int(b0frames[j])-int(b0frames[j-1]) > 1]
    if int(b0frames[-1]) < int(nvolume)-1:
        segframes.append( range(1+int(b0frames[-1]), int(nvolume) ))

    return segframes

def get_frames_fram_bval(filename, threshold_min=90, threshold_max=400, return_verbose=False, verbose=True):
    fin = open(filename)
    frames = [int(value) for value in fin.readline().strip().split()]
    target_frames = []
    rtn_str = '# '
    for i in range(len(frames)):
        if threshold_min <= frames[i] <= threshold_max:
            target_frames.append(i)
            rtn_str += str('(%s):%s ' % (i,frames[i]))
        else:
            rtn_str += str('%s:%s ' % (i,frames[i]))
        if (i+1)%5 == 0:
            rtn_str += '\n# '
    rtn_str += '\n'

    fin.close()
    if verbose:
        print rtn_str
    if return_verbose:
        return target_frames, rtn_str
    return target_frames


