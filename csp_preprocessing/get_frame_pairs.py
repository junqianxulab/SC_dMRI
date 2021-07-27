#!/usr/bin/env python

import sys
import os

if len(sys.argv) < 5:
    sys.stderr.write('Usage: %s dwi_filename frame_1 frame_2 n_frame [n_dwi]\n' % os.path.basename(sys.argv[0]))
    sys.exit(-1)

fn = sys.argv[1]
frame_1 = int(sys.argv[2])
frame_2 = int(sys.argv[3])
n_frame = int(sys.argv[4])

if len(sys.argv) > 5:
    n_dwi = int(sys.argv[5])
else:
    import nibabel as nib
    n_dwi = int(nib.load(sys.argv[1]).shape[-1] / n_frame)

to_extract = [frame_1, frame_2]
merged_info = [ i%2 for i in range(n_dwi) ]

merged_frame_to_extract_pairs = [
        i*n_frame + to_extract[merged_info[i]] for i in range(n_dwi)
        ]

print(merged_frame_to_extract_pairs)

bn = os.path.basename(fn).split('.')[0]

cmds = [
        'fslroi %s %s_rep_%s %s 1' % (fn, bn, i, merged_frame_to_extract_pairs[i]) for i in range(n_dwi)
        ]

cmds += [
        'fslmerge -t %s_rep %s' % (bn, ' '.join(
            ['%s_rep_%s' % (bn, i) for i in range(n_dwi)]))
        ]

#cmd += ['fslmaths %s_rep -mas %s %s_rep_masked' % (bn, fn_mask, bn)]

for cmd in cmds:
    print cmd


