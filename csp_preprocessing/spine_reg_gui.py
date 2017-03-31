#!/usr/bin/env python

from Tkinter import Tk, Entry, TOP, BOTH, RIGHT, RAISED, X, LEFT, N, E, W, S, EW, NS, NSEW, RIDGE, Toplevel
from Tkinter import Grid, BooleanVar, Spinbox
from ttk import Frame, Button, Style, Label, Entry, Checkbutton
import tkFileDialog
import tkMessageBox
import prepare_dwi_gui
from dwi_utils import extname, filename_wo_ext, run_command, generate_dti_maps, create_md_mask, calculate_voxel_values
import os
import sys
import parameter
import time
import spine_reg

#filenameDialog_text = prepare_dwi_gui.filenameDialog_text
#dirnameDialog_text = prepare_dwi_gui.dirnameDialog_text

def append_text(fn, text, header=''):
    with open(fn, 'a') as fout:
        fout.write('#%s\n' % header)
        fout.write('%s\n\n' % text)


class SpineRegGui(Frame):

    def __init__(self, parent, param_filename=None):
        Frame.__init__(self, parent)

        self.title = 'SpineDTI'
        self.log = os.path.join(os.path.abspath('.'), 'log')
        self.parent = parent
        self.parameter_values = {}
        self.param = parameter.Parameter()
        self.initUI()
        self.init_variable()
        self.update_param_from_text()
        if param_filename is not None:
            self.read_param(param_filename)

    def init_variable(self):
        #self.txt_multiband.insert(0, '2')
        self.reset_entry(self.spin_multiband, 2)
        #self.txt_nitr.insert(0, '3')
        self.reset_entry(self.spin_nitr, 3)
        self.chk_nosearch.invoke()
        self.reset_entry(self.txt_b0_threshold, 90)

        dirname = os.path.dirname(os.path.realpath(__file__))
        self.reset_entry(self.txt_topup_config, os.path.join(dirname, 'b02b0_ON_1_3mm.cnf'))
        self.reset_entry(self.txt_reg_config, os.path.join(dirname, 'xytrans.sch'))
        self.reset_entry(self.txt_topup_app, 'topup_2d')
        self.reset_entry(self.txt_eddy_app, 'eddy')

        if self.txt_working.get() == '':
            self.reset_entry(self.txt_working, os.path.abspath('.'), force=True)

    def update_parameter(self):
        if self.parameter_values.has_key('subject'):
            self.reset_entry(self.txt_subject, self.parameter_values['subject'])
        if self.parameter_values.has_key('dwi'):
            self.reset_entry(self.txt_dwi, self.parameter_values['dwi'])
            self.reset_entry(self.txt_bval, filename_wo_ext(self.parameter_values['dwi'])+'.bval')
            self.reset_entry(self.txt_bvec, filename_wo_ext(self.parameter_values['dwi'])+'.bvec')
        if self.parameter_values.has_key('b0'):
            self.reset_entry(self.txt_b0, self.parameter_values['b0'])
        if self.parameter_values.has_key('working'):
            self.reset_entry(self.txt_working, self.parameter_values['working'], force=True)

    def reset_entry(self, entry, text, force=False):
        entry.delete(0, len(entry.get()))

        dirname = self.txt_working.get()
        len_dirname = len(dirname)
        # print dirname, text
        if (not force) and type(text) == type('') and len_dirname > 1 and dirname == text[:len_dirname]:
            entry.insert(0, text[len_dirname:].strip(os.path.sep))
        else:
            entry.insert(0, text)

    def filenameDialog_text(self, text):
        fl = tkFileDialog.askopenfilename()
        if fl != '':
            self.reset_entry(text, fl)

    def dirnameDialog_text(self, text):
        fl = tkFileDialog.askdirectory(initialdir=text.get())
        if fl != '':
            self.reset_entry(text, fl)

    def prefix(self):
        if self.txt_subject.get() == '':
            return ''
        else:
            return self.txt_subject.get() + '_'

    def update_param_from_text(self):
        self.param.subject              = self.txt_subject.get()
        self.param.working_dir          = self.txt_working.get()
        self.param.fn_bval              = self.txt_bval.get()
        self.param.fn_bvec              = self.txt_bvec.get()
        self.param.b0_threshold         = self.txt_b0_threshold.get()

        self.param.topup_app            = self.txt_topup_app.get()
        self.param.topup_config         = self.txt_topup_config.get()
        self.param.eddy_app             = self.txt_eddy_app.get()

        self.param.fn_b0                = self.txt_b0.get()
        self.param.fn_b0_mask           = self.txt_b0mask.get()
        self.param.fn_dwi               = self.txt_dwi.get()
        self.param.fn_dwi_mask          = self.txt_dwi_mask.get()
        self.param.fn_dwi_eddy_out      = self.txt_output_eddy.get()

        self.param.multiband            = self.spin_multiband.get()
        self.param.niter                = self.spin_nitr.get()
        self.param.nosearch             = self.nosearch.get()
        self.param.noresample           = self.noresample.get()
        self.param.schedule             = self.txt_reg_config.get()

        self.param.fn_dwi_eddy          = self.txt_dwi_eddy.get()
        self.param.fn_reg_mask          = self.txt_reg_mask.get()
        self.param.fn_reg_outlier       = self.txt_reg_outlier.get()
        self.param.fn_reg_out           = self.txt_output_reg.get()
        self.param.fn_dti               = self.txt_dti.get()
        self.param.fn_dwi_dti           = self.txt_dwi_dti.get()

        self.param.fn_dwi_dti_mask      = self.txt_md_mask.get()
        self.param.fn_dwi_dti_outlier   = self.txt_dti_outlier.get()

        self.param.fn_dti_value         = self.txt_dti_value.get()
        self.param.fn_dti_value_roi     = self.txt_dti_roi.get()
        self.param.fn_dti_value_out     = self.txt_output_dti_value.get()

        self.param.set_variable_types()

    def update_text_from_param(self):
        self.reset_entry(self.txt_subject         , self.param.subject)
        self.reset_entry(self.txt_working         , self.param.working_dir, force=True)
        self.reset_entry(self.txt_bval            , self.param.fn_bval)
        self.reset_entry(self.txt_bvec            , self.param.fn_bvec)
        self.reset_entry(self.txt_b0_threshold    , self.param.b0_threshold)
                                                        
        self.reset_entry(self.txt_topup_app       , self.param.topup_app)
        self.reset_entry(self.txt_topup_config    , self.param.topup_config)
        self.reset_entry(self.txt_eddy_app        , self.param.eddy_app)
                                                        
        self.reset_entry(self.txt_b0              , self.param.fn_b0)
        self.reset_entry(self.txt_b0mask          , self.param.fn_b0_mask)
        self.reset_entry(self.txt_dwi             , self.param.fn_dwi)
        self.reset_entry(self.txt_dwi_mask        , self.param.fn_dwi_mask)
        self.reset_entry(self.txt_output_eddy     , self.param.fn_dwi_eddy_out)
                                                        
        self.reset_entry(self.spin_multiband      , self.param.multiband)
        self.reset_entry(self.spin_nitr           , self.param.niter)
        self.nosearch.set(self.param.nosearch)
        self.noresample.set(self.param.noresample)
        self.reset_entry(self.txt_reg_config      , self.param.schedule)
                                                        
        self.reset_entry(self.txt_dwi_eddy        , self.param.fn_dwi_eddy)
        self.reset_entry(self.txt_reg_mask        , self.param.fn_reg_mask)
        self.reset_entry(self.txt_reg_outlier     , self.param.fn_reg_outlier)
        self.reset_entry(self.txt_output_reg      , self.param.fn_reg_out)
                                                        
        self.reset_entry(self.txt_dwi_dti         , self.param.fn_dwi_dti)
        self.reset_entry(self.txt_dti             , self.param.fn_dti)
        self.reset_entry(self.txt_md_mask         , self.param.fn_dwi_dti_mask)
        self.reset_entry(self.txt_dti_outlier     , self.param.fn_dwi_dti_outlier)

        self.reset_entry(self.txt_dti_value       , self.param.fn_dti_value)
        self.reset_entry(self.txt_dti_roi         , self.param.fn_dti_value_roi)
        self.reset_entry(self.txt_output_dti_value, self.param.fn_dti_value_out)
       

    def read_param(self, filename=None):
        if filename is None or not os.path.isfile(filename):
            fn = tkFileDialog.askopenfilename()
        else:
            fn = filename
        if os.path.isfile(fn):
            print 'Reading %s' % fn
            self.param.read(fn)
            self.update_text_from_param()
            self.reset_entry(self.txt_param, fn)

        else:
            sys.stderr.write('File not exist: %s\n' % fn)

    def save_param(self):
        self.update_param_from_text()
        if self.txt_param.get() != '':
            if os.path.dirname(self.txt_param.get()) == '':
                filename = tkFileDialog.asksaveasfilename(initialdir=self.txt_working.get(), initialfile=os.path.basename(self.txt_param.get()))
            else:
                filename = tkFileDialog.asksaveasfilename(initialdir=os.path.dirname(self.txt_param.get()), initialfile=os.path.basename(self.txt_param.get()))
        else:
            filename = tkFileDialog.asksaveasfilename(initialdir=self.txt_working.get())

        if filename != '':
            self.param.save(filename)
            self.reset_entry(self.txt_param, filename)

    def initUI(self):
        self.parent.title(self.title)

        self.nosearch = BooleanVar(self)
        self.noresample = BooleanVar(self)

        i = 0
        #frame_top = Frame(self)
        #frame_top.grid(row=i, column=0, columnspan=6, sticky=EW)
        #Label(frame_top, text='topframe').grid(row=0, column=0, sticky=EW)
        #Label(frame_top, text='topframe').grid(row=0, column=1, sticky=EW)
        #Label(frame_top, text='topframe').grid(row=1, column=0, sticky=EW)
        #Label(frame_top, text='topframe').grid(row=1, column=1, sticky=EW)
        #Entry(frame_top, text='topframe').grid(row=0, column=0, sticky=EW)
        #Entry(frame_top, text='topframe').grid(row=0, column=1, sticky=EW)
        #Entry(frame_top, text='topframe').grid(row=1, column=0, sticky=EW)
        #Entry(frame_top, text='topframe').grid(row=1, column=1, sticky=EW)
        #frame_top.grid_columnconfigure(0, weight=1)
        #frame_top.grid_columnconfigure(1, weight=1)

        dddWidth = 4
        labelWidth = 10
        i += 1;  Label(self, text='Working dir'  , width=labelWidth).grid(row=i, column=0)
        self.txt_working = Entry(self); self.txt_working.grid(row=i, column=1, sticky=EW)
        btn_working = Button(self, text='...', width=dddWidth, command=lambda:self.dirnameDialog_text(self.txt_working))
        btn_working.grid(row=i, column=2, sticky=W)
        i += 1;  Label(self, text='Parameters', width=labelWidth).grid(row=i, column=0)
        self.txt_param = Entry(self); self.txt_param.grid(row=i, column=1, sticky=EW)
        btn_param_read = Button(self, text='Read', command=self.read_param); btn_param_read.grid(row=i, column=2)
        btn_param_save = Button(self, text='Save', command=self.save_param); btn_param_save.grid(row=i, column=3)

        i += 1
        frm_prepa = Frame(self)
        ii  = 0; Label(frm_prepa, text='Subject'      , width=labelWidth).grid(row=ii, column=0)
        self.txt_subject = Entry(frm_prepa); self.txt_subject.grid(row=ii, column=1, sticky=EW)
        btn_prepa = Button(frm_prepa, text='Prepare DWIs', command=self.run_prepare); btn_prepa.grid(row=ii, column=3)
        #btn_update = Button(frm_prepa, text='Update filenames', command=self.update_parameter); btn_update.grid(row=ii, column=4)
        btn_update = Button(frm_prepa, text='Update filenames', command=self.update_text_from_param); btn_update.grid(row=ii, column=4)

        ii += 1; Label(frm_prepa, text='bval'         , width=labelWidth).grid(row=ii, column=0)
        self.txt_bval = Entry(frm_prepa); self.txt_bval.grid(row=ii, column=1, sticky=EW)
        btn_bval = Button(frm_prepa, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_bval));
        btn_bval.grid(row=ii, column=2, sticky=W)

        Label(frm_prepa, text='B0 threshold').grid(row=ii, column=3, sticky=E)
        self.txt_b0_threshold = Entry(frm_prepa, width=10); self.txt_b0_threshold.grid(row=ii, column=4)

        ii += 1; Label(frm_prepa, text='bvec'         , width=labelWidth).grid(row=ii, column=0)
        self.txt_bvec = Entry(frm_prepa); self.txt_bvec.grid(row=ii, column=1, sticky=EW)
        btn_bvec = Button(frm_prepa, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_bvec));
        btn_bvec.grid(row=ii, column=2, sticky=W)

        frm_prepa.grid_columnconfigure(1, weight=1)
        frm_prepa.grid(row=i, rowspan=ii+1, column=0, columnspan=6, sticky=NSEW)
        ii += 1

        i += ii
        # TOPUP/EDDY or Rigid registration
        i += 1; Label(self, text='  TOPUP/EDDY or Rigid registration').grid(row=i, column=0, columnspan=2, sticky=W)
        i += 1
        frm_eddy = Frame(self)
        width = 10
        ii  = 0
        jj  = 0; Label(frm_eddy, text='TOPUP', width=width, justify=RIGHT, anchor=E).grid(row=ii, column=jj, sticky=E)
        jj += 1; self.txt_topup_app = Entry(frm_eddy, width=width); self.txt_topup_app.grid(row=ii, column=jj)
        jj += 1; Label(frm_eddy, text='Configure', width=width, justify=RIGHT, anchor=E).grid(row=ii, column=jj, sticky=E)
        jj += 1; self.txt_topup_config = Entry(frm_eddy); self.txt_topup_config.grid(row=ii, column=jj, sticky=EW)
        jj += 1; btn_topup_config = Button(frm_eddy, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_topup_config))
        btn_topup_config.grid(row=ii, column=jj, sticky=W)
        jj += 1; Label(frm_eddy, text='EDDY', width=width, justify=RIGHT, anchor=E).grid(row=ii, column=jj, sticky=E)
        jj += 1; self.txt_eddy_app = Entry(frm_eddy, width=width); self.txt_eddy_app.grid(row=ii, column=jj)

        frm_eddy.grid(row=i, rowspan=ii+1, column=0, columnspan=6, sticky=NSEW)
        frm_eddy.grid_columnconfigure(3, weight=1)

        i += ii

        i += 1; Label(self, text='B0'           , width=labelWidth).grid(row=i, column=0)
        self.txt_b0 = Entry(self); self.txt_b0.grid(row=i, column=1, sticky=EW)
        btn_b0 = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_b0))
        btn_b0.grid(row=i, column=2, sticky=W)
        btn_topup = Button(self, text='TOPUP', command=self.run_topup); btn_topup.grid(row=i, column=4, rowspan=2, sticky=NS)
        btn_rigid = Button(self, text='Rigid\nReg', command=self.run_eddy); btn_rigid.grid(row=i, column=5, rowspan=4, sticky=NS)

        i += 1; Label(self, text='B0 mask'      , width=labelWidth).grid(row=i, column=0)
        self.txt_b0mask = Entry(self); self.txt_b0mask.grid(row=i, column=1, sticky=EW)
        btn_b0mask = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_b0mask))
        btn_b0mask.grid(row=i, column=2, sticky=W)
        btn_b0mask_make = Button(self, text='Make', command=lambda:self.make_mask(self.txt_b0.get(), mask4d=True, rtn=self.txt_b0mask))
        btn_b0mask_make.grid(row=i, column=3)

        i += 1; Label(self, text='DWI'          , width=labelWidth).grid(row=i, column=0)
        self.txt_dwi = Entry(self); self.txt_dwi.grid(row=i, column=1, sticky=EW)
        btn_dwi_eddy = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_dwi))
        btn_dwi_eddy.grid(row=i, column=2, sticky=W)
        btn_eddy = Button(self, text='EDDY', command=self.run_eddy); btn_eddy.grid(row=i, column=4, rowspan=2, sticky=NS)

        i += 1; Label(self, text='DWI mask'     , width=labelWidth).grid(row=i, column=0)
        self.txt_dwi_mask = Entry(self); self.txt_dwi_mask.grid(row=i, column=1, sticky=EW)
        btn_dwi_mask = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_dwi_mask))
        btn_dwi_mask.grid(row=i, column=2, sticky=W)
        btn_dwi_mask_make = Button(self, text='Make', command=lambda:self.make_mask(self.txt_dwi.get(), rtn=self.txt_dwi_mask))
        btn_dwi_mask_make.grid(row=i, column=3)

        i += 1; Label(self, text='Output'       , width=labelWidth).grid(row=i, column=0)
        self.txt_output_eddy = Entry(self); self.txt_output_eddy.grid(row=i, column=1, sticky=EW)
        btn_output_eddy = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_output_eddy))
        btn_output_eddy.grid(row=i, column=2, sticky=W)
        btn_output_eddy_auto = Button(self, text='Auto', command=self.gen_eddy_outputname); btn_output_eddy_auto.grid(row=i, column=3)


        # 2D registration
        i += 1; Label(self, text='  2D registration (X-Y translation)').grid(row=i, column=0, columnspan=2, sticky=W)

        i += 1
        frm_param = Frame(self)
        width = 12
        width_half = 6
        ii  = 0
        jj  = 0; Label(frm_param, text='Multi band'   , width=width, justify=RIGHT, anchor=E).grid(row=ii, column=jj, sticky=E)
        #jj += 1; self.txt_multiband = Entry(frm_param, width=width_half); self.txt_multiband.grid(row=ii, column=jj)
        jj += 1; self.spin_multiband = Spinbox(frm_param, from_=0, to=1000, increment=1, width=width_half); self.spin_multiband.grid(row=ii, column=jj)
        jj += 1; Label(frm_param, text='#iteration'   , width=width, justify=RIGHT, anchor=E).grid(row=ii, column=jj, sticky=E)
        #jj += 1; self.txt_nitr = Entry(frm_param, width=width_half); self.txt_nitr.grid(row=ii, column=jj)
        jj += 1; self.spin_nitr = Spinbox(frm_param, from_=0, to=1000, increment=1, width=width_half); self.spin_nitr.grid(row=ii, column=jj)
        jj += 1; Label(frm_param, text=''   , width=width, justify=RIGHT, anchor=E).grid(row=ii, column=jj, sticky=E)
        jj += 1; self.chk_nosearch = Checkbutton(frm_param, text='nosearch', variable=self.nosearch, width=width)
        self.chk_nosearch.grid(row=ii, column=jj)
        jj += 1; self.chk_noresample = Checkbutton(frm_param, text='noresample', variable=self.noresample, width=width)
        self.chk_noresample.grid(row=ii, column=jj)
        ii += 1
        jj  = 0; Label(frm_param, text='Schedule').grid(row=ii, column=jj, sticky=E)
        jj += 1; self.txt_reg_config = Entry(frm_param); self.txt_reg_config.grid(row=ii, column=jj, sticky=EW, columnspan=4)
        jj += 4; btn_reg_config = Button(frm_param, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_reg_config))
        btn_reg_config.grid(row=ii, column=jj, sticky=W)

        frm_param.grid_columnconfigure(4, weight=1)
        frm_param.grid(row=i, rowspan=ii+1, column=0, columnspan=6, sticky=NSEW)

        i += ii

        i += 1; Label(self, text='DWI'          , width=labelWidth).grid(row=i, column=0)
        self.txt_dwi_eddy = Entry(self); self.txt_dwi_eddy.grid(row=i, column=1, sticky=EW)
        btn_dwi_eddy = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_dwi_eddy))
        btn_dwi_eddy.grid(row=i, column=2, sticky=W)
        btn_dwi_eddy_copy = Button(self, text='Copy', command=lambda:self.reset_entry(self.txt_dwi_eddy, self.txt_output_eddy.get()))
        btn_dwi_eddy_copy.grid(row=i, column=3)
        btn_dwi_eddy = Button(self, text='Save\nParam', command=self.save_reg_param); btn_dwi_eddy.grid(row=i, column=4, rowspan=2, sticky=NSEW)
        btn_dwi_eddy = Button(self, text='XY-Reg', command=self.run_xy_reg); btn_dwi_eddy.grid(row=i, column=5, rowspan=2, sticky=NSEW)

        i += 1; Label(self, text='Mask'         , width=labelWidth).grid(row=i, column=0)
        self.txt_reg_mask = Entry(self); self.txt_reg_mask.grid(row=i, column=1, sticky=EW)
        btn_reg_mask = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_reg_mask))
        btn_reg_mask.grid(row=i, column=2, sticky=W)
        btn_reg_mask_make = Button(self, text='Make', command=lambda:self.make_mask(self.txt_dwi_eddy.get(), rtn=self.txt_reg_mask)); btn_reg_mask_make.grid(row=i, column=3)

        i += 1; Label(self, text='Outlier'      , width=labelWidth).grid(row=i, column=0)
        self.txt_reg_outlier = Entry(self); self.txt_reg_outlier.grid(row=i, column=1, sticky=EW)
        btn_reg_outlier = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_reg_outlier))
        btn_reg_outlier.grid(row=i, column=2, sticky=W)
        btn_reg_outlier_calculate = Button(self, text='Calculate', command=self.run_xy_reg_outlier)
        btn_reg_outlier_calculate.grid(row=i, column=3)
        btn_reg_apply = Button(self, text='Applywarp', command=self.run_applywarp); btn_reg_apply.grid(row=i, column=4, columnspan=2, sticky=EW)

        i += 1; Label(self, text='Output'       , width=labelWidth).grid(row=i, column=0)
        self.txt_output_reg = Entry(self); self.txt_output_reg.grid(row=i, column=1, sticky=EW)
        btn_output_reg = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_output_reg))
        btn_output_reg.grid(row=i, column=2, sticky=W)
        btn_output_reg_auto = Button(self, text='Auto', command=self.gen_reg_outputname); btn_output_reg_auto.grid(row=i, column=3)


        # DTI
        i += 1; Label(self, text='  DTI map with outlier rejection').grid(row=i, column=0, columnspan=2, sticky=W)

        i += 1; Label(self, text='DWI'          , width=labelWidth).grid(row=i, column=0)
        self.txt_dwi_dti = Entry(self); self.txt_dwi_dti.grid(row=i, column=1, sticky=EW)
        btn_dwi_dti = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_dwi_dti))
        btn_dwi_dti.grid(row=i, column=2, sticky=W)
        btn_dwi_dti_copy = Button(self, text='Copy', command=lambda:self.reset_entry(self.txt_dwi_dti, self.txt_output_reg.get()))
        btn_dwi_dti_copy.grid(row=i, column=3)
        btn_dti_maps = Button(self, text='Generate DTI maps', command=self.run_generate_dti_maps); btn_dti_maps.grid(row=i, column=4, columnspan=2, sticky=EW)

        i += 1; Label(self, text='DTI'          , width=labelWidth).grid(row=i, column=0)
        self.txt_dti = Entry(self); self.txt_dti.grid(row=i, column=1, sticky=EW)
        btn_dti = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_dti))
        btn_dti.grid(row=i, column=2, sticky=W)
        btn_dti_auto = Button(self, text='Auto', command=self.gen_dti_name); btn_dti_auto.grid(row=i, column=3)
        btn_dti_outlier_calc = Button(self, text='Calculate Outlier', command=self.run_dti_outlier); btn_dti_outlier_calc.grid(row=i, column=4, columnspan=2, rowspan=2, sticky=NSEW)

        i += 1; Label(self, text='Mask'         , width=labelWidth).grid(row=i, column=0)
        self.txt_md_mask = Entry(self); self.txt_md_mask.grid(row=i, column=1, sticky=EW)
        btn_md_mask = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_md_mask))
        btn_md_mask.grid(row=i, column=2, sticky=W)
        btn_md_mask_make = Button(self, text='Make', command=self.run_md_mask); btn_md_mask_make.grid(row=i, column=3)

        i += 1; Label(self, text='Outlier'      , width=labelWidth).grid(row=i, column=0)
        self.txt_dti_outlier = Entry(self); self.txt_dti_outlier.grid(row=i, column=1, sticky=EW)
        btn_dti_outlier = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_dti_outlier))
        btn_dti_outlier.grid(row=i, column=2, sticky=W)
        btn_dti_maps_itr = Button(self, text='Generate DTI maps', command=lambda:self.run_generate_dti_maps(outlier=True)); btn_dti_maps_itr.grid(row=i, column=4, columnspan=2, sticky=EW)

        #i += 1
        #Label(self, text='DWI'          , width=labelWidth).grid(row=i, column=0)
        #self.txt_dwi_dti_itr = Entry(self); self.txt_dwi_dti_itr.grid(row=i, column=1, sticky=EW)
        #btn_dwi_dti_itr = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_dwi_dti_itr))
        #btn_dwi_dti_itr.grid(row=i, column=2, sticky=W)

        #i += 1; Label(self, text='Output'       , width=labelWidth).grid(row=i, column=0)
        #self.txt_output_dti = Entry(self); self.txt_output_dti.grid(row=i, column=1, sticky=EW)
        #btn_output_dti = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_output_dti))
        #btn_output_dti.grid(row=i, column=2, sticky=W)
        #btn_output_dti_auto = Button(self, text='Auto', command=self.gen_dti_outputname); btn_output_dti_auto.grid(row=i, column=3)


        # DTI Values
        i += 1; Label(self, text='  DTI values').grid(row=i, column=0, columnspan=2, sticky=W)

        i += 1; Label(self, text='DTI'          , width=labelWidth).grid(row=i, column=0)
        self.txt_dti_value = Entry(self); self.txt_dti_value.grid(row=i, column=1, sticky=EW)
        btn_dti_value = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_dti_value))
        btn_dti_value.grid(row=i, column=2, sticky=W)
        btn_dwi_dti_value_copy = Button(self, text='Copy', command=lambda:self.reset_entry(self.txt_dti_value, self.txt_dti.get()))
        btn_dwi_dti_value_copy.grid(row=i, column=3)
        btn_dti_value_get = Button(self, text='Get DTI\nValues', command=self.run_get_dti_values); btn_dti_value_get.grid(row=i, column=4, rowspan=2, columnspan=2, sticky=NSEW)

        i += 1; Label(self, text='ROI'          , width=labelWidth).grid(row=i, column=0)
        self.txt_dti_roi = Entry(self); self.txt_dti_roi.grid(row=i, column=1, sticky=EW)
        btn_dti_roi = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_dti_roi))
        btn_dti_roi.grid(row=i, column=2, sticky=W)
        btn_dti_roi = Button(self, text='Make'); btn_dti_roi.grid(row=i, column=3)

        i += 1; Label(self, text='Output'       , width=labelWidth).grid(row=i, column=0)
        self.txt_output_dti_value = Entry(self); self.txt_output_dti_value.grid(row=i, column=1, sticky=EW)
        btn_output_dti_value = Button(self, text='...', width=dddWidth, command=lambda:self.filenameDialog_text(self.txt_output_dti_value))
        btn_output_dti_value.grid(row=i, column=2, sticky=W)
        btn_output_dti_value_auto = Button(self, text='Auto', command=self.gen_dti_value_outputname); btn_output_dti_value_auto.grid(row=i, column=3)


        self.grid_columnconfigure(1, weight=1)
        ni = i + 1
        for i in range(ni):
            self.grid_rowconfigure(i, weight=1, minsize=20)

        self.pack(fill=BOTH, expand=True)

    def initUI_1(self):

        self.parent.title("Review")
        self.pack(fill=BOTH, expand=True)

        frame1 = Frame(self)
        frame1.pack(fill=X)

        lbl1 = Label(frame1, text="Title", width=6)
        lbl1.pack(side=LEFT, padx=5, pady=5)

        entry1 = Entry(frame1)
        entry1.pack(fill=X, padx=5, expand=True)

        frame2 = Frame(self)
        frame2.pack(fill=X)
        #frame2.pack(fill=BOTH, expand=True)

        lbl2 = Label(frame2, text="Author", width=6)
        lbl2.pack(side=LEFT, padx=5, pady=5)

        entry2 = Entry(frame2)
        entry2.pack(fill=X, padx=5, expand=True)

        frame3 = Frame(self)
        frame3.pack(fill=BOTH, expand=True)

        lbl3 = Label(frame3, text="Review", width=6)
        lbl3.pack(side=LEFT, anchor=N, padx=5, pady=5)

        txt = Entry(frame3)
        txt.pack(fill=BOTH, pady=5, padx=5, expand=True)


    def onOpen(self):
        ftypes = [('Python files', '*.py'), ('All files', '*')]
        dlg = tkFileDialog.Open(self, filetypes=ftypes)
        fl = dlg.show()

        if fl != '':
            #text = self.readFile(fl)
            print fl

    def gen_eddy_outputname(self):
        fn = self.txt_output_eddy.get()
        if fn == '' or '_rigid.nii.gz' in fn:
            fn = filename_wo_ext(self.txt_dwi.get()) + '_eddy_unwarped.nii.gz'
        else:
            fn = filename_wo_ext(self.txt_dwi.get()) + '_rigid.nii.gz'
        self.reset_entry(self.txt_output_eddy, fn)

    def gen_reg_outputname(self):
        fn = filename_wo_ext(self.txt_dwi_eddy.get()) + '_xenc.nii.gz'
        self.reset_entry(self.txt_output_reg, fn)

    #def gen_dti_outputname(self):
    #    fn = filename_wo_ext(self.txt_output_reg.get()) + '_dti.nii.gz'
    #    self.reset_entry(self.txt_output_dti, fn)

    def gen_dti_name(self):
        self.reset_entry(self.txt_dti, filename_wo_ext(self.txt_dwi_dti.get()) + '_dti')

    def gen_dti_value_outputname(self):
        fn = filename_wo_ext(self.txt_dti_value.get()) + '.csv'
        self.reset_entry(self.txt_output_dti_value, fn)

    def check_dependency(self, dependency):
        dirname = self.txt_working.get()
        for entry in dependency:
            if not os.path.isfile(os.path.join(dirname, entry.get())):
                sys.stderr.write('File not exist: %s\n' % entry.get())
                return False
        return True

    def run_prepare(self):
        #prep = Tk()
        #app = prepare_dwi_gui.PrepareDWI(prep, subject=self.txt_subject.get(), output_dir=self.txt_working.get(), return_value=self.parameter_values, make_param=False)
        #prep.mainloop()
        self.prep_window = Toplevel(self.parent)
        self.prep_app = prepare_dwi_gui.PrepareDWI(self.prep_window, subject=self.txt_subject.get(), output_dir=self.txt_working.get(), return_value=self.parameter_values, obj_return_value=self, make_param=False)
        #app.mainloop()
        #self.update_parameter()

        #tkMessageBox.showinfo(self.title, 'Done: DWI preparation')

    def make_mask(self, filename=None, mask4d=False, rtn=None):
        import gui_mask
        root = Tk()
        #current_pwd = os.path.abspath(os.path.curdir())
        working = self.txt_working.get()
        #os.chdir(working)
        if filename is None or filename == '':
            app = gui_mask.CreateMask(root, filename=None, dirname=working, mask4d=mask4d, obj_return_value=rtn)
        else:
            app = gui_mask.CreateMask(root, filename=os.path.join(working, filename), dirname=working, mask4d=mask4d, obj_return_value=rtn)
        root.mainloop()
        #os.chdir(current_pwd)

        #tkMessageBox.showinfo(self.title, 'Done: Making mask')

    def run_topup(self):
        dependency = [self.txt_b0, self.txt_b0mask, ]
        if not self.check_dependency(dependency):
            return

        working = self.txt_working.get()
        fn_b0 = os.path.join(working, self.txt_b0.get())
        fn_b0mask = os.path.join(working, self.txt_b0mask.get())
        fn_b0masked = filename_wo_ext(fn_b0) + '_masked.nii.gz'
        fn_out = os.path.join(working, self.prefix() + 'topup')
        run_command('fslmaths %s -mas %s %s' % (fn_b0, fn_b0mask, fn_b0masked))

        topup_command = [self.txt_topup_app.get()]
        if self.txt_topup_config.get() != '':
            topup_command.append('--config=%s' % self.txt_topup_config.get())
        topup_command.append('--datain=%s' % (os.path.join(working, self.prefix() + 'acqparams.txt')))
        topup_command.append('--imain=%s' % fn_b0masked)
        topup_command.append('--out=%s' % fn_out)
        topup_command.append('--iout=%s_warped' % fn_out)
        topup_command.append('--fout=%s_field' % fn_out)
        topup_command.append('--logout=%s_log' % fn_out)
        topup_command.append('1> %s_topup_log.out' % self.txt_subject.get())

        cmd = ' '.join(topup_command)
        append_text(self.log, cmd, 'run_topup')
        run_command(cmd)
        print 'Done: TOPUP'

        tkMessageBox.showinfo(self.title, 'Done: TOPUP')

    def run_eddy(self):
        dependency = [self.txt_dwi, self.txt_dwi_mask, ]
        if not self.check_dependency(dependency):
            return

        working = self.txt_working.get()
        eddy_command = ['eddy']
        eddy_command.append('--flm=quadratic')
        eddy_command.append('--acqp=%s' % (os.path.join(working, self.prefix() + 'acqparams.txt')))
        eddy_command.append('--bvals=%s' % (os.path.join(working, self.txt_bval.get())))
        eddy_command.append('--bvecs=%s' % (os.path.join(working, self.txt_bvec.get())))
        eddy_command.append('--imain=%s' % (os.path.join(working, self.txt_dwi.get())))
        eddy_command.append('--index=%s' % (os.path.join(working, self.prefix() + 'index.txt')))
        eddy_command.append('--mask=%s' % (os.path.join(working, self.txt_dwi_mask.get())))
        eddy_command.append('--topup=%s' % (os.path.join(working, self.prefix() + 'topup')))
        eddy_command.append('1> %s_eddy_log.out' % self.txt_subject.get())

        if self.txt_output_eddy.get() == '':
            self.gen_eddy_outputname()

        eddy_command.append('--out=%s' % (os.path.join(working, self.txt_output_eddy.get())))

        cmd = ' '.join(eddy_command)
        append_text(self.log, cmd, 'run_eddy')
        run_command(cmd)
        print 'Done: Eddy'

        self.reset_entry(self.txt_dwi_eddy, self.txt_output_eddy.get())
        tkMessageBox.showinfo(self.title, 'Done: EDDY')

    def run_rigid_reg(self):
        dependency = [self.txt_b0, self.txt_dwi_mask,  ]
        if not self.check_dependency(dependency):
            return
        tkMessageBox.showinfo(self.title, 'Non implemented yet')
        #tkMessageBox.showinfo(self.title, 'Done: Rigid registration')

    def save_reg_param(self):
        dependency = [self.txt_dwi_eddy, self.txt_reg_mask, ]
        if not self.check_dependency(dependency):
            return

        filename = tkFileDialog.asksaveasfilename(initialdir=self.txt_working.get(), initialfile='reg_2d.params')
        if filename == '':
            return

        self.update_param_from_text()
        param_2d = parameter.Parameter_reg2d(param_object=self.param)
        param_2d.save_param(filename)

    def run_xy_reg(self):
        dependency = [self.txt_dwi_eddy, self.txt_reg_mask ]
        if not self.check_dependency(dependency):
            return

        self.update_param_from_text()
        param_2d = parameter.Parameter_reg2d(param_object=self.param)

        spine_reg.set_print_cmd(os.path.join(self.txt_working.get(), self.prefix() + 'xy_reg_cmd'))

        cmd = '\n'.join([str(vars(param_2d)), 'spine_reg.run_registration(param_2d)'])
        append_text(self.log, cmd, 'run_xy_reg')

        fn_out = spine_reg.run_registration(param_2d)
        self.reset_entry(self.txt_dwi_dti, fn_out)
        os.chdir(self.txt_working.get())
        print 'Done: XY-Reg'
        tkMessageBox.showinfo(self.title, 'Done: XY-registration')

    def run_applywarp(self):
        dependency = [self.txt_dwi_eddy, self.txt_reg_mask, ]
        if not self.check_dependency(dependency):
            return

        param_2d = parameter.Parameter_reg2d(param_object=self.param)
        spine_reg.set_print_cmd(os.path.join(self.txt_working.get(), self.prefix() + 'xy_reg_cmd'))
        fn_out = spine_reg.run_applywarp(param_2d)
        self.reset_entry(self.txt_dwi_dti, fn_out)
        print 'Done: Applywarp'
        tkMessageBox.showinfo(self.title, 'Done: Applywarp')

    def run_xy_reg_outlier(self, param_2d=None):
        dependency = [self.txt_dwi_eddy, self.txt_reg_mask, self.txt_bval]
        #if not self.check_dependency(dependency):
        #    return

        self.update_param_from_text()

        if param_2d is None:
            param_2d = parameter.Parameter_reg2d(param_object=self.param)
        spine_reg.generate_xy_trans(param_2d)
        os.chdir(self.txt_working.get())

        print 'Done: XY-Reg outlier calculation'

    def run_generate_dti_maps(self, outlier=False):
        cmd = 'generate_dti_maps(%s)' % (', '.join([self.txt_dwi_dti.get(),
                                                    'bval_bvec=None',
                                                    'bval=%s' % os.path.join(self.txt_working.get(), self.txt_bval.get()),
                                                    'bvec=%s' % os.path.join(self.txt_working.get(), self.txt_bvec.get()),
                                                    'prefix=%s' %self.prefix(),
                                                    'outlier=%s' % outlier]))
        append_text(self.log, cmd, 'run_generate_dti_maps')
        fn_out = generate_dti_maps(self.txt_dwi_dti.get(), bval_bvec=None, bval=os.path.join(self.txt_working.get(), self.txt_bval.get()), bvec=os.path.join(self.txt_working.get(), self.txt_bvec.get()), prefix=self.prefix(), outlier=outlier)
        print 'Done: Generate DTI maps'

        self.reset_entry(self.txt_dti_value, fn_out)
        tkMessageBox.showinfo(self.title, 'Done: Generating DTI maps')

    def run_md_mask(self):
        if not os.path.isfile(os.path.join(self.txt_working.get(), self.txt_dwi_dti.get())):
            sys.stderr.write('File not exist: %s\n' % self.txt_dwi_dti.get())
            return
        working = self.txt_working.get()

        if os.path.isfile(os.path.join(working, self.txt_dti_value.get())):
            fn_md = filename_wo_ext(self.txt_dti_value.get()) + '_MD.nii.gz'
        else:
            fn_md = filename_wo_ext(self.txt_dwi_dti.get()) + '_dti_MD.nii.gz'

        fn_mask = create_md_mask(os.path.join(working, fn_md), fn_out=None, thr_min=0.1, thr_max=1.1)
        print 'Done: Making MD mask'

        self.reset_entry(self.txt_md_mask, os.path.basename(fn_mask))
        #tkMessageBox.showinfo(self.title, 'Done: MD mask')

    def run_dti_outlier(self):
        import outlier_4dfp_gui
        root = Tk()
        working = self.txt_working.get()

        filename     = os.path.join(working, filename_wo_ext(self.txt_dwi_dti.get()) + '_res.nii.gz')
        filename_roi = os.path.join(working, self.txt_md_mask.get())
        filename_csv = os.path.join(working, self.txt_dti_outlier.get())

        if filename_csv == '':
            filename_csv = filename_wo_ext(filename) + '.csv'

        print filename, filename_roi, filename_csv

        app = outlier_4dfp_gui.Outlier4dfp(root, filename=filename, filename_roi=filename_roi, filename_csv=filename_csv, dirname=working, obj_return_value=self.txt_dti_outlier)
        root.mainloop()

        #tkMessageBox.showinfo(self.title, 'Done: Making mask')

    def run_get_dti_values(self):

        filenames = [ '%s_%s' % (self.txt_dti_value.get(), tmp) for tmp in ['MD', 'FA', 'RD', 'RA', 'Eig1', 'Eig2', 'Eig3']]
        filename_roi = self.txt_dti_roi.get()
        filename_csv = self.txt_output_dti_value.get()
        if filename_csv == '':
            self.gen_dti_value_outputname()
            filename_csv = self.txt_output_dti_value.get()

        value_matrix = calculate_voxel_values(filenames, filename_roi, filename_csv=filename_csv)

        print 'Done: Getting DTI values'
        #tkMessageBox.showinfo(self.title, 'Done: DTI values')

def main():

    root = Tk()

    filename = None
    if len(sys.argv) > 1:
        if os.path.isfile(sys.argv[1]):
            filename = sys.argv[1]

    app = SpineRegGui(root, param_filename=filename)
    root.mainloop()


if __name__ == '__main__':
    main()

