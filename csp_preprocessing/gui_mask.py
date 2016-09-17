#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')

import numpy as np

import scipy.ndimage as ndimage

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigureCanvas
from matplotlib.figure import Figure

from Tkinter import Tk, Entry, TOP, BOTH, RAISED, N, E, W, S, EW, NS, NW, NSEW
from Tkinter import RIDGE, TRUE, FALSE, VERTICAL, LEFT, RIGHT, X, Y, BOTTOM
from Tkinter import Grid, BooleanVar, StringVar, Canvas, Scrollbar, Spinbox
from ttk import Frame, Button, Style, Label, Entry, Checkbutton, Radiobutton
import tkFileDialog

import create_mask
import nibabel as nib
import functools
import parameter
import os

class MplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=4, height=2, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure()

        #
        FigureCanvas.__init__(self, fig, master=parent)
        self.show()
        #self.get_tk_widget().pack(side=TOP, fill=BOTH, expand=TRUE)

        #FigureCanvas.setSizePolicy(self,
        #                           QtGui.QSizePolicy.Expanding,
        #                           QtGui.QSizePolicy.Expanding)
        #FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass

class ImshowMplCanvas(MplCanvas):
    """Simple canvas with a sine plot."""
    def drawData(self, data, mask, drawMask):
        shape = (data.shape[1], data.shape[0], 3)
        dataRGB = np.ndarray(shape, dtype=np.float16)
        dataRGB[:,:,0] = dataRGB[:,:,1] = dataRGB[:,:,2] = data.T.copy()
        dataRGB /= dataRGB.max()

        if drawMask:
            seed = (mask == 1).nonzero()
            ijs = zip(seed[0], seed[1])
            for (j,i) in ijs:
                dataRGB[i,j,0] = 1.0
                dataRGB[i,j,1] *= 2.0/3
                dataRGB[i,j,2] *= 2.0/3
        
        #self.axes.imshow(data.T, origin='lower', cmap='Greys_r')
        self.axes.imshow(dataRGB, origin='lower')
        self.axes.axis('off')
        self.draw()

class MySpinBox(Spinbox):
    def __init__(self, frame, k):
        Spinbox.__init__(self, frame, from_=0, to=100000, increment=50, width=8)
        self.k = k

class CreateMask(Frame):
    def __init__(self, parent, filename=None, dirname=None, mask4d=False, obj_return_value=None):
        Frame.__init__(self, parent)
        self.parent = parent
        self.obj_return_value = obj_return_value

        if dirname is None:
            self.dirname = ''
        else:
            self.dirname = dirname

        if filename is None or not os.path.isfile(filename):
            if dirname is not None:
                filename = tkFileDialog.askopenfilename(initialdir=dirname)
            else:
                filename = tkFileDialog.askopenfilename()

        if not os.path.isfile(filename):
            parent.destroy()
            return

        self.readImage(filename)
        _, _, filename_mask = create_mask.set_filenames(filename)
        self.filename_mask = filename_mask + '.nii.gz'

        if mask4d:
            if len(self.shape) > 3:
                self.mask_base = np.zeros(self.shape, dtype=np.int8)
                self.mask = self.mask_base[:,:,:,0]
                self.mask4d = True
        else:
            self.mask = np.zeros(self.shape[:3], dtype=np.int8)
            self.mask4d = False
        self.volume = 0

        self.initUI()

        self.modifyUI()
        self.drawSlice()

    def initUI(self):
        self.frame_bottom = Frame(self)
        self.frame_main = Frame(self)

    def modifyUI(self):
        frame_bottom = self.frame_bottom
        frame_main = self.frame_main
        self.sliceView = [ ImshowMplCanvas(frame_main) for k in range(self.shape[2]) ]
        self.sliceSpin = [ MySpinBox(frame_main, k) for k in range(self.shape[2]) ]

        #Layout
        self.parent.title('CreateMask')

        max_spin = 0
        if self.mask4d:
            max_spin = self.shape[3] - 1
        self.spin_volume = Spinbox(frame_bottom, from_=0, to=max_spin, increment=1, command=self.change_volume, width=8)
        self.reset_box(self.spin_volume, 0)

        buttonMask = Button(frame_bottom, text='B0 from Bval...', command=self.createBaseImage)
        buttonMaskCreate = Button(frame_bottom, text='Create', command=self.createMaskAll)
        buttonSave  = Button(frame_bottom, text='Save', command=self.saveMask)

        self.drawmask = BooleanVar()
        buttonDrawMask = Checkbutton(frame_bottom, text='Draw Mask', variable=self.drawmask, command=self.drawSlice)

        num_col = 3
        num_row = 6
        for col in range(self.shape[2]/num_row):
            for k in range(num_row):
                ind = col*num_row + k
                if ind >= self.shape[2]:
                    break
                self.sliceView[ind].get_tk_widget().grid(row=k, column=col*num_col, sticky=NSEW)
                self.sliceSpin[ind].grid(row=k, column=col*num_col+1)
            self.frame_main.grid_columnconfigure(col*num_col, weight=1)

        for k in range(num_row):
            self.frame_main.grid_rowconfigure(k, weight=1)

        self.spin_volume.grid(row=0, column=0)
        buttonMask.grid(row=0, column=1)
        buttonMaskCreate.grid(row=0, column=2)
        buttonSave.grid(row=0, column=3)
        buttonDrawMask.grid(row=0, column=5)

        frame_bottom.pack(side=BOTTOM)
        frame_main.pack(fill=BOTH, expand=TRUE)
        self.pack(fill=BOTH, expand=True)

    def change_volume(self):
        if not self.mask4d:
            return

        self.volume = int(self.spin_volume.get())
        self.data = self.data_base[:,:,:,self.volume]
        self.mask = self.mask_base[:,:,:,self.volume]
        self.drawSlice()


    def saveMask(self):
        print 'save mask to %s' % self.filename_mask
        if self.mask4d:
            img_out = nib.Nifti1Image(self.mask_base, self.img.get_affine(), self.img.get_header())
        else:
            img_out = nib.Nifti1Image(self.mask, self.img.get_affine(), self.img.get_header())
        nib.save(img_out, os.path.join(self.dirname, self.filename_mask))

        if self.obj_return_value is not None:
            self.obj_return_value.delete(0, len(self.obj_return_value.get()))
            self.obj_return_value.insert(0, os.path.basename(self.filename_mask))


    def readImage(self, filename):
        self.filename = filename
        self.img = nib.load(filename)
        self.shape = self.img.shape
        self.data_base = self.img.get_data()
        if len(self.data_base.shape) > 3:
            for f in range(self.data_base.shape[3]):
                for z in range(self.data_base.shape[2]):
                    self.data_base[:,:,z,f] = ndimage.gaussian_filter(self.data_base[:,:,z,f], sigma=0.5)
            self.data = self.data_base[:,:,:,0]
        else:
            for z in range(self.data_base.shape[2]):
                self.data_base[:,:,z] = ndimage.gaussian_filter(self.data_base[:,:,z], sigma=0.5)
            self.data = self.data_base

    def createBaseImage(self):
        bvalFilename = tkFileDialog.askopenfilename()
        if bvalFilename:
            b0frames = parameter.get_b0_from_bval(bvalFilename)
            print self.filename, b0frames
            filename, filename_fltd, filename_mask = create_mask.set_filenames(self.filename)
            create_mask.make_base_image(filename, filename_fltd+'.nii.gz', b0frames)

            self.filename_mask = filename_mask + '.nii.gz'
            self.readImage(filename_fltd + '.nii.gz')
            self.drawSlice()

    def setSpinBoxConnect(self):
        for k in range(self.shape[2]):
            #self.sliceSpin[k].valueChanged.connect(lambda: self.createMaskSlice(k))
            #self.sliceSpin[k].configure(command=lambda:self.createMaskSlice(k, self.sliceSpin[k]))
            self.sliceSpin[k].configure(command=functools.partial(self.createMaskSlice, k, self.sliceSpin[k]) )


    def drawSlice(self):
        for k in range(self.shape[2]):
            self.sliceView[k].drawData(self.data[:,:,k], self.mask[:,:,k], self.drawmask.get())

    def reset_box(self, box, text):
        box.delete(0, len(box.get()))
        box.insert(0, text)

    def createMaskSlice(self, k, value=None):
        self.mask[:,:,k] = 0
        if value is None:
            value_i = int(self.sliceSpin[k].get())
        else:
            value_i = int(value.get())
        #create_mask.mask_from_threshold(self.com_x[k], self.com_y[k], self.sliceSpin[k].value(), self.data, self.mask, k)
        create_mask.mask_from_threshold(self.com_x[k], self.com_y[k], value_i, self.data, self.mask, k)
        self.sliceView[k].drawData(self.data[:,:,k], self.mask[:,:,k], self.drawmask.get())

    def createMaskAll(self):
        self.drawmask.set(True)
        hdr = self.img.get_header()
        zooms = hdr.get_zooms()
        shape = self.shape

        self.com_x, self.com_y, roi_x, roi_y = create_mask.calculate_com(self.data, zooms)

        for k in range(shape[2]):
            dat_subflt = self.data[roi_x[k][0]:roi_x[k][1], roi_y[k][0]:roi_y[k][1], k].flatten()

            #threshold = np.mean(dat_subflt)
            threshold = np.mean(dat_subflt) + np.std(dat_subflt)
            self.reset_box(self.sliceSpin[k], int(threshold))
            self.createMaskSlice(k)
            #create_mask.mask_from_threshold(self.com_x[k], self.com_y[k], threshold, self.data, self.mask, k)

        #self.drawSlice()
        self.setSpinBoxConnect()


if __name__ == '__main__':
    import sys

    root = Tk()

    app = CreateMask(root)
    root.mainloop()


