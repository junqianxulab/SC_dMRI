#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')

import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigureCanvas
from matplotlib.figure import Figure

from Tkinter import Tk, Entry, TOP, BOTH, RAISED, N, E, W, S, EW, NS, NW, NSEW
from Tkinter import RIDGE, TRUE, FALSE, VERTICAL, LEFT, RIGHT, X, Y, BOTTOM
from Tkinter import Grid, BooleanVar, StringVar, Canvas, Scrollbar, Spinbox, IntVar
from Tkinter import Checkbutton
#from ttk import Frame, Button, Style, Label, Entry, Checkbutton, Radiobutton
from ttk import Frame, Button, Style, Label, Entry, Radiobutton
import tkFileDialog

import create_mask
import nibabel as nib
import functools
import parameter
import os
from dwi_utils import filename_wo_ext
import struct

class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=4, height=2, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)
        FigureCanvas.__init__(self, fig, master=parent)
        self.show()

class OutlierReg(Frame):
    def __init__(self, parent, filename=None, filename_csv=None, dirname=None, obj_return_value=None):
        Frame.__init__(self, parent)
        self.parent = parent
        self.obj_return_value = obj_return_value

        # check filenames
        if filename is None or not os.path.isfile(filename):
            if dirname is not None:
                filename = tkFileDialog.askopenfilename(initialdir=dirname)
            else:
                filename = tkFileDialog.askopenfilename()

        if not os.path.isfile(filename):
            parent.destroy()
            return

        if dirname is None:
            self.dirname = os.path.dirname(filename)
        else:
            self.dirname = dirname

        if filename_csv is None:
            filename_csv = os.path.join(self.dirname, filename_wo_ext(os.path.basename(filename)) + '.csv')

        self.filename = filename
        self.filename_csv = filename_csv

        n_slices = 0
        n_frames = 0

        with open(filename) as f:
            line = f.readline()
            while True:
                line = f.readline()
                if not line:
                    break
                slice, frame = [int(tmp) for tmp in line.strip().split(',')[:2]]
                if n_slices < slice:
                    n_slices = slice
                if n_frames < frame:
                    n_frames = frame

        n_slices += 1
        n_frames += 1

        self.shape = (n_slices, n_frames)

        #self.distance_x = np.zeros( self.shape, dtype=float )
        #self.distance_y = np.zeros( self.shape, dtype=float )
        #self.distance = np.zeros( self.shape, dtype=float )
        self.values = np.zeros( (n_slices, n_frames, 3), dtype=float )
        self.badenc = np.zeros( self.shape, dtype=np.int16 )
        self.prev = np.zeros( self.shape, dtype=np.int16 )
        self.z = 0

        self.initUI()
        self.run()

        #if filename_csv is not None and os.path.isfile(filename_csv):
        #    self.run_read()

    def reset(self):
        self.badenc[:,:] = np.zeros( self.shape, dtype=np.int16 )

    def make_checkbox_all(self, frame, width=4, ncol=20):
        self.lst_checkbox_slices_values = []
        self.lst_checkbox_slices = []
        ii = 0

        for z in range(self.shape[0]):
            self.lst_checkbox_slices_values.append([BooleanVar() for f in range(self.shape[1])])
            boxes = [Checkbutton(frame, text=str('%s' % f), variable=self.lst_checkbox_slices_values[z][f], width=width) for f in range(self.shape[1])]
            jj = 0
            for f in range(self.shape[1]):
                btn = boxes[f]
                btn.grid(row=z, column=f)
                jj += 1
                if ncol is not None and ncol <= jj:
                    ii += 1
                    jj = 0

            self.lst_checkbox_slices.append(boxes)
            if jj > 0:
                ii += 1

    def make_checkbox(self, frame, width=4, ncol=20):
        #self.lst_checkbox_slices_values = [BooleanVar() for f in range(self.shape[1])]
        self.lst_checkbox_slices_values = [IntVar() for f in range(self.shape[1])]
        self.lst_checkbox_slices = [Checkbutton(frame, text=str('%s' % f), variable=self.lst_checkbox_slices_values[f], width=width, command=functools.partial(self.click_check, f)) for f in range(self.shape[1])]

        ii = 0
        jj = 0
        for f in range(self.shape[1]):
            btn = self.lst_checkbox_slices[f]
            btn.grid(row=ii, column=jj)
            jj += 1
            if ncol is not None and ncol <= jj:
                jj = 0
                ii += 1

        if jj > 0:
            ii += 1


    def click_check(self, f):
        value =  self.lst_checkbox_slices_values[f].get()
        if value == 1:
            value = 0
        else:
            value = 1
        self.lst_checkbox_slices_values[f].set(value)
        self.badenc[self.z, f] = value
        print self.z, f, self.badenc[self.z, f]

    def set_z(self, z=None):
        if z is None:
            z = self.z
        else:
            self.z = z

        for i in range(len(self.lst_checkbox_slices_values)):
            chkbox = self.lst_checkbox_slices_values[i]
            # FIXME
            # avail chkbox
            if self.prev[z,i] > 0 or  self.badenc[z,i] > 0:
                to_check = True
                to_check = 1
                self.lst_checkbox_slices[i].select()
            else:
                to_check = False
                to_check = 0
                self.lst_checkbox_slices[i].deselect()

            chkbox.set(to_check)
            # avail

    def run(self):
        dz = self.shape[0]
        df = self.shape[1]

        #self.values = np.zeros( (dz, df), dtype=float )
        with open(self.filename) as f:
            line = f.readline()
            while True:
                line = f.readline()
                if not line:
                    break
                slice, frame, d_x, d_y, d = line.strip().split(',')[:5]
                if d_x == 'NA':
                    d_x, d_y, d = np.nan, np.nan, np.nan

                self.values[int(slice), int(frame),:] = float(d_x), float(d_y), float(d)

        self.draw_slice()

    def draw_slice(self):
        z = self.z
        df = self.shape[1]

        print z
        #r = self.frame_graph.axes.boxplot(self.values.T)

        sorted_values = np.sort(self.values[z,:])
        q1 = sorted_values[df/4]
        q2 = sorted_values[df/2]
        q3 = sorted_values[df-df/4]
        iqr = q3 - q1
        #if1 = q1 - 1.5*iqr
        #if2 = q3 + 1.5*iqr
        of1 = q1 - 3*iqr
        of2 = q3 + 3*iqr

        xx = np.arange(df)
        z_mean = self.values[z,:].mean()
        z_std  = self.values[z,:].std()
        #z_min  = self.values[z,:].min()
        #z_max  = self.values[z,:].max()

        self.frame_graph.axes.hold(False)
        self.frame_graph.axes.plot([0, df], [z_mean, z_mean], 'k-')
        self.frame_graph.axes.hold(True)
        self.frame_graph.axes.plot([0, df], [z_mean + 1*z_std, z_mean + 1*z_std], 'y--')
        self.frame_graph.axes.plot([0, df], [z_mean + 2*z_std, z_mean + 2*z_std], 'g--')
        self.frame_graph.axes.plot([0, df], [z_mean + 3*z_std, z_mean + 3*z_std], 'b--')
        self.frame_graph.axes.plot([0, df], [of2, of2], 'r-')

        for f in range(df):
            if self.prev[z,f] > 0:
                self.frame_graph.axes.plot(f, self.values[z,f], 'ko')
                self.frame_graph.axes.text(f+0.2, self.values[z,f], str(f))
            elif of1 < self.values[z,f] < of2:
            #if self.values[z,f] < z_mean + 3*z_std:
                self.frame_graph.axes.plot(f, self.values[z,f], 'bo')
                if self.values[z,f] > z_mean + 3*z_std:
                    self.frame_graph.axes.text(f+0.2, self.values[z,f], str(f))
            else:
                self.frame_graph.axes.plot(f, self.values[z,f], 'ro')
                self.frame_graph.axes.text(f+0.2, self.values[z,f], str(f))

        self.frame_graph.draw()
        self.set_z(z)

    def run_read(self):
        filename = tkFileDialog.askopenfilename(initialdir=self.dirname)
        if filename == '':
            return

        with open(filename) as f:
            values = [ [ int(tmp) for tmp in line.strip().split(',') ] for line in f.readlines() ]

        self.prev[:,:] = values

    def run_save(self):
        filename = tkFileDialog.asksaveasfilename(initialdir=os.path.dirname(self.filename_csv), initialfile=os.path.basename(self.filename_csv))
        if filename == '':
            return

        dz = self.shape[0]
        df = self.shape[1]

        badenc = self.badenc.copy()
        badenc[self.prev > 0] = 1

        with open(filename, 'w') as f:
            for z in range(dz):
                f.write('%s\n' % (','.join([str(tmp) for tmp in badenc[z,:]])))

        if False:
            with open(os.path.join(self.dirname, 'badenc.dat'), 'wb') as fout:
                if False:
                    for z in range(dz):
                        row = struct.pack('i'*df, *badenc[z,:])
                        fout.write(row)

                fout.write('%s %s\n' % (dz, df))
                for z in range(dz):
                    row = ' '.join(str(tmp) for tmp in badenc[z,:])
                    fout.write('%s\n' % row)

        if self.obj_return_value is not None:
            self.obj_return_value.delete(0, len(self.obj_return_value.get()))
            self.obj_return_value.insert(0, filename)

    def initUI(self):
        self.frame_top = Frame(self)
        self.frame_graph = MplCanvas(self)
        self.frame_bottom = Frame(self)

        Label(self.frame_top, text='Z = ').pack(side=LEFT)
        self.spin_z = Spinbox(self.frame_top, from_=0, to=self.shape[0]-1, increment=1, command=self.change_z)
        self.spin_z.pack()
        self.make_checkbox(self.frame_bottom)

        button_reset = Button(self.frame_top, text='Reset', command=self.reset).pack()

        Label(self.frame_top, text='CSV').pack()
        self.txt_filename_csv = Entry(self.frame_top)
        self.txt_filename_csv.pack()
        self.button_read = Button(self.frame_top, text='Read', command=self.run_read)
        self.button_read.pack()
        self.button_save = Button(self.frame_top, text='Save', command=self.run_save)
        self.button_save.pack()

        self.frame_top.pack(side=TOP)
        self.frame_graph.get_tk_widget().pack(fill=BOTH, expand=TRUE)
        self.frame_bottom.pack(fill=BOTH, expand=TRUE)
        self.pack(fill=BOTH, expand=True)

    def change_z(self):
        self.z = int(self.spin_z.get())
        self.draw_slice()

    def reset_box(self, box, text):
        box.delete(0, len(box.get()))
        box.insert(0, text)

if __name__ == '__main__':
    import sys

    root = Tk()

    app = OutlierReg(root)
    root.mainloop()

