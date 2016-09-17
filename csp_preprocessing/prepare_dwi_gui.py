from Tkinter import Tk, Entry, TOP, BOTH, RAISED, N, E, W, S, EW, NS, NW, NSEW
from Tkinter import RIDGE, TRUE, FALSE, VERTICAL, LEFT, RIGHT, X, Y
from Tkinter import Grid, BooleanVar, StringVar, Canvas, Scrollbar
from ttk import Frame, Button, Style, Label, Entry, Checkbutton, Radiobutton
import tkFileDialog
import os
import shutil
from dwi_utils import filename_wo_ext, extname, create_merge, create_acqparams, create_index, abspath_to_relpath
import parameter

COPY_AS_SYMBOLIC_LINK = True

def filenameDialog_text(text):
    fl = tkFileDialog.askopenfilename()
    if fl != '':
        text.delete(0, len(text.get()))
        text.insert(0, fl)

def dirnameDialog_text(text):
    fl = tkFileDialog.askdirectory(initialdir=text.get())
    if fl != '':
        text.delete(0, len(text.get()))
        text.insert(0, fl)

class VerticalScrolledFrame(Frame):
    """A pure Tkinter scrollable frame that actually works!
    * Use the 'interior' attribute to place widgets inside the scrollable frame
    * Construct and pack/place/grid normally
    * This frame only allows vertical scrolling

    """
    def __init__(self, parent, *args, **kw):
        Frame.__init__(self, parent, *args, **kw)

        # create a canvas object and a vertical scrollbar for scrolling it
        vscrollbar = Scrollbar(self, orient=VERTICAL)
        vscrollbar.pack(fill=Y, side=RIGHT, expand=FALSE)
        canvas = Canvas(self, bd=0, highlightthickness=0,
                        yscrollcommand=vscrollbar.set)
        canvas.pack(side=LEFT, fill=BOTH, expand=TRUE)
        vscrollbar.config(command=canvas.yview)

        # reset the view
        canvas.xview_moveto(0)
        canvas.yview_moveto(0)

        # create a frame inside the canvas which will be scrolled with it
        self.interior = interior = Frame(canvas)
        interior_id = canvas.create_window(0, 0, window=interior,
                                           anchor=NW)

        # track changes to the canvas and frame width and sync them,
        # also updating the scrollbar
        def _configure_interior(event):
            # update the scrollbars to match the size of the inner frame
            size = (interior.winfo_reqwidth(), interior.winfo_reqheight())
            canvas.config(scrollregion="0 0 %s %s" % size)
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # update the canvas's width to fit the inner frame
                canvas.config(width=interior.winfo_reqwidth())
        interior.bind('<Configure>', _configure_interior)

        def _configure_canvas(event):
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # update the inner frame's width to fill the canvas
                canvas.itemconfigure(interior_id, width=canvas.winfo_width())
        canvas.bind('<Configure>', _configure_canvas)

class PrepareDWI(Frame):

    def __init__(self, parent, subject='', output_dir='', return_value=None, obj_return_value=None, make_param=False):
        Frame.__init__(self, parent)
        self.make_param = make_param
        self.return_value = return_value
        self.obj_return_value = obj_return_value
        self.parent = parent
        self.subject_name = subject
        self.output_dir = output_dir
        self.initUI()

        fn_dwi_filenames = os.path.join(self.txt_output.get(), '%sDWI.txt' % self.prefix())
        if os.path.isfile(fn_dwi_filenames):
            self.read_dwi_filenames(fn_dwi_filenames)

    def read_dwi_filenames(self, filename):
        with open(filename) as f:
            for a_name in f.readlines():
                self.lst_dwi.add_one(a_name.strip())

    def save_dwi_filenames(self, filename):
        with open(filename, 'w') as f:
            for a_dwi in self.lst_dwi.lst_dwi:
                f.write('%s\n' % a_dwi.text_from.get())

    class List_DWI:
        def __init__(self, parent, frame):
            self.parent = parent
            self.frame = frame
            self.lst_dwi = []
            self.len = 0

        def size(self):
            return len(self.lst_dwi)

        def set_filenames_to(self, middle=None):
            for a_dwi in self.lst_dwi:
                a_dwi.set_filename_to(middle=middle)

        def get_dw_filenames(self, with_bval_bvec=True):
            rtn = []
            for a_dwi in self.lst_dwi:
                a_name = a_dwi.get_dwi_filenames()
                if with_bval_bvec:
                    rtn += a_name
                else:
                    rtn += a_name[:1]
            return rtn

        def get_b0_filenames(self):
            rtn = []
            for a_dwi in self.lst_dwi:
                rtn += a_dwi.get_b0_filename()
            return rtn

        def prefix(self):
            return self.parent.prefix()

        def fn_b0(self, text=None):
            filename = '%sB0Ref.nii.gz' % (self.prefix())
            if text is not None:
                text.delete(0, len(text.get()))
                text.insert(0, filename)
            return filename

        def fn_dw(self, text=None):
            has_ap = False
            has_pa = False
            n_dwi = 0
            for a_dwi in self.lst_dwi:
                if not a_dwi.is_dw.get():
                    continue
                n_dwi += 1
                #if 'ap' in a_dwi.text_to.get().lower():
                if 'AP' == a_dwi.direction.get():
                    has_ap = True
                #if 'pa' in a_dwi.text_to.get().lower():
                if 'PA' == a_dwi.direction.get():
                    has_pa = True
            if has_ap and not has_pa:
                middle = 'AP_'
            elif has_pa and not has_ap:
                middle = 'PA_'
            else:
                middle = 'DWI_'

            filename = '%s%s%smerged.nii.gz' % (self.prefix(), middle, n_dwi)
            if text is not None:
                text.delete(0, len(text.get()))
                text.insert(0, filename)

            return filename

        def get_b0_directions(self):
            return  [ a_dwi.direction.get() for a_dwi in self.lst_dwi if a_dwi.is_b0.get()]
        def get_dw_directions(self):
            return  [ a_dwi.direction.get() for a_dwi in self.lst_dwi if a_dwi.is_dw.get()]

        def add_one(self, text=''):
            a_dwi = self.A_DWI(self, self.frame, text=text, row=self.len, column=0)
            fn_bval = filename_wo_ext(text) + '.bval'
            if os.path.isfile(fn_bval):
                a_dwi.is_b0.set(False)
                a_dwi.is_dw.set(True)
            else:
                a_dwi.is_b0.set(True)
                a_dwi.is_dw.set(False)
            a_dwi.set_appa()
            a_dwi.set_filename_to()
            self.len += 1
            self.lst_dwi.append(a_dwi)

        def remove(self, row):
            a_dwi = self.lst_dwi.pop(row)
            self.len -= 1
            for i in range(row, self.len):
                self.lst_dwi[i].change_row(i)
            del a_dwi

        def up(self, row):
            if row == 0:
                return
            self.lst_dwi[row].change_row(row-1)
            self.lst_dwi[row-1].change_row(row)
            self.lst_dwi[row], self.lst_dwi[row-1] = self.lst_dwi[row-1], self.lst_dwi[row]

        def dn(self, row):
            if row >= self.len-1:
                return
            self.lst_dwi[row].change_row(row+1)
            self.lst_dwi[row+1].change_row(row)
            self.lst_dwi[row], self.lst_dwi[row+1] = self.lst_dwi[row+1], self.lst_dwi[row]

        class A_DWI:
            def __init__(self, container, frame, label='', text='', row=0, column=0):
                self.container = container
                self.is_b0 = BooleanVar(container.parent)
                self.is_dw = BooleanVar(container.parent)
                self.column = column
                self.direction = StringVar(container.parent)

                self.label_from = Label(frame, text='from')
                self.text_from = Entry(frame)
                self.text_from.insert(0, text)
                self.button_file_from = Button(frame, text='...', command=lambda:filenameDialog_text(self.text_from))
                self.button_rm = Button(frame, text='remove', command=self.click_remove)
                self.radio_ap = Radiobutton(frame, text='AP', variable=self.direction, value='AP', command=self.set_direction)
                self.radio_pa = Radiobutton(frame, text='PA', variable=self.direction, value='PA', command=self.set_direction)

                self.label_to = Label(frame, text='to')
                self.text_to = Entry(frame)
                #self.text_to.insert(0, text)
                self.button_file_to = Button(frame, text='Gen', command=self.set_filename_to)
                self.check_b0 = Checkbutton(frame, text='B0',   variable=self.is_b0)
                self.check_dw = Checkbutton(frame, text='DWI',  variable=self.is_dw)

                self.button_up = Button(frame, text='up', width=3, command=self.click_up)
                self.button_dn = Button(frame, text='down', width=3, command=self.click_dn)

                self.row = -1
                self.change_row(row)
                if text != '':
                    self.set_appa()
                    self.set_filename_to()

            def prefix(self):
                return self.container.parent.prefix()

            def set_direction(self):
                pass

            def get_dwi_filenames(self):
                '''
                :return: [('from', 'to'), ('from', 'to')]
                '''
                filename_from = self.text_from.get()
                filename_to = self.text_to.get()
                if self.is_dw.get():
                    rtn = [ [filename_from, filename_to] ]
                    filename_b_from = filename_wo_ext(filename_from)
                    filename_b_to = filename_wo_ext(filename_to)
                    rtn.append( [filename_b_from+'.bval', filename_b_to+'.bval'] )
                    rtn.append( [filename_b_from+'.bvec', filename_b_to+'.bvec'] )
                    return rtn
                return []

            def get_b0_filename(self):
                '''
                :return: [('from', 'to')]
                '''
                filename_from = self.text_from.get()
                filename_to = self.text_to.get()
                ext = extname(filename_to)
                if self.is_b0.get():
                    if self.is_dw.get():
                        filename_to = '%s_B0%s' % (filename_wo_ext(filename_to), ext)
                    return [ [filename_from, filename_to] ]
                return []

            def set_appa(self):
                filename_from = self.text_from.get()
                basename_from = os.path.basename(filename_from)
                basename_b_from = filename_wo_ext(basename_from)
                if 'pa' in basename_b_from.lower():
                    self.direction.set('PA')
                elif 'ap' in basename_b_from.lower():
                    self.direction.set('AP')
                else:
                    pass

            def set_filename_to(self, middle=None):
                filename_from = self.text_from.get()
                basename_from = os.path.basename(filename_from)
                number = os.path.dirname(filename_from).split('/')[-1].split('_')[0]
                if number == '':
                    number = str(self.row)
                else:
                    try:
                        int(number)
                    except:
                        number = str(self.row)


                ext = extname(basename_from)
                intermediate = self.direction.get()
                if intermediate == '':
                    intermediate = 'DWI'

                if self.is_b0.get() and not self.is_dw.get():
                    intermediate += '_B0'

                self.text_to.delete(0, len(self.text_to.get()))
                self.text_to.insert(0, '%s%s_%s%s' % (self.prefix(), number, intermediate, ext))

            def change_row(self, row):
                if self.row == row:
                    return
                self.row = row
                i = 2*row
                j = self.column
                j += 0; self.button_up.grid(row=i, column=j)
                j += 1; self.label_from.grid(row=i, column=j)
                j += 1; self.text_from.grid(row=i, column=j, sticky=EW)
                j += 1; self.button_file_from.grid(row=i, column=j)
                j += 1; self.button_rm.grid(row=i, column=j)
                j += 1; self.radio_ap.grid(row=i, column=j)
                j += 1; self.radio_pa.grid(row=i, column=j)
                i += 1
                j = 0
                j += 0; self.button_dn.grid(row=i, column=j)
                j += 1; self.label_to.grid(row=i, column=j)
                j += 1; self.text_to.grid(row=i, column=j, sticky=EW)
                j += 1; self.button_file_to.grid(row=i, column=j)
                j += 1
                j += 1; self.check_b0.grid(row=i, column=j)
                j += 1; self.check_dw.grid(row=i, column=j)

            def click_remove(self):
                self.container.remove(self.row)

                self.button_up.destroy()
                self.label_from.destroy()
                self.text_from.destroy()
                self.button_file_from.destroy()
                self.button_rm.destroy()
                self.radio_ap.destroy()
                self.radio_pa.destroy()
                self.button_dn.destroy()
                self.label_to.destroy()
                self.text_to.destroy()
                self.button_file_to.destroy()
                self.check_b0.destroy()
                self.check_dw.destroy()

            def click_up(self):
                self.container.up(self.row)
            def click_dn(self):
                self.container.dn(self.row)

    def add_dwi_list(self):
        while True:
            filenames = tkFileDialog.askopenfilenames()
            if len(filenames) < 1:
                break
            for filename in filenames:
                if extname(filename) == '.nii.gz' or extname(filename) == '.nii':
                    self.lst_dwi.add_one(text=filename)

    def subject(self):
        return self.txt_subject.get()

    def initUI(self):
        self.parent.title('Prepare DWI')
        frm_output = Frame(self)
        frm_button = Frame(self)
        frm_list = VerticalScrolledFrame(self)
        self.frm_list_sub = frm_list.interior

        self.lst_dwi = self.List_DWI(self, self.frm_list_sub)

        i = 0
        frm_output.grid(row=i, column=0, columnspan=6, sticky=EW)
        frm_output.grid_columnconfigure(1, weight=1)
        ii  = 0; Label(frm_output, text='Subject').grid(row=ii, column=0)
        self.txt_subject = Entry(frm_output); self.txt_subject.grid(row=ii, column=1, sticky=EW)
        self.txt_subject.insert(0, self.subject_name)
        ii += 1; Label(frm_output, text='Output Directory').grid(row=ii, column=0)
        self.txt_output = Entry(frm_output); self.txt_output.grid(row=ii, column=1, sticky=EW)
        self.txt_output.insert(0, self.output_dir)
        btn_output = Button(frm_output, text='...', command=lambda:dirnameDialog_text(self.txt_output)); btn_output.grid(row=ii, column=2)
        ii += 1; Label(frm_output, text='Merged B0').grid(row=ii, column=0)
        self.txt_output_b0 = Entry(frm_output); self.txt_output_b0.grid(row=ii, column=1, sticky=EW)
        btn_output_b0 = Button(frm_output, text='Gen', command=lambda:self.lst_dwi.fn_b0(self.txt_output_b0)); btn_output_b0.grid(row=ii, column=2)
        ii += 1; Label(frm_output, text='Merged DWI').grid(row=ii, column=0)
        self.txt_output_dw = Entry(frm_output); self.txt_output_dw.grid(row=ii, column=1, sticky=EW)
        btn_output_dw = Button(frm_output, text='Gen', command=lambda:self.lst_dwi.fn_dw(self.txt_output_dw)); btn_output_dw.grid(row=ii, column=2)

        i += ii
        i += 1
        frm_button.grid(row=i, column=0, columnspan=6, sticky=EW)

        ii = 0
        btn_add = Button(frm_button, text='ADD Nifti1 file', command=self.add_dwi_list); btn_add.grid(row=ii, column=0)
        btn_run = Button(frm_button, text='Run', command=self.run); btn_run.grid(row=ii, column=1)

        i += ii
        i += 1
        i_frm_list = i

        self.frm_list_sub.grid_columnconfigure(2, weight=1)
        frm_list.grid(row=i, column=0, columnspan=6, sticky=NSEW)
        self.grid_rowconfigure(i_frm_list, weight=1, minsize=20)

        #self.lst_dwi.add_one('test1')
        #self.lst_dwi.add_one('test2')


        self.grid_columnconfigure(0, weight=1)

        ni = i + 1
        #for i in range(i_frm_list, ni):
        #    self.grid_rowconfigure(i, weight=1, minsize=20)

        self.pack(fill=BOTH, expand=True)

    #def configure_interior(self, event):
    def configure_interior(self, canvas):
        #size = self.frm_list_sub.winfo_reqwidth(), self.frm_list_sub.winfo_reqheight()
        #self.canvas.config(scrollregion="0 0 %s %s" % size)
        #if self.frm_list_sub.winfo_reqwidth() != self.canvas.winfo_width():
        #    self.canvas.config(width=self.frm_list_sub.winfo_width())
        canvas.configure(scrollregion=canvas.bbox("all"))

    def configure_canvas(self, event):
        '''Reset the scroll region to encompass the inner frame'''
        if self.frm_list_sub.winfo_reqwidth() != self.canvas.winfo_width():
            self.canvas.itemconfigure(self.interior_id, width=self.canvas.winfo_width())

    def prefix(self):
        prefix = self.subject()
        if prefix != '':
            prefix = prefix + '_'
        return prefix

    def run(self):
        if self.lst_dwi.size() == 0:
            return

        dir_out = os.path.abspath(self.txt_output.get())
        os.chdir(dir_out)

        # copy Nifti1 files
        filenames_b0 = self.lst_dwi.get_b0_filenames()
        for a_pair in filenames_b0:
            if COPY_AS_SYMBOLIC_LINK:
                rel_path = abspath_to_relpath(dir_out, os.path.abspath(os.path.dirname(a_pair[0])))
                if os.path.lexists(a_pair[1]):
                    os.remove(a_pair[1])
                os.symlink(os.path.join(rel_path, os.path.basename(a_pair[0])), a_pair[1])
                a_pair[1] = os.path.join(dir_out, a_pair[1])
            else:
                a_pair[1] = os.path.join(dir_out, a_pair[1])
                shutil.copy(*a_pair)
        #print filenames_b0

        filenames_dw = self.lst_dwi.get_dw_filenames()
        for a_pair in filenames_dw:
            if COPY_AS_SYMBOLIC_LINK:
                rel_path = abspath_to_relpath(dir_out, os.path.abspath(os.path.dirname(a_pair[0])))
                if os.path.lexists(a_pair[1]):
                    os.remove(a_pair[1])
                os.symlink(os.path.join(rel_path, os.path.basename(a_pair[0])), a_pair[1])
                a_pair[1] = os.path.join(dir_out, a_pair[1])
            else:
                a_pair[1] = os.path.join(dir_out, a_pair[1])
                shutil.copy(*a_pair)
        #print filenames_dw

        fn_b0s = [a_pair[1] for a_pair in filenames_b0]
        fn_dws = [a_pair[1] for a_pair in filenames_dw if extname(a_pair[1])[-5:-2] != '.bv']
        print fn_dws

        dir_b0s = self.lst_dwi.get_b0_directions()
        dir_dws = self.lst_dwi.get_dw_directions()

        # merge Nifti1, bval, bvec files
        fn_out_b0 = self.txt_output_b0.get()
        if fn_out_b0 == '':
            fn_out_b0 = self.lst_dwi.fn_b0(self.txt_output_b0)
        fn_out_dw = self.txt_output_dw.get()
        if fn_out_dw == '':
            fn_out_dw = self.lst_dwi.fn_dw(self.txt_output_dw)

        create_merge(os.path.join(dir_out, fn_out_b0), fn_b0s, verbose=True)
        create_merge(os.path.join(dir_out, fn_out_dw), fn_dws, is_dwi=True, verbose=True)

        # create acqparams.txt, index.txt
        create_acqparams(os.path.join(dir_out,'%sacqparams.txt' % self.prefix()),
                         fn_b0s, direction=dir_b0s)
        create_index(os.path.join(dir_out,'%sindex.txt' % self.prefix()),
                     fn_b0s, fn_dws, dir_b0s, dir_dws)

        fn_dwi_filenames = os.path.join(dir_out, '%sDWI.txt' % self.prefix())
        self.save_dwi_filenames(fn_dwi_filenames)

        print 'Done'

        if self.make_param:
            fn_param = os.path.join(dir_out, '%sparams' % self.prefix())

            param = parameter.Parameter()
            if os.path.isfile(fn_param):
                param.read(fn_param)
            param.fn_b0 = os.path.join(dir_out, fn_out_b0)
            param.fn_dwi = os.path.join(dir_out, fn_out_dw)
            param.subject = self.subject()
            param.working_dir = dir_out
            param.save(fn_param)

            print 'Saved : %s' % fn_param

        if self.return_value is not None:
            self.return_value['subject'] = self.subject()
            self.return_value['b0'] = fn_out_b0
            self.return_value['dwi'] = fn_out_dw
            self.return_value['output_directory'] = dir_out

        if self.obj_return_value is not None:
            self.obj_return_value.update_parameter()
            self.parent.destroy()
            return

        return self.subject(), dir_out, fn_out_b0, fn_out_dw


def main():

    root = Tk()
    #root.geometry("400x600+100+100")
    app = PrepareDWI(root)
    root.mainloop()

if __name__ == '__main__':
    main()
