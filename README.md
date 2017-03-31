# SC_dMRI
Spinal cord diffusion MRI preprocessing tool

## Dependency
#### python 2.7:
Only tested on 2.7.
I recommend [pip](https://pip.pypa.io/en/stable/installing/) for the following module installation. 
- numpy
- scipy
- nibabel
- matplotlib

#### [python-tk](https://wiki.python.org/moin/TkInter):
- Ubuntu/Debian: `sudo apt-get install python-tk`
- Fedora: `yum install tkinter`

#### [FSL 5.0.x](https://fsl.fmrib.ox.ac.uk/fsldownloads/fsldownloadmain.html)

#### 4dfp
- for DTI map generation with outlier rejection

## Installation
We don't have any fancy installation method yet.
1. Download to any directory: `git clone git@github.com:junqianxulab/SC_dMRI.git`
2. Give executable permission to the main file: `chmod +x {download_directory}/SC_dMRI/csp_preprocessing/spine_reg_gui.py`
3. Make a symbolic link in a user directory (e.g. `ln -s {download_directory}/SC_dMRI/csp_preprocessing/spine_reg_gui.py ${HOME}/bin`) or a system directory (e.g. `ln -s {download_directory}/SC_dMRI/csp_preprocessing/spine_reg_gui.py /local/bin`).
- Note that the directory containing the symbolic link should be in a PATH environment.
4. Try `spine_reg_gui.py`

## Documentation
- [Manual](manual/manual.md)

#### under construction
