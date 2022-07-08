import matplotlib.pyplot as plt
import numpy as np
import Arriero as ar
import Eduviges.extractor as extr
import matplotlib.cm as cm
import matplotlib

outfile = './steady_state_test.h5'

plots_folder = '/home/zach/Documents/Code_Projects/paramo/Plots/'


def run_bw1D_afterglow():
    rr = ar.Runner(comp_kw={'OMP': True, 'HDF5': True})
    rr.run_Aglow(clean=False)


run_bw1D_afterglow()