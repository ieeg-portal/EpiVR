# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: tcarnold@seas.upenn.edu

"""fileConvert

FileConvert is a series of simple functions for transferring data between
programming languages. Currently it can be used to generate:
.mat -> .npz
.npz -> .mat

Example:
    import fileConvert as fc
    fc.npz2mat('filename.npz')

Todo:
    * add csv functionality (read in and write out)

Alterations:
    TCA 2/21/19 - initialized fileConvert (npz2mat)

"""
import numpy as np
import scipy.io
import h5py


def npz2mat(filename):

    data = np.load(filename)  # load in
    fileout = filename[:-3] + 'mat'  # append with proper suffix
    scipy.io.savemat(fileout, data)


def mat2npz(filename):

    data = scipy.io.loadmat(filename)  # load in
    data = np.array(data)  # For converting to numpy array
    data.pop('__header__')  # remove irrelevant variables from dict (loadmat)
    data.pop('__version__')
    data.pop('__globals__')
    fileout = filename[:-4]  # get filename without suffix
    np.savez(fileout, **data)

def mat2npz_v73(filename):

    # Try h5py first
    f = h5py.File(filename, 'r')
    data = {}
    for k, v in f.items():
        print(k, ' ', v)
        data[k] = np.array(v)

    fileout = filename[:-4]  # get filename without suffix
    np.savez(fileout, **data)
