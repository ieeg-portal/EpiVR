# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: tcarnold@seas.upenn.edu

"""fileConvert

fileconvert is a series of simple function for transfer of data between
programming languages. Currently it can be used to generate .mat files or to
load in .mat files into python

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


def npz2mat(filename):

    data = np.load(filename)  # load in
    fileout = filename[:-3] + 'mat'  # append with proper suffix
    scipy.io.savemat(fileout, data)


def mat2npz(filename):

    data = scipy.io.loadmat(filename)
    fileout = filename[:-3] + 'npz'  # append with proper suffix
    np.save(fileout, data)
