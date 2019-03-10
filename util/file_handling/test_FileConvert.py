# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: tcarnold@seas.upenn.edu
import os
import scipy.io
import numpy as np
import FileConvert as fc


def test_npz2mat():
    fc.npz2mat('tests/data/test_npz2mat.npz')
    data = scipy.io.loadmat('tests/data/test_npz2mat.mat')
    assert data['A'][1][1] == 0.34011903728438
    os.system('rm tests/data/test_npz2mat.mat')


def test_mat2npz():
    fc.mat2npz('tests/data/test_mat2npz.mat')
    data = np.load('tests/data/test_mat2npz.npz')
    assert data['A'][1][1] == 0.34011903728438
    os.system('rm tests/data/test_mat2npz.npz')


def test_mat2npz_v73():
    fc.mat2npz_v73('tests/data/test_mat2npz_v73.mat')
    data = np.load('tests/data/test_mat2npz_v73.npz')
    assert data['A'][1][1] == 0.34011903728438
    os.system('rm tests/data/test_mat2npz_v73.npz')
