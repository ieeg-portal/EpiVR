# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: tcarnold@seas.upenn.edu

import FileConvert as fc
import scipy.io
import numpy as np
import pytest


def test_npz2mat():

    fc.npz2mat('test_npz2mat.npz')
    data = scipy.io.loadmat('test_npz2mat.mat')
    assert data['A'][1][1] == 0.34011903728438


def test_mat2npz():

    fc.mat2npz('test_mat2npz.mat')
    data = np.load('test_mat2npz.npz')
    assert data['A'][1][1] == 0.34011903728438

def test_mat2npz_v73():

    fc.mat2npz('test_mat2npz_v73.mat')
    data = np.load('test_mat2npz_v73.npz')
    assert data['A'][1][1] == 0.34011903728438
