# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: tcarnold@seas.upenn.edu

from RadiologicImage import RadiologicImage as ri
import pytest


# test creating class with valid input
def test_ri_good_input():
    ex = ri("MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
    assert ex.radiologic_type == "MRI"
    assert ex.contrast == "T1"
    assert ex.acquisition_date == "20000123"
    assert ex.image_dimensions == ["x", "y", "z"]
    assert ex.image_resolution == [2, 2, 2]
    assert ex.image_units == ['mm', 'mm', 'mm']


# test that incorrect user type raises ValueError
def test_ri_bad_type():
    with pytest.raises(TypeError, match=r'.* is not a valid image type..*'):
        ex = ri("Garbage", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
