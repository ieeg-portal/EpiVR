# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: tcarnold@seas.upenn.edu

"""RadiologicImage.py unit test
python -m pytest -v test_RadiologicImage.py
"""

from RadiologicImage import RadiologicImage as ri
import pytest


# test creating class with valid input
def test_ri_good_input():
    ex = ri("MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
    assert ex.radiologic_type == "MRI"
    assert ex.contrast == "T1"
    assert ex.date_of_acquisition == "20000123"
    assert ex.image_dimensions == ["x", "y", "z"]
    assert ex.image_resolution == [2, 2, 2]


# test that incorrect user type raises ValueError
def test_ri_bad_type():
    with pytest.raises(ValueError, match=r'.* is not a valid image type..*'):
        ex = ri("Garbage", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
