# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: lkini@pennmedicine.upenn.edu

from RadiologicImage import RadiologicImage as ri
from RadiologicImageSet import RadiologicImageSet as ris
import pytest


# test creating class with valid input
def test_ri_good_input():
    image1 = ri("MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
    image2 = ri("MRI", "T2", "20000123", ["x", "y", "z"], [0.5, 0.5, 4])
    image3 = ri("MRI", "T1+GAD", "20000123", ["x", "y", "z"], [1, 1, 1.5])
    image_set = ris([image1, image2, image3])

    # Manually create filtered set after filtering by contrast for T2 images only
    filtered_set2 = ris([image2])

    assert image_set.radiology_images == [image1, image2, image3]
    assert image_set.filter_by_radiology_type('MRI').radiology_images == image_set.radiology_images
    assert image_set.filter_by_contrast('T2').radiology_images == filtered_set2.radiology_images
    assert image_set.filter_by_acquisition_date('20000123').radiology_images == image_set.radiology_images


# test filter_by_radiology_type raises error on inappropriate radiology type
def test_ris_radiology_type_filter():
    with pytest.raises(TypeError, match=r'*not a proper image type*'):
        image1 = ri("MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
        image2 = ri("MRI", "T2", "20000123", ["x", "y", "z"], [0.5, 0.5, 4])
        image3 = ri("MRI", "T1+GAD", "20000123", ["x", "y", "z"], [1, 1, 1.5])
        image_set = ris([image1, image2, image3])
        image_set.filter_by_radiology_type('SPECT')


# test filter_by_contrast raises error on inappropriate contrast
def test_ris_contrast_filter():
    with pytest.raises(TypeError, match=r'*contrast is not a string*'):
        image1 = ri("MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
        image2 = ri("MRI", "T2", "20000123", ["x", "y", "z"], [0.5, 0.5, 4])
        image3 = ri("MRI", "T1+GAD", "20000123", ["x", "y", "z"], [1, 1, 1.5])
        image_set = ris([image1, image2, image3])
        image_set.filter_by_contrast(1337)  # 1337 is not a proper contrast


# test filter_by_acquisition_date raises error on non-string acquisition date
def test_ris_acquisition_date_str_filter():
    with pytest.raises(TypeError, match=r'*date is not a string*'):
        image1 = ri("MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
        image2 = ri("MRI", "T2", "20000123", ["x", "y", "z"], [0.5, 0.5, 4])
        image3 = ri("MRI", "T1+GAD", "20000123", ["x", "y", "z"], [1, 1, 1.5])
        image_set = ris([image1, image2, image3])
        image_set.filter_by_acquisition_date(1337)  # 1337 is not a proper string format for date


# test filter_by_acquisition_date raises error on incorrectly formatted acquisition date
def test_ris_acquisition_date_format_filter():
    with pytest.raises(TypeError, match=r'*should follow YYYYMMDD date format*'):
        image1 = ri("MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
        image2 = ri("MRI", "T2", "20000123", ["x", "y", "z"], [0.5, 0.5, 4])
        image3 = ri("MRI", "T1+GAD", "20000123", ["x", "y", "z"], [1, 1, 1.5])
        image_set = ris([image1, image2, image3])
        image_set.filter_by_acquisition_date('01/23/00')  # 01/23/00 is not a proper string format for date
    with pytest.raises(TypeError, match=r'*should follow YYYYMMDD date format*'):
        image1 = ri("MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
        image2 = ri("MRI", "T2", "20000123", ["x", "y", "z"], [0.5, 0.5, 4])
        image3 = ri("MRI", "T1+GAD", "20000123", ["x", "y", "z"], [1, 1, 1.5])
        image_set = ris([image1, image2, image3])
        image_set.filter_by_acquisition_date('012300')  # 012300 is not a proper string format for date
