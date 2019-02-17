# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: tcarnold@seas.upenn.edu

from RadiologicImage import RadiologicImage as ri
from ComputedImage import ComputedImage as ci
import pytest


# test creating class with valid input
def test_ci_good_input():
    ri_ex = ri('tests/data/mri.nii.gz', "MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])

    # Test RadiologicImage passed as reference
    ci_ex = ci('tests/data/seg.nii.gz', ri_ex, 'SegmentationImage',
               command_used_to_generate='run.sh',
               computation_properties={'notes': 'just a test case'})
    assert ci_ex.computed_type == 'SegmentationImage'
    assert ci_ex.reference_radiologic_image.radiologic_type == "MRI"
    assert ci_ex.reference_radiologic_image.contrast == "T1"
    assert ci_ex.reference_radiologic_image.acquisition_date == "20000123"
    assert ci_ex.reference_radiologic_image.image_dimensions == ["x", "y", "z"]
    assert ci_ex.reference_radiologic_image.image_resolution == [2, 2, 2]
    assert ci_ex.reference_radiologic_image.image_units == ['mm', 'mm', 'mm']

    # Test without command_used_to_generate or computation_properties parameter
    ci_ex = ci('tests/data/seg.nii.gz', ri_ex, 'SegmentationImage')


# test that incorrect user type raises ValueError
def test_ci_bad_type():
    ri_good = ri('tests/data/mri.nii.gz', "MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
    with pytest.raises(IOError):  # Test image_location format
        ci_ex = ci('tests/data/blah.nii.gz', ri_good, 'SegmentationImage')
    with pytest.raises(TypeError):  # Reference radiologic image parameter with incorrect RadiologicImage input format
        ci_ex = ci('tests/data/seg.nii.gz', 3.0, 'SegmentationImage')
    with pytest.raises(TypeError):  # Invalid computational type
        ci_ex = ci('tests/data/seg.nii.gz', ri_good, 'Segmntation')
    with pytest.raises(TypeError):  # Invalid computation properties type
        ci_ex = ci('tests/data/seg.nii.gz', ri_good, 'SegmentationImage', computation_properties=3.0)
