# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: lkini@pennmedicine.upenn.edu

from RadiologicImage import RadiologicImage as ri
from ComputedImage import ComputedImage as ci
from ComputedImageSet import ComputedImageSet as cis
import pytest


# test creating class with valid input
def test_cis_good_input():
    ri_ex = ri('tests/data/mri.nii.gz', 'MRI', 'T1', '20000123', ['x', 'y', 'z'], [2, 2, 2])

    # Test RadiologicImage passed as reference
    ci_ex1 = ci('tests/data/seg1.nii.gz', ri_ex, 'SegmentationImage', 'NA', 'NA',
                command_used_to_generate='run1.sh')
    ci_ex2 = ci('tests/data/seg2.nii.gz', ri_ex, 'RegistrationImage', 'MRI', 'T1',
                command_used_to_generate='run2.sh')
    image_set = cis([ci_ex1, ci_ex2])

    # Manually create filtered set after filtering by contrast for T2 images only
    filtered_set = cis(ci_ex2)
    assert image_set.computed_images == [ci_ex1, ci_ex2]
    assert image_set.filter_by_radiology_type('MRI').computed_images == filtered_set.computed_images
    assert image_set.filter_by_computed_type('RegistrationImage').computed_images == filtered_set.computed_images
    assert image_set.filter_by_contrast('T1').computed_images == filtered_set.computed_images


# test creating class with invalid inputs
def test_cis_bad_input():
    with pytest.raises(ValueError):
        image_set = cis(3.0)


# test filter_by_radiology_type raises error on inappropriate radiology type
def test_cis_radiology_type_filter():
    with pytest.raises(TypeError):
        ri_ex = ri('tests/data/mri.nii.gz', 'MRI', 'T1', '20000123', ['x', 'y', 'z'], [2, 2, 2])
        ci_ex1 = ci('tests/data/seg1.nii.gz', ri_ex, 'SegmentationImage', 'NA', 'NA',
                    command_used_to_generate='run1.sh')
        ci_ex2 = ci('tests/data/seg2.nii.gz', ri_ex, 'RegistrationImage', 'MRI', 'T1',
                    command_used_to_generate='run2.sh')
        image_set = cis([ci_ex1, ci_ex2])
        image_set.filter_by_radiology_type('NIRS')


# test filter_by_contrast raises error on inappropriate contrast
def test_cis_contrast_filter():
    with pytest.raises(TypeError):
        ri_ex = ri('tests/data/mri.nii.gz', 'MRI', 'T1', '20000123', ['x', 'y', 'z'], [2, 2, 2])
        ci_ex1 = ci('tests/data/seg1.nii.gz', ri_ex, 'SegmentationImage', 'NA', 'NA',
                    command_used_to_generate='run1.sh')
        ci_ex2 = ci('tests/data/seg2.nii.gz', ri_ex, 'RegistrationImage', 'MRI', 'T1',
                    command_used_to_generate='run2.sh')
        image_set = cis([ci_ex1, ci_ex2])
        image_set.filter_by_contrast(1337)  # 1337 is not a proper contrast


# test filter_by_computed_type raises error on inappropriate computed type
def test_cis_computed_type_filter():
    with pytest.raises(TypeError):
        ri_ex = ri('tests/data/mri.nii.gz', 'MRI', 'T1', '20000123', ['x', 'y', 'z'], [2, 2, 2])
        ci_ex1 = ci('tests/data/seg1.nii.gz', ri_ex, 'SegmentationImage', 'NA', 'NA',
                    command_used_to_generate='run1.sh')
        ci_ex2 = ci('tests/data/seg2.nii.gz', ri_ex, 'RegistrationImage', 'MRI', 'T1',
                    command_used_to_generate='run2.sh')
        image_set = cis([ci_ex1, ci_ex2])
        image_set.filter_by_computed_type('Segmntation')
