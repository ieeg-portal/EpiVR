# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: tcarnold@seas.upenn.edu

"""MatchElectrodes

This function takes in three files.
    1. resection/ablation roi image
    2. electrodes rois image
    3. electrode labels file

The output of the function is a list of the electrode labels that have been
resected.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * For module TODOs
    * You have to also use ``sphinx.ext.todo`` extension

.. _Google Python Style Guide:
   http://google.github.io/styleguide/pyguide.html

"""

import os
import nibabel as nib
import pandas as pd


def read_electrode_snaplabel(eleclabels):

    # labels for the columns in the file. Current column labels are
    # contained in the header, we could write the function such that it
    # pulls  them out rather than hard-coding them in.
    varnames = ['IDX', '-R-', '-G-', '-B-', '-A--', 'VIS', 'MSH', 'LABEL']

    # open up labels file and extract pandas table
    with open(eleclabels) as f:
        label_table = pd.read_csv(f, delim_whitespace=True,lineterminator='\n',
                             comment='#',names=varnames)

    return label_table


def intensity2label(intersect,label_table):

    for m,n in intersect:
        print(label_table['LABEL'][n])

    return

def imgdiff(img1,img2):

    img1 = nib.load('C:\\Users\\tca11\\Desktop\\CNT\\Coregistration'
                   '\\resection_imgs\\imaging_campbell_mostRecents'
                   '\\RID142_6mth_postablation_MR'
                   '\\T00_RID142_ablation_roi_old.nii.gz')
    img2 = nib.load('C:\\Users\\tca11\\Desktop\\CNT\\Coregistration'
                    '\\resection_imgs\\imaging_campbell_mostRecents'
                    '\\RID142_6mth_postablation_MR'
                    '\\T00_RID142_ablation_roi_new2.nii.gz')
    diff = img2.get_fdata() - img1.get_fdata() # take difference
    diff[diff < 0] = 0
    diff_img = nib.Nifti1Image(diff, img2.affine, img2.header)
    nib.save(diff_img, 'C:\\Users\\tca11\\Desktop\\CNT\\Coregistration'
                   '\\resection_imgs\\imaging_campbell_mostRecents'
                   '\\RID142_6mth_postablation_MR'
                   '\\T00_RID142_ablation_roi2.nii.gz')

    return posdiff

def roi_intersect(roi1, roi2):
    """

    Function parameters should be documented in the ``Args`` section. The name
    of each parameter is required. The type and description of each parameter
    is optional, but should be included if not obvious.

    If \*args or \*\*kwargs are accepted,
    they should be listed as ``*args`` and ``**kwargs``.

    The format for a parameter is::

        name (type): description
            The description may span multiple lines. Following
            lines should be indented. The "(type)" is optional.

            Multiple paragraphs are supported in parameter
            descriptions.

    Args:
        roifile (str): points to a nifti file of the resection/ablation
        elecfile (str): points to a nifti file of the electrode ROIs
        eleclabels (str): point to a txt file corresponding elecfile
            intensity values with clinical labels

    Returns:
        bool: True if successful, False otherwise.

        The return type is optional and may be specified at the beginning of
        the ``Returns`` section followed by a colon.

        The ``Returns`` section may span multiple lines and paragraphs.
        Following lines should be indented to match the first line.

        The ``Returns`` section supports any reStructuredText formatting,
        including literal blocks::

            {
                'param1': param1,
                'param2': param2
            }

    Raises:
        AttributeError: The ``Raises`` section is a list of all exceptions
            that are relevant to the interface.
        ValueError: If `param2` is equal to `param1`.

    """

    # reshape matrices into arrays
    roi1 = roi1.reshape(-1)
    roi2 = roi2.reshape(-1)

    # get set of all unique image intensity pairs
    pairs = set(tuple(sorted([m, n])) for m, n in zip(roi1, roi2)) #

    # add non-zero pairs to an output variable
    intersect = []
    for m, n in pairs:
        if m != 0:
            intersect.append((m, n))

    return intersect


