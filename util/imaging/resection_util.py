# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: lkini@pennmedicine.upenn.edu

"""Example Google style docstrings. A docstring of this format should be placed
at the top of all python scripts.

This module demonstrates documentation as specified by the `Google Python
Style Guide`_. Docstrings may extend over multiple lines. Sections are created
with a section header and a colon followed by a block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

"""

import os
import nibabel as nib

from app.RadiologicImage import RadiologicImage as ri
from app.RadiologicImageSet import RadiologicImageSet as ris

try:
    DATA_DIR = os.environ['DATA_DIR']
except KeyError:
    raise KeyError('DATA_DIR has not been properly setup. Did you setup .env file correctly?')


def load_image_manifest(patient_id):
    """This is an example of a module level function.

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
        patient_id (str): Patient ID.
        param2 (:obj:`str`, optional): The second parameter. Defaults to None.
            Second line of description should be indented.
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

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
    # Check patient ID
    if not os.path.isdir(os.path.join(DATA_DIR, patient_id)) or not os.path.isfile(
            os.path.join(DATA_DIR, patient_id, 'image_manifest.csv')):
        raise IOError(
            'Patient ID %s does not have an appropriately set up data directory in %s' % (patient_id, DATA_DIR))

    # Set image manifest path
    image_manifest_path = os.path.join(DATA_DIR, patient_id, 'image_manifest.csv')
    image_manifest = open(image_manifest_path, 'r').readlines()

    # Create list of all radiologic images specified in image manifest file
    ri_list = []
    for image_line in image_manifest:
        date, modality, modality_type, filepath = image_line.replace('\n', '').replace('\r', '').split(',')

        # Check filepath
        if not os.path.isfile(filepath):
            raise IOError('File %s does not exist.' % filepath)

        # Check file type
        if filepath.endswith('.nii') or filepath.endswith('.nii.gz'):
            img = nib.load(filepath)
            img_hdr = img.get_header()
            img_resolution = img_hdr.get_zooms()
            if len(img_resolution) == 3:
                img_dimensions = ["x", "y", "z'"]
            elif len(img_resolution) == 2:
                img_dimensions = ["x", "y"]
            else:
                raise NotImplementedError(
                    'Handling files with more than 3 dimensions like %s has not been implemented yet.' % os.path.basename(
                        filepath))
        else:
            raise NotImplementedError(
                'Handling files that end with the extension in filename %s has not been implemented yet.' % os.path.basename(
                    filepath))

        # Initialize RadiologicImage
        radiologic_image = ri(filepath, modality, modality_type, date, img_dimensions, img_resolution)
        ri_list.append(radiologic_image)

    # Generate the RadiologicImageSet
    radiologic_image_set = ris(ri_list)

    return radiologic_image_set
