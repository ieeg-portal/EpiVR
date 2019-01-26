# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/

"""RadiologicImage class

RadiologicImage is an abstract class designed to be used across multiple
imaging modalities. It provides a unified class structure that allows easy
manipulation across modalities.

Example:
    import RadiologicImage as ri
    ex = ri.RadiologicImage("MRI", "T1", "20000123", ["x","y","z"], [2,2,2])
    ex.unitTest()

Attributes:
    radiologic_type (str): The imaging modality for your data
        (Ex: "MRI","CT","PET")
    contrast (str): Pulse sequence
        (Ex: "T1","T2")
    date_of_acquisition (str): Date should come from image_manifest.json
        (Ex: "20000123")
    image_dimensions (str): Dimensions in human interpretable form
        (Ex: ["x","y","z"], ["x","y","time"], etc.)
    image_resolution (int): Scale of each correspond dimension, currently in mm
        (Ex: [2,2,2])

Todo:
    * image resolution: should we mandate a scale, such as measurements must be
     in mm? (think about PET being on cm scale)
    * get date of acquisition from from image manifest
    * consider limiting radiologic_type
    Ex: RADIOLOGICY_TYPE_OPTIONS = ['MRI','CT','PET','SPECT','XRAY']
    * Consider using datetime instead of str for date of acquisition

Alterations:
    TCA 1/23/19 - initialized class & wrote unit test
"""

import RadiologicImage as ri

ex = ri.RadiologicImage("MRI", "T1", "20000123", ["x", "y", "z"], [2, 2, 2])
for property, value in vars(ex).items():
    print(property, ": ", value)
