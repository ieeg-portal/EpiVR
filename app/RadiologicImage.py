# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: tcarnold@seas.upenn.edu

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
    acquisition_date (str): Date should come from image_manifest.json
        (Ex: "20000123")
    image_dimensions (list of str): Dimensions in human interpretable form
        (Ex: ["x","y","z"], ["x","y","time"], etc.)
    image_resolution (int): Scale of each corresponding dimension
        (Ex: [2,2,2])
    image_units (list): Units of each corresponding dimension
        (Ex: ['mm', 'cm', 'mm'])

Todo:
    * image resolution: should we mandate a scale, such as measurements must be
     in mm? (think about PET being on cm scale)
    * get date of acquisition from from image manifest
    * Consider using datetime instead of str for date of acquisition

Alterations:
    TCA 1/23/19 - initialized class & wrote unit test
    TCA 1/25/19 - updated unit test to pytest format, added TYPE restrictions
"""


class RadiologicImage:
    # standards for data
    TYPE = ('MRI', 'CT', 'PET', 'SPECT', 'XRAY')
    UNITS = ('mm', 'cm', 'm', 'in', 'ft', 's', 'sec', 'min', 'hr')

    def __init__(self, radiologic_type, contrast, acquisition_date,
                 image_dimensions, image_resolution, image_units=None):
        """

        :param radiologic_type:
        :param contrast:
        :param acquisition_date:
        :param image_dimensions:
        :param image_resolution:
        :param image_units:
        """

        # enforce entry standards
        if radiologic_type not in self.TYPE:
            raise TypeError('%s is not a valid image type. Valid Types: %s' %
                            (radiologic_type, self.TYPE))

        # apply input values to attributes of class instance
        self.radiologic_type = radiologic_type
        self.contrast = contrast
        self.acquisition_date = acquisition_date
        self.image_dimensions = image_dimensions
        self.image_resolution = image_resolution

        # Autopopulate to mm units if not passed in during instantiation
        if image_units is None:
            image_units = ['mm'] * len(self.image_dimensions)
        # Check image units meet standards and save property to instance
        if len(image_units) != len(self.image_dimensions):
            raise TypeError('The length of image dimension units should match the number of dimensions.')
        for image_unit in image_units:
            if image_unit not in self.UNITS:
                raise TypeError('The unit %s is not a proper unit of measurement.' % image_unit)
        self.image_units = image_units
