# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: lkini@pennmedicine.upenn.edu

"""ComputedImage class

ComputedImage is an abstract class designed to abstract computational images such as segmentations and registrations.
It provides a unified class structure that allows easy manipulation across image pipelines.

Example:
    import RadiologicImage as ri
    import ComputedImage as ci
    ri = ri.RadioligcImage("/tmo/original_mri.nii", "MRI", "T1", "20000123", ["x","y","z"], [2,2,2])
    ex = ci.ComputedImage("/tmp/registered_image.nii", ri)
    ex.unitTest()

Attributes:
    reference_radiologic_image (RadiologicImage or str): Reference RadiologicImage instance or string file path
        to reference radiology image
    computed_type (str): Type of computation
        (Ex: "RegisteredImage","SegmentationImage","Tensor")
"""
import os
from RadiologicImage import RadiologicImage


class ComputedImage:
    # standards for data
    TYPE = ('MRI', 'CT', 'PET', 'SPECT', 'XRAY')
    COMP_TYPE = ('RegisteredImage', 'SegmentationImage', 'TensorImage', 'TemplateImage')
    UNITS = ('mm', 'cm', 'm', 'in', 'ft', 's', 'sec', 'min', 'hr')

    def __init__(self, image_location, reference_radiologic_image, computed_type):
        """

        :param reference_radiologic_image:
        :param computed_type:
        """

        # enforce entry standards
        if not os.path.exists(os.path.expanduser(image_location)) or not os.path.isfile(
                os.path.expanduser(image_location)):
            raise IOError('%s is not a valid image location.' % image_location)
        if not isinstance(reference_radiologic_image, RadiologicImage):
            raise TypeError('Input %s is not a valid reference radiology image.' % reference_radiologic_image)
        if computed_type not in self.COMP_TYPE:
            raise TypeError('%s is not a valid computational image type. Valid Types: %s' %
                            (computed_type, self.COMP_TYPE))

        # apply input values to attributes of class instance
        self.image_location = image_location
        self.reference_radiologic_image = reference_radiologic_image
        self.computed_type = computed_type
