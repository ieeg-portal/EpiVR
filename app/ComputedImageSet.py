# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: lkini@pennmedicine.upenn.edu

"""ComputedImageSet class

ComputedImageSet is an abstract class designed to be wrap around multiple instances of ComputedImage instances.
It provides a unified class structure that allows easy manipulation across a group of radiological images. For instance,
we should be easily able to filter a group radiology image instances based on MRI contrast, date of acquisition, etc.

Example:
    import RadiologicImage as ri
    import ComputedImage as ci
    import ComputedImageSet as cis
    ref_image1 = ri.RadiologicImage("path/to/mri1.nii.gz", "MRI", "T1", "20000123", ["x","y","z"], [2,2,2])
    comp_image1 = ci.ComputedImage("path/to/gray_segmentation.nii.gz", ref_image1, "SegmentationImage")
    comp_image2 = ci.ComputedImage("path/to/white_segmentation.nii.gz", ref_image1, "SegmentationImage")
    ex = cis.ComputedImageSet([comp_image1, comp_image2])

"""
from ComputedImage import ComputedImage


class ComputedImageSet:
    def __init__(self, computed_image_list):
        # enforce entry standards
        if type(computed_image_list) is not list and not isinstance(computed_image_list, ComputedImage):
            raise ValueError('%s is not a valid set of radiology images.')

        if isinstance(computed_image_list, ComputedImage):
            self.computed_images = [computed_image_list]  # Only one ComputedImage passed
        else:
            self.computed_images = computed_image_list

    def copy(self):
        """
        create a new instance of ComputedImageSet, with the same data as this instance.
        """
        return ComputedImageSet(self.computed_images)

    def filter_by_computed_type(self, computed_type):
        """This module level function filters by type of computation (e.g. Segmentations vs Registrations).

        Args:
            computed_type (str): Type of computational image passed as string.
                Type can only be one of the computational types in ComputedImage class.

        Returns:
            ComputedImageSet instance with only the appropriate ComputedImages

        Raises:
            TypeError: If computed_type is not one of the ComputedImage TYPES

        """
        if computed_type not in ComputedImage.COMP_TYPE:
            raise TypeError('Computated type %s is a not a proper image type.' % computed_type)
        filtered_cis = self.copy()
        filtered_cis.computed_images = [computed_image for computed_image in self.computed_images if
                                        computed_image.computed_type == computed_type]
        return filtered_cis

    def filter_by_radiology_type(self, radiology_type):
        """This module level function filters by type of image (e.g. MRI, CT, PET).

        Args:
            radiology_type (str): Type of radiology image passed as string.
                Type can only be one of the radiology types in RadiologicImage class.

        Returns:
            ComputedImageSet instance with only the appropriate ComputedImages

        Raises:
            TypeError: If radiology_type is not one of the RadiologicImage TYPES

        """
        if radiology_type not in ComputedImage.TYPE:
            raise TypeError('Radiology type %s is a not a proper image type.' % radiology_type)
        filtered_cis = self.copy()
        filtered_cis.computed_images = [computed_image for computed_image in self.computed_images if
                                        computed_image.radiologic_type == radiology_type]
        return filtered_cis

    def filter_by_contrast(self, contrast):
        """This module level function filters by contrast (both intrinsic and extrinsic).

        For instance, for MRIs can support filtering by T1-weighted, T2-weighted, FLAIR, as well as Gadolinium enhanced
            images.

        Args:
            contrast (str): Type of contrast passed as string.

        Returns:
            ComputedImageSet instance with only the appropriate RadiologyImages

        Raises:
            TypeError: If contrast is not a string

        """
        if type(contrast) != str:
            raise TypeError('Image Contrast passed to filter_by_contrast is not a string.')
        filtered_cis = self.copy()
        filtered_cis.computed_images = [computed_image for computed_image in self.computed_images if
                                        computed_image.contrast == contrast]
        return filtered_cis

    def filter_by_computation_property(self, property_name, property_value):
        pass

    def filter_by_patient(self, patient_id):
        """

        :param patient_id:
        :return:
        """
        pass

    def filter_by_event(self, event_type):
        """

        :param event_type:
        :return:
        """
        pass
