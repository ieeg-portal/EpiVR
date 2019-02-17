# -*- coding: utf-8 -*-
# Copyright (c) 2019 Center for Neuroengineering and Therapeutics.
# Copyright (c) 2019 University of Pennsylvania.
# License: MIT, see https://mit-license.org/
# Contact: lkini@pennmedicine.upenn.edu

"""RadiologicImageSet class

RadiologicImageSet is an abstract class designed to be wrap around multiple instances of RadiologicImage instances.
It provides a unified class structure that allows easy manipulation across a group of radiological images. For instance,
we should be easily able to filter a group radiology image instances based on MRI contrast, date of acquisition, etc.

Example:
    import RadiologicImageSet as ris
    image1 = ri.RadiologicImage("path/to/mri1.nii.gz", "MRI", "T1", "20000123", ["x","y","z"], [2,2,2])
    image2 = ri.RadiologicImage("path/to/mri2.nii.gz", "MRI", "T2", "20000123", ["x","y","z"], [0.5,0.5,4])
    image3 = ri.RadiologicImage("path/to/mri3.nii.gz", "MRI", "T1+GAD", "20000123", ["x","y","z"], [1,1,1.5])
    ex = ris.RadiologicImageSet([image1, image2, image3])

Attributes:
    radiologic_image_list (list): List of RadiologicImage instances or just one RadiologicImage instance
"""
from RadiologicImage import RadiologicImage


class RadiologicImageSet:
    def __init__(self, radiologic_image_list):
        # enforce entry standards
        if type(radiologic_image_list) is not list and type(radiologic_image_list) is not RadiologicImage:
            raise ValueError('%s is not a valid set of radiology images.')

        if type(radiologic_image_list) is RadiologicImage:
            self.radiology_images = [radiologic_image_list]  # Only one RadiologicImage passed
        else:
            self.radiology_images = radiologic_image_list

    def copy(self):
        """
        create a new instance of RadiologicImageSet, with the same data as this instance.
        """
        return RadiologicImageSet(self.radiology_images)

    def filter_by_radiology_type(self, radiology_type):
        """This module level function filters by type of image (e.g. MRI, CT, PET).

        Args:
            radiology_type (str): Type of radiology image passed as string.
                Type can only be one of the radiology types in RadiologicImage class.

        Returns:
            RadiologicImageSet instance with only the appropriate RadiologyImages

        Raises:
            TypeError: If radiology_type is not one of the RadiologicImage TYPES

        """
        if radiology_type not in RadiologicImage.TYPE:
            raise TypeError('Radiology type %s is a not a proper image type.' % radiology_type)
        filtered_ris = self.copy()
        filtered_ris.radiology_images = [radiology_image for radiology_image in self.radiology_images if
                                         radiology_image.radiologic_type == radiology_type]
        return filtered_ris

    def filter_by_contrast(self, contrast):
        """This module level function filters by contrast (both intrinsic and extrinsic).

        For instance, for MRIs can support filtering by T1-weighted, T2-weighted, FLAIR, as well as Gadolinium enhanced
            images.

        Args:
            contrast (str): Type of contrast passed as string.

        Returns:
            RadiologicImageSet instance with only the appropriate RadiologyImages

        Raises:
            TypeError: If contrast is not a string

        """
        if type(contrast) != str:
            raise TypeError('Image Contrast passed to filter_by_contrast is not a string.')
        filtered_ris = self.copy()
        filtered_ris.radiology_images = [radiology_image for radiology_image in self.radiology_images if
                                         radiology_image.contrast == contrast]
        return filtered_ris

    def filter_by_acquisition_date(self, acquisition_date):
        """This module level function filters by date of acquisition.

        This can be useful to get all radiology images acquired during a study imaging session.

        Args:
            acquisition_date (str): Date of acquisition in the form "YYYYMMDD".

        Returns:
            RadiologicImageSet instance with only the appropriate RadiologyImages

        Raises:
            TypeError: If acquisition_date is not a string.
            ValueError: If acquisition_date is not the appropriate format.

        """
        if type(acquisition_date) != str:
            raise TypeError('Acquisition date passed to filter_by_acquisition_date is not a string.')
        if len(acquisition_date) != 8 or not acquisition_date.isdigit():
            raise ValueError(
                'Acquisition date passed to filter_by_acquisition_date should follow YYYYMMDD date format.')

        filtered_ris = self.copy()
        filtered_ris.radiology_images = [radiology_image for radiology_image in self.radiology_images if
                                         radiology_image.acquisition_date == acquisition_date]
        return filtered_ris

    def filter_by_event(self):
        pass
