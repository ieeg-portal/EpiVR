import os
import nibabel as nib
import ElectrodeLocalization as el
import pandas as pd

roifile = 'C:\\Users\\tca11\\Desktop\\CNT\\Coregistration' \
          '\\test_roi_intersect\\T00_RID274_ablation_roi.nii.gz'
elecfile = 'C:\\Users\\tca11\\Desktop\\CNT\\Coregistration' \
          '\\test_roi_intersect\\T00_RID274_mprageelectrodelabels_spheres' \
          '.nii.gz'
eleclabels = 'C:\\Users\\tca11\\Desktop\\CNT\\Coregistration' \
          '\\test_roi_intersect\\electrode_snaplabel.txt'

# load in the ablation ROI
img = nib.load(roifile)
roi1 = img.get_fdata()

# load in the electrode ROIs
img = nib.load(elecfile)
roi2 = img.get_fdata()


intersect = el.roi_intersect(roi1, roi2)
label_table = el.read_electrode_snaplabel(eleclabels)
el.intensity2label(intersect,label_table)

maindir = 'C:\\Users\\tca11\\Desktop\\CNT\\Coregistration\\resection_imgs' \
          '\\imaging_campbell_mostRecents'

subdirs = os.listdir(maindir)

for sub in subdirs:
    try:
        roifile = os.path.join(maindir,sub,'T00_' + sub[:6]
                               +'_ablation_roi.nii.gz')
        elecfile = os.path.join(maindir, sub,'T00_' + sub[:6] +'_mprageelectrodelabels_spheres.nii.gz')
        eleclabels = os.path.join(maindir,sub,'electrode_snaplabel.txt')

        # load in the ablation ROI
        img = nib.load(roifile)
        roi1 = img.get_fdata()

        # load in the electrode ROIs
        img = nib.load(elecfile)
        roi2 = img.get_fdata()

        intersect = el.roi_intersect(roi1, roi2)
        label_table = el.read_electrode_snaplabel(eleclabels)
        print(sub)
        el.intensity2label(intersect, label_table)
    except:
        print(roifile)
        print(elecfile)
        print(eleclabels)

