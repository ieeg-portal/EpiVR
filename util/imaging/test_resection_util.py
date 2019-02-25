""" Unit tests for resection_util.py
"""

from settings import *
import resection_util

from app.RadiologicImageSet import RadiologicImageSet as ris

result_ris = resection_util.load_image_manifest('HUP064')
assert len(result_ris.radiology_images) == len(open(os.path.join(DATA_DIR, 'HUP064', 'image_manifest.csv'), 'r').readlines())
