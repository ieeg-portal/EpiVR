""" Unit tests for resection_util.py
"""

from settings import *
import resection_util
import pytest

# Test working example of load_image_manifest
result_ris = resection_util.load_image_manifest('HUP064')
assert len(result_ris.radiology_images) == len(
    open(os.path.join(DATA_DIR, 'HUP064', 'image_manifest.csv'), 'r').readlines())

# Test incorrect image manifest setup
with pytest.raises(IOError):
    result_ris = resection_util.load_image_manifest('NOT_A_PATIENT')
