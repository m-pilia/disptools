from typing import *
import logging
import sys
import numpy as np
import SimpleITK as sitk
import platform
import os
import site
from glob import glob

import _disptools

try:
    import itk
except ImportError as e:
    # ITK is an optional dependency. Do not raise a warning here.
    pass

__author__    = 'Martino Pilia'
__copyright__ = 'Copyright (C) 2018 Martino Pilia'
__license__   = 'MIT'
__version__   = '0.1'

# Set the floating point type
# OBS! This must be the same used in the C library
if _disptools.get_float_type_size() == 32:
    np_float_type = np.float32
    sitk_float_type = sitk.sitkFloat32
    sitk_vector_float_type = sitk.sitkVectorFloat32
    itk_float_type = itk.F if 'itk' in sys.modules else None
elif _disptools.get_float_type_size() == 64:
    np_float_type = np.float64
    sitk_float_type = sitk.sitkFloat64
    sitk_vector_float_type = sitk.sitkVectorFloat64
    itk_float_type = itk.D if 'itk' in sys.modules else None
else:
    logging.critical('Unsupported floating point size (%d bit)'
                     % _disptools.get_float_type_size())
    raise Exception('Unsupported floating point type size')

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

