import logging
import os
import tempfile
import unittest
import numpy as np
import SimpleITK as sitk
import warnings
import importlib

import disptools.predict as predict

class Test_Predict(unittest.TestCase):

    def __init__(self, *args, **kwargs):
            super(Test_Predict, self).__init__(*args, **kwargs)
            # ITK messes with the warnings...
            importlib.reload(warnings)


if __name__ == '__main__':
    unittest.main()
