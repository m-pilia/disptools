import logging
import os
import tempfile
import unittest
import numpy as np
import SimpleITK as sitk
import warnings
import importlib

import disptools.simulatrophy as simulatrophy

class Test_Simulatrophy(unittest.TestCase):

    def __init__(self, *args, **kwargs):
            super(Test_Simulatrophy, self).__init__(*args, **kwargs)
            # ITK messes with the warnings...
            importlib.reload(warnings)


if __name__ == '__main__':
    unittest.main()
