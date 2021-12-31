import importlib
import sys
import warnings

import unittest

import numpy as np
import SimpleITK as sitk

import disptools.drawing as drawing
import disptools.simulatrophy as simulatrophy

class Test_Simulatrophy(unittest.TestCase):

    def __init__(self, *args, **kwargs):
            super(Test_Simulatrophy, self).__init__(*args, **kwargs)
            # ITK messes with the warnings...
            importlib.reload(warnings)

    def test_mask_to_simulatrophy_mask_radius_none(self):
        size = 11
        input_mask = np.zeros([size, size, size]).astype(drawing.np_float_type)
        input_mask[4:7, 4:7, 4:7] = 1.0  # 3x3x3 ROI centred in the volume
        radius = None
        input_image = sitk.GetImageFromArray(input_mask)
        oracle = 2 * input_mask + 1  # ROI = 3 (WM), background = 1 (CSF)

        result = simulatrophy.mask_to_simulatrophy_mask(input_image, radius)

        result_array = sitk.GetArrayFromImage(result)
        self.assertTrue(np.array_equal(result_array, oracle))


if __name__ == '__main__':
    unittest.main()
