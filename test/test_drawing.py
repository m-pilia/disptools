import logging
import os
import sys
import tempfile
import unittest
import numpy as np
import SimpleITK as sitk
import warnings
import importlib

try:
    import itk
except ImportError as e:
    print("Warning: cannot import 'itk' module. " +
          "Some functionalities depending upon it may be unavailable.")
try:
    import skimage.draw
except ImportError as e:
    print("Warning: cannot import 'skimage.draw' module. " +
          "Some functionalities depending upon it may be unavailable.")

import disptools.drawing  as drawing

class Test_Drawing(unittest.TestCase):

    def __init__(self, *args, **kwargs):
            super(Test_Drawing, self).__init__(*args, **kwargs)
            # ITK messes with the warnings...
            importlib.reload(warnings)

    def test_resize_image(self):

        size_1 = (100, 100, 100)
        size_2 = (125, 125, 125)

        img = sitk.Image(size_1, drawing.sitk_float_type)

        res = drawing.resize_image(img, size_2)

        self.assertTrue(res.GetSize() == size_2)

    def test_polar_conversion(self):

        c = np.array([10,10,10])
        r = 10
        norm = 2

        points = [
            c,
            c + np.array([r, 0, 0]),
            c + np.array([0, r, 0]),
            c + np.array([0, 0, r]),
            c - np.array([r, 0, 0]),
            c - np.array([0, r, 0]),
            c - np.array([0, 0, r])
        ]

        oracles = [
            (0, 0,       0       ),
            (1, np.pi/2, 0       ),
            (1, np.pi/2, np.pi/2 ),
            (1, 0,       0       ),
            (1, np.pi/2, np.pi   ),
            (1, np.pi/2, -np.pi/2),
            (1, np.pi,   0       )
        ]

        for point, oracle in zip(points, oracles):
            with self.subTest(point=point):
                result = drawing.polar_conversion(norm, c, r, *point)

                self.assertTrue(np.allclose(result, oracle))

    @unittest.skipIf('skimage.draw' is sys.modules, "skimage.draw required for this feature")
    def test_create_sphere(self):

        r = 20
        sphere = drawing.create_sphere(r)
        data = sitk.GetArrayViewFromImage(sphere)
        oracle = skimage.draw.ellipsoid(r, r, r).astype(drawing.np_float_type)

        self.assertTrue(np.array_equal(data, oracle))

        r = 20
        size = 50
        c = size // 2
        oracle = np.zeros((size, size, size), dtype=drawing.np_float_type)
        oracle[c-r:c+r+1, c-r:c+r+1, c-r:c+r+1] = 1.0
        sphere = drawing.create_sphere(r, size=size, centre=(c, c, c), norm='max')
        data = sitk.GetArrayViewFromImage(sphere)

        self.assertTrue(np.array_equal(data, oracle))

    def test_mask(self):

        n = 50
        data = np.random.rand(n, n, n)
        index = np.random.rand(n, n, n) < 0.2

        image = sitk.GetImageFromArray(data)
        mask = sitk.GetImageFromArray(index.astype(drawing.np_float_type))

        res = drawing.mask(image, mask)

        a = sitk.GetArrayViewFromImage(res)

        self.assertFalse(np.any(a[np.logical_not(index)]))
        self.assertTrue(np.array_equal(data[index], a[index]))

        res = drawing.mask(image, mask, jacobian=True)

        a = sitk.GetArrayViewFromImage(res)

        self.assertFalse(np.any(a[np.logical_not(index)] - 1.0))
        self.assertTrue(np.array_equal(data[index], a[index]))

    def test_float_dilate(self):

        n = 100
        dilations = list(range(1, 5))

        for dilation in dilations:
            with self.subTest(dilation=dilation):
                data = (np.random.rand(n, n, n) < 0.5).astype(np.uint8)

                binary_image = sitk.GetImageFromArray(data)
                dilated = sitk.BinaryDilate(binary_image, dilation)

                image = sitk.Cast(binary_image, drawing.sitk_float_type)
                oracle = sitk.Cast(dilated, drawing.sitk_float_type)

                result = drawing.float_dilate(image, dilation)

                self.assertTrue(np.allclose(sitk.GetArrayViewFromImage(result),
                                            sitk.GetArrayViewFromImage(oracle)))

    @unittest.skipIf('itk' not in sys.modules, "itk module required for this feature")
    def test_sitk_to_itk(self):

        n = 100
        spacing = (1.0, 2.0, 3.0)
        origin = (2.0, 3.0, 4.0)
        size = (n, n, n)
        data_arrays = [
            np.random.rand(*size).astype(drawing.np_float_type),
            np.random.rand(*size, 3).astype(drawing.np_float_type),
            np.random.rand(*np.random.randint(20, 150, size=3), 3).astype(drawing.np_float_type)
        ]

        for data in data_arrays:
            with self.subTest(shape=data.shape):

                image_sitk = sitk.GetImageFromArray(data)
                image_sitk.SetSpacing(spacing)
                image_sitk.SetOrigin(origin)

                result = drawing.sitk_to_itk(image_sitk)

                self.assertTrue(result.GetOrigin() == origin)
                self.assertTrue(result.GetSpacing() == spacing)
                self.assertTrue(np.array_equal(itk.GetArrayViewFromImage(result),
                                               data))


if __name__ == '__main__':
    unittest.main()
