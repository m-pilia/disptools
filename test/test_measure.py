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
    import skimage.draw
except ImportError as e:
    print("Warning: cannot import 'skimage.draw' module. " +
          "Some functionalities depending upon it may be unavailable.")

import disptools.measure as measure
import disptools.drawing as drawing

class Test_Measure(unittest.TestCase):

    def __init__(self, *args, **kwargs):
            super(Test_Measure, self).__init__(*args, **kwargs)
            # ITK messes with the warnings...
            importlib.reload(warnings)

    def test_reject_outliers(self):

        n = 200
        k = 2.0
        mu = 0.0
        sigma = 1.0

        data = np.random.normal(mu, sigma, size=n)
        deviations = data - np.median(data)
        deviations = np.abs(deviations)
        mad = np.median(deviations)
        ind = deviations < mad * k

        result = measure.reject_outliers(data, k)

        self.assertTrue(np.array_equal(np.sort(data[ind]),
                                       np.sort(result)))

    @unittest.skipIf('skimage.draw' is sys.modules, "skimage.draw required for this feature")
    def test_cubicity(self):

        n = 100
        r = 20

        # cube, uniform grid
        data = np.zeros((n, n, n), dtype=np.uint8)
        data[n//2-r:n//2+r, n//2-r:n//2+r, n//2-r:n//2+r] = 1
        image = sitk.GetImageFromArray(data)

        result = measure.cubicity(image)

        self.assertEqual(1.0, result)

        # cube, anisotropic grid
        data = np.zeros((n, n, n), dtype=np.uint8)
        data[n//2-2*r:n//2+2*r, n//2-r:n//2+r, n//2-r:n//2+r] = 1
        image = sitk.GetImageFromArray(data)
        image.SetSpacing((1.0, 1.0, 0.5))

        result = measure.cubicity(image)

        self.assertEqual(1.0, result)

        # sphere
        data = skimage.draw.ellipsoid(r, r, r).astype(np.uint8)
        image = sitk.GetImageFromArray(data)

        result = measure.cubicity(image)

        self.assertTrue(result < 0.51)

    @unittest.skipIf('skimage.draw' is sys.modules, "skimage.draw required for this feature")
    def test_sphericity(self):

        n = 100
        r = 20

        data = np.zeros((n, n, n), dtype=np.uint8)
        data[n//2-r:n//2+r, n//2-r:n//2+r, n//2-r:n//2+r] = 1
        image = sitk.GetImageFromArray(data)

        result = measure.sphericity(image)

        self.assertTrue(result < 0.91)

        data = skimage.draw.ellipsoid(r, r, r).astype(np.uint8)
        image = sitk.GetImageFromArray(data)

        result = measure.sphericity(image)

        self.assertTrue(result > 0.99)

    @unittest.skipIf('skimage.draw' is sys.modules, "skimage.draw required for this feature")
    def test_average_jacobian_error(self):

        n = 100
        r = 20

        data_arrays = [
            np.ones((n, n, n), dtype=measure.np_float_type),
            skimage.draw.ellipsoid(r, r, r).astype(measure.np_float_type),
            1.0 + skimage.draw.ellipsoid(r, r, r).astype(measure.np_float_type),
            1.0 + 2.0 * skimage.draw.ellipsoid(r, r, r).astype(measure.np_float_type)
        ]
        masks = [
            None,
            skimage.draw.ellipsoid(r, r, r).astype(measure.np_float_type),
            skimage.draw.ellipsoid(r, r, r).astype(measure.np_float_type),
            skimage.draw.ellipsoid(r, r, r).astype(measure.np_float_type)
        ]
        oracles = [
            0.0,
            0.0,
            1.0,
            2.0
        ]

        item = 0
        for data, mask_data, oracle in zip(data_arrays, masks, oracles):
            item += 1

            with self.subTest(item=item):

                image = sitk.GetImageFromArray(data)
                mask = sitk.GetImageFromArray(mask_data) if mask_data is not None else None

                result = measure.average_jacobian_error(image, mask)

                self.assertEqual(oracle, result)

    @unittest.skipIf('skimage.draw' is sys.modules, "skimage.draw required for this feature")
    def test_minkowski_compactness(self):

        a = np.zeros((300, 300, 300), dtype=np.uint8)
        a[10:290,10:290,10:290] = 1

        cube = sitk.GetImageFromArray(a)
        sphere = sitk.GetImageFromArray(skimage.draw.ellipsoid(100, 100, 100).astype(np.uint8)) > 0
        diamond = drawing.create_sphere(40, size=100, centre=(50, 50, 50), norm=1) > 0

        self.assertTrue(measure.minkowski_compactness(cube, norm='inf') > 0.99)
        self.assertTrue(measure.minkowski_compactness(sphere, norm='inf') < 0.97)
        self.assertTrue(measure.minkowski_compactness(diamond, norm='inf') < 0.91)

        self.assertTrue(measure.minkowski_compactness(cube, norm=1) < 0.91)
        self.assertTrue(measure.minkowski_compactness(sphere, norm=1) < 0.98)
        self.assertTrue(measure.minkowski_compactness(diamond, norm=1) > 0.99)


if __name__ == '__main__':
    unittest.main()
