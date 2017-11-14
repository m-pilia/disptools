import logging
import os
import tempfile
import unittest
import numpy as np
import SimpleITK as sitk
import warnings
import importlib

import disptools.displacements as displacements
import disptools.drawing as drawing

class Test_Displacements(unittest.TestCase):

    def __init__(self, *args, **kwargs):
            super(Test_Displacements, self).__init__(*args, **kwargs)
            # ITK messes with the warnings...
            importlib.reload(warnings)

    def test_regularise(self):

        n = 100
        epsilon = 0.05
        original_data = np.random.rand(n, n, n).astype(drawing.np_float_type) + epsilon

        mask = np.random.randint(5, size=(n, n, n)) % 4 == 0
        not_mask = np.logical_not(mask)
        perturbed_data = np.copy(original_data)
        perturbed_data[mask] = epsilon - np.random.rand()

        self.assertTrue(np.all(original_data >= epsilon))
        self.assertTrue(np.any(perturbed_data[mask] < epsilon))
        self.assertTrue(np.array_equal(original_data[not_mask], perturbed_data[not_mask]))

        img = sitk.GetImageFromArray(perturbed_data)

        res = displacements.regularise(img)

        a = sitk.GetArrayViewFromImage(res)
        self.assertTrue(np.array_equal(original_data[not_mask], a[not_mask]))
        self.assertTrue(np.all(a[mask] == epsilon))

    def test_jacobian(self):

        n = 10
        img = drawing.sin_vector_field(n)

        jacobian_pw = displacements.jacobian(img)
        jacobian_itk = sitk.DisplacementFieldJacobianDeterminant(img)

        # On the border, ITK is behaving weird
        jacobian = sitk.GetArrayViewFromImage(jacobian_pw)[1:n-1, 1:n-1, 1:n-1]
        oracle = sitk.GetArrayViewFromImage(jacobian_itk)[1:n-1, 1:n-1, 1:n-1]

        self.assertTrue(np.allclose(jacobian, oracle, 1e-3))

    def test_redistribute_volume_change(self):

        size = (100, 100, 100)
        data = np.random.rand(*size) * 2.0
        index = np.random.rand(*size) < 0.2

        img = sitk.GetImageFromArray(data)
        mask = sitk.GetImageFromArray(index.astype(displacements.np_float_type))

        result = displacements.redistribute_volume_change(img, mask)

        a = sitk.GetArrayViewFromImage(result)

        self.assertTrue(np.array_equal(a[index], data[index]))
        self.assertTrue(np.sum(a) - a.size < 1e-5)

    def test_average_jacobian(self):

        n = 10
        d = 100
        size = (d, d, d)
        spacing = (1.0, 2.0, 3.0)
        ext = 'nii'

        with tempfile.TemporaryDirectory() as tmp_dir:

            avg = np.ones(size, dtype=displacements.np_float_type)

            for i in range(0, n):
                data = np.random.rand(*size) + 0.5
                avg *= data
                img = sitk.GetImageFromArray(data)
                img.SetSpacing(spacing)
                sitk.WriteImage(img, "%s/%d.%s" % (tmp_dir, i, ext))

            oracle = avg ** (1.0 / n)

            average_jacobian = displacements.average_jacobian("%s/*.%s" % (tmp_dir, ext))

        result_data = sitk.GetArrayFromImage(average_jacobian)

        self.assertTrue(np.allclose(oracle, result_data))
        self.assertTrue(average_jacobian.GetSpacing() == spacing)

    def test_deformation_to_displacement(self):

        n = 50
        deformation = np.empty((n, n, n, 3))
        for z, y, x in np.ndindex((n, n, n)):
            deformation[z, y, x, 0] = x
            deformation[z, y, x, 1] = y
            deformation[z, y, x, 2] = z
        field = sitk.GetImageFromArray(deformation)

        res = displacements.deformation_to_displacement(field)

        self.assertFalse(np.any(sitk.GetArrayViewFromImage(res)))

        displacement = 2 * np.random.rand(n, n, n, 3) - 1.0
        deformation = deformation + displacement
        field = sitk.GetImageFromArray(deformation)

        res = displacements.deformation_to_displacement(field)

        self.assertTrue(np.allclose(sitk.GetArrayFromImage(res), displacement))

    def test_compose_displacements(self):

        # Test composition in a inverse consistency setting

        n = 100
        absolute_tolerance = 5e-2
        displacement = drawing.sin_vector_field(n)
        inverse_displacement = sitk.InvertDisplacementField(displacement)

        field = displacements.compose_displacements(displacement, inverse_displacement)

        image = sitk.GetImageFromArray(np.random.rand(n, n, n))
        warped = sitk.Warp(image, field)
        control = sitk.Warp(image, displacement)

        # Crop two voxels along the boundary in order to ignore boundary
        # effects due to sitk.InvertDisplacementField
        a_image = sitk.GetArrayViewFromImage(image)[3:n-2, 3:n-2, 3:n-2]
        a_warped = sitk.GetArrayViewFromImage(warped)[3:n-2, 3:n-2, 3:n-2]
        a_control = sitk.GetArrayViewFromImage(control)[3:n-2, 3:n-2, 3:n-2]

        self.assertTrue(np.allclose(a_image, a_warped, atol=absolute_tolerance))
        self.assertFalse(np.allclose(a_image, a_control, atol=absolute_tolerance))

    def test_field_zero_padding(self):

        n = 20
        d = 2
        data = np.random.rand(n, n, n, 3)
        image = sitk.GetImageFromArray(data)

        for size_x in np.ndindex(d, d):
            for size_y in np.ndindex(d, d):
                for size_z in np.ndindex(d, d):
                    with self.subTest(size=(size_x, size_y, size_z)):

                        result = displacements.field_zero_padding(image, size_x, size_y, size_z)

                        nx = size_x[0]
                        ny = size_y[0]
                        nz = size_z[0]
                        a = sitk.GetArrayFromImage(result)

                        self.assertTrue(np.array_equal(a[nx:nx+n,ny:ny+n,nz:nz+n,:], data))

                        a[nx:nx+n,ny:ny+n,nz:nz+n,:] = 0.0

                        self.assertFalse(np.any(a))

    def test_warp_points_by_displacement(self):

        n = 100
        m = 100
        data = 10.0 * np.random.rand(n, n, n, 3) - 5.0
        displacement = sitk.GetImageFromArray(data)
        points = np.random.randint(0, n-1, size=(m, 3))

        result = displacements.warp_points_by_displacement(points, displacement)

        oracle = np.empty(points.shape)

        for i in range(0, points.shape[0]):
            oracle[i,:] = points[i,:] + data[points[i,2], points[i,1], points[i,0], :]

        self.assertTrue(np.allclose(result, oracle))


if __name__ == '__main__':
    unittest.main()
