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
    import vtk
except ImportError as e:
    print("Warning: cannot import 'vtk' module. " +
          "Some functionalities depending upon it may be unavailable.")

import disptools.io as io

class Test_IO(unittest.TestCase):

    def __init__(self, *args, **kwargs):
            super(Test_IO, self).__init__(*args, **kwargs)
            # ITK messes with the warnings...
            importlib.reload(warnings)

    @unittest.skipIf(io.np_float_type != np.float64, "works only in double precision")
    def test_read_rvf(self):

        n = 100
        size = (n, n, n)
        spacing = (1.0, 2.0, 3.0)
        data = np.random.rand(n, n, n, 3)

        with tempfile.NamedTemporaryFile() as tmp_file:
            tmp_file.write(("%d %d %d\n%f %f %f\n" % (*size, *spacing)).encode('utf-8'))
            data.tofile(tmp_file)
            tmp_file.flush()

            field = io.read_rvf(tmp_file.name)

            self.assertTrue(np.array_equal(sitk.GetArrayViewFromImage(field), data))
            self.assertTrue(field.GetSpacing() == spacing)
            self.assertTrue(field.GetSize() == size)

    def test_make_unique_directory(self):

        with tempfile.TemporaryDirectory() as tmp_dir:

            name = os.path.join(tmp_dir, 'foo')

            self.assertFalse(os.path.isdir(name))

            result = io.make_unique_directory(name)

            self.assertTrue(name == result)
            self.assertTrue(os.path.isdir(result))

            result = io.make_unique_directory(os.path.join(tmp_dir, name))

            self.assertEqual(name + '_1', result)
            self.assertTrue(os.path.isdir(result))

            result = io.make_unique_directory(os.path.join(tmp_dir, name))

            self.assertEqual(name + '_2', result)
            self.assertTrue(os.path.isdir(result))

    @unittest.skipIf('vtk' not in sys.modules, "vtk module required for this feature")
    def test_write_vtk_points(self):

        n = 500
        m = 3

        points = np.random.rand(n, m)

        with tempfile.NamedTemporaryFile() as tmp_file:
            io.write_vtk_points(points, tmp_file.name)

            reader = vtk.vtkPolyDataReader()
            reader.SetFileName(tmp_file.name)
            reader.Update()
            polydata = reader.GetOutput()

            for i in range(0, points.shape[0]):
                with self.subTest(i=i):
                    self.assertTrue(np.allclose(points[i,:], polydata.GetPoint(i)))

    def test_read_elastix_points(self):

        n = 100
        points = np.random.rand(n, 3)

        with tempfile.NamedTemporaryFile(mode='w') as tmp_file:
            tmp_file.write("point\n%d\n" % points.shape[0])
            for i in range(0, points.shape[0]):
                tmp_file.write("%.10f %.10f %.10f\n" % tuple([*points[i,:]]))
            tmp_file.flush()

            result = io.read_elastix_points(tmp_file.name)

            self.assertTrue(np.allclose(points, result))

    def test_write_elastix_points(self):

        n = 100
        points = np.random.rand(n, 3)

        with tempfile.NamedTemporaryFile(mode='w') as tmp_file:

            io.write_elastix_points(points, tmp_file.name)
            result = io.read_elastix_points(tmp_file.name)

            self.assertTrue(np.allclose(points, result))


if __name__ == '__main__':
    unittest.main()
