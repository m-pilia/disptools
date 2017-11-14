import numpy as np
import os
import re
import SimpleITK as sitk
from typing import *
from disptools import *

try:
    import vtk
except ImportError as e:
    print("Warning: cannot import 'vtk' module. " +
          "Some functionalities depending upon it may be unavailable.")


def make_unique_directory(path: str) -> str:
    r""" Create a unique directory.

    If a directory with the given name already exists,
    a suffix in the form `_x` with `x` integer will be
    added at the end of its name.

    Parameters
    ----------
    path : str
        Relative or absolute path for the new directory.

    Returns
    -------
    str
        A string containing the path of the new directory.
    """

    if not os.path.exists(path):
        os.makedirs(path)
        return path

    base_path, name = os.path.split(path)
    name, extension = os.path.splitext(name)

    i = 1
    while True:
        path = os.path.join(base_path, '%s_%d%s' % (name, i, extension))
        if not os.path.exists(path):
            os.makedirs(path)
            return path
        i += 1


def read_rvf(filename: str) -> sitk.Image:
    r""" Read an RVF file.

    Read image data from an RVF file.

    Parameters
    ----------
    filename : str
        Filename for the input vector field.

    Returns
    -------
    result : sitk.Image
        Vector field.
    """

    with open(filename, 'rb') as f:
        size = [int(i) for i in f.readline().split()]
        spacing = [float(i) for i in f.readline().split()]
        data = np.fromfile(f, dtype=np_float_type)

    image = sitk.GetImageFromArray(data.reshape((*size, 3), order='C'))
    image.SetSpacing(tuple(spacing))

    return image


def write_vtk_points(points: np.ndarray, filename: str) -> None:
    r""" Write a set of points to a VTK PolyData file.

    The points are given as a numpy bidimensional array, with a
    row for each point, and three columns, one per component.

    Parameters
    ----------
    points : np.ndarray
        A `n × m` array containing `n` points with `m` components.
    filename : str
        Output file.
    """

    if 'vtk' not in sys.modules:
        raise Exception('write_vtk_points: vtk module is required to use this feature.')

    vtk_points = vtk.vtkPoints()

    for i in range(0, points.shape[0]):
        vtk_points.InsertNextPoint(points[i,0], points[i,1], points[i,2])

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polydata)
    writer.SetFileName(filename)
    writer.Update()


def read_elastix_points(filename: str) -> np.ndarray:
    r""" Read a set of points from a file in Elastix format.

    The points are returned as a two-dimensional numpy array, with a
    row for each point and three columns, one per coordinate.

    Parameters
    ----------
    filename : str
        Output file.

    Returns
    -------
    np.ndarray
        A `n × m` array containing `n` points with `m` components.
    """

    with open(filename, 'r') as f:

        line = f.readline()
        if not re.match(r'\s*point|index\s*', line):
            raise Exception('Invalid point file format')

        line = f.readline()
        try:
            n = int(re.match(r'([0-9]+)', line).group(1))
        except:
            raise Exception('Invalid point file format')

        points = np.empty([n, 3])

        i = 0
        for line in f:
            try:
                m = re.match(r'([0-9.e+-]+)\s+([0-9.e+-]+)\s+([0-9.e+-]+)', line)
                points[i,:] = [float(x) for x in m.groups()]
            except:
                raise Exception('Invalid point file format')

            i += 1

        return points


def write_elastix_points(
        points: np.ndarray,
        filename: str,
        point_format: str = 'point'
        ) -> None:
    r""" Write a set of points to a file in Elastix format.

    The points are passed as a two-dimensional numpy array, with a row
    for each point and three columns, one per coordinate.

    Parameters
    ----------
    points : np.ndarray
        A `n × m` array containing `n` points with `m` components.
    filename : str
        Output file.
    point_format : str
        One of 'point' (default) or 'index'.
    """

    if point_format not in ['point', 'index']:
        raise Exception('Unsupported point format %s' % point_format)

    if len(points.shape) != 2 and points.shape[1] != 3:
        raise Exception('Invalid point array')

    with open(filename, 'w') as f:
        f.write('%s\n' % point_format)
        f.write('%d\n' % points.shape[0])
        for i in range(0, points.shape[0]):
            f.write('%.10e\t%.10e\t%.10e\n' % (points[i,0], points[i,1], points[i,2]))
