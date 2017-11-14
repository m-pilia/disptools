import SimpleITK as sitk
import numpy as np
from disptools import *

def create_target_volume(image: sitk.Image, atrophy_rate: float):
    r""" Create a target volume map for the PREDICT tool.

    Given an input segmentation map, create mask target volume map for the
    PREDICT tool, with the given atrophy rate. The volume change in
    :math:`(x,y,z)` is defined over a cubic pach formed by the voxel
    :math:`(x,y,z)`, its three successors along the axis and their respective
    successors to close the cube. Hence, if the input mask has size
    :math:`m \times n \times s`, the volume map has size
    :math:`(m-1) \times (n-1) \times (s-1)`.

    .. note::
        The output must be in np.float32, since the PREDICT tool uses
        the C `float` type as hard-coded type.

    Parameters
    ----------
    image : sitk.Image
        Input segmentation map.
    atrophy_rate : float
        Target atropy rate for the segmented ROI.

    Returns
    -------
    sitk.Image
        A SimpleITK image object with the target volume map.
    """

    (x,y,z) = tuple([x - 1 for x in image.GetSize()])

    volumes = 1.0 - atrophy_rate * sitk.GetArrayViewFromImage(I)[0:x,0:y,0:z]

    result = sitk.GetImageFromArray(volumes)
    result.CopyInformation(image)

    return result


def read_deformation(filename: str, size: Tuple[int, int, int]) -> sitk.Image:
    r""" Read a deformation field from a file in the PREDICT format

    The deformation field is stored as binary uncompressed 32-float data.

    PREDICT uses a different convention for the coordinates (with rows
    along the :math:`y` axis, columns along the :math:`x` axis, and slices
    along the :math:`z` axis), so a permutation of the components in the
    result is required. The identity is subtracted in order to convert the
    deformation to a displacement field.

    .. note::
        The size of the volume is not stored within PREDICT files, so it needs
        to be known and passed as an argument.

    Parameters
    ----------
    filename : str
        Input file.
    size : (int, int, int)
        Size of the volume `(x,y,z)`.

    Returns
    -------
    sitk.Image
        A SimpleITK image object containing the corresponding displacement field.
    """

    with open(filename, 'rb') as f:
        a = np.fromfile(f, dtype=np.float32)

    a = a.reshape((*size, 3))
    b = np.empty(a.shape)

    # Convert from PREDICT's coordinate system to ITK's.
    # Also convert from deformation field to displacement field, i.e.
    # subtract the identity transformation.
    for x, y, z in np.ndindex(b.shape[0:3]):
        b[x,y,z,0] = a[x,y,z,1] - (z+1)
        b[x,y,z,1] = a[x,y,z,0] - (y+1)
        b[x,y,z,2] = a[x,y,z,2] - (x+1)

    return sitk.GetImageFromArray(b)


def read_img(filename: str, size: Tuple[int, int, int]) -> sitk.Image:
    r""" Read an image from file in the PREDICT format.

    PREDICT uses a different convention for the coordinates (with rows
    along the :math:`y` axis, columns along the :math:`x` axis, and slices
    along the :math:`z` axis), so a permutation of the components in the
    result is required.

    .. note::
        The size of the volume is not stored within PREDICT files, so it needs
        to be known and passed as an argument.

    Parameters
    ----------
    filename : str
        Input file.
    size : (int, int, int)
        Size of the volume `(x,y,z)`.

    Returns
    -------
    sitk.Image
        a SimpleITK image object containing the corresponding image
    """

    with open(filename, 'rb') as f:
        a = np.fromfile(f, dtype=np.float32)

    a = a.reshape(size)
    b = np.empty(a.shape)

    # Convert the coordinates of the points from the tool's
    # coordinate system to ITK's
    for x, y in np.ndindex(a.shape[0:2]):
        b[x,y,:] = a[y,x,:]

    return sitk.GetImageFromArray(b)


def write_img(image: sitk.Image, filename: str) -> None:
    r""" Write an image to file in the PREDICT format.

    PREDICT uses a different convention for the coordinates (with rows
    along the :math:`y` axis, columns along the :math:`x` axis, and slices
    along the :math:`z` axis), so a permutation of the components in the
    result is required.

    .. note::
        The size of the volume will not be stored within the file.

    Parameters
    ----------
    filename : str
        Input file.
    size : (int, int, int)
        Size of the volume `(x,y,z)`.

    Returns
    -------
    sitk.Image
        A SimpleITK image object containing the corresponding image.


    """

    a = sitk.GetArrayViewFromImage(image)
    b = np.empty(a.shape)

    # Convert the coordinates of the points
    for x, y in np.ndindex(b.shape[0:2]):
        b[x,y,:] = a[y,x,:]

    with open(filename, 'wb') as f:
        b.tofile(f)


