import SimpleITK as sitk
import numpy as np
import math
import functools
from typing import *
from disptools import *
import disptools.drawing as drawing

import _disptools

try:
    import itk
except ImportError as e:
    print("Warning: cannot import 'itk' module. " +
          "Some functionalities depending upon it may be unavailable.")


def reject_outliers(data: np.ndarray, k: float = 2.0, absolute=True) -> np.ndarray:
    r""" Reject outliers from an array, based on the median deviation (MD).

    Reject as outliers all values whose deviation from the median is higher
    than :math:`k \cdot MD`.

    Parameters
    ----------
    data : np.ndarray
        Array of values.
    k : float
        Number of MDs for outlier rejection.
    absolute : bool
        If True, use absolute deviations.

    Returns
    -------
    np.ndarray
        An array with the outliers removed.
    """

    dev = data - np.median(data)
    if absolute:
        dev = np.abs(dev)
    md = np.median(dev)
    return data[(dev / md if md else 0.0) < k]


def cubicity(image: sitk.Image) -> float:
    r""" Measure the cubicity of an object.

    The cubicity is defined [4]_ as the ratio between the volume of the
    object and the volume of its bounding cube.
    A cube has cubicity equal to 1, a sphere has cubicity pi/6, and in
    general non-cubic objects have cubicity strictly lesser than 1.

    Here the bounding cube is estimated as the cube whose side is
    equal to the longest side of the oriented bounding box of the
    object.

    References
    ----------
    .. [4] O'Flannery, LJ and O'Mahony, MM, "Precise shape grading of
           coarse aggregate". Magazine of Concrete Research 51.5 (1999),
           pp. 319-324.

    Parameters
    ----------
    image : sitk.Image
        Input binary (sitkUInt8) image.

    Returns
    -------
    float
        A floating point value of cubicity in the interval [0, 1].
    """

    if image.GetPixelID() != sitk.sitkUInt8:
        raise Exception('Unsupported %s pixel type' % image.GetPixelIDTypeAsString())

    # NOTE: the size of the bounding box is already in image space
    # units, while the volume (number of voxels) needs to be multiplied
    # by the voxel size
    dv = functools.reduce(lambda x, y: x * y, image.GetSpacing(), 1.0)

    try:
        lssif = sitk.LabelShapeStatisticsImageFilter()
        lssif.ComputeOrientedBoundingBoxOn()
        lssif.Execute(image)

        (s1, s2, s3) = lssif.GetOrientedBoundingBoxSize(1)
        volume = lssif.GetNumberOfPixels(1) * dv

    except AttributeError:
        # Use ITK as a fallback if the method is not available in
        # SimpleITK

        if 'itk' not in sys.modules:
            raise Exception('sitk_to_itk: itk module is required to use this feature.')

        itk_image = drawing.sitk_to_itk(image)
        li2slmf = itk.LabelImageToShapeLabelMapFilter.IUC3LM3.New(itk_image)
        li2slmf.ComputeOrientedBoundingBoxOn()

        statistics = li2slmf()[0][1]

        #  FIXME GetOrientedBoundingBoxSize() seems to be broken
        #  (s1, s2, s3) = statistics.GetOrientedBoundingBoxSize() #
        (s1, s2, s3) = statistics.GetBoundingBox().GetSize()

        # FIXME GetBoundingBox is in voxels, GetOrientedBoundingBox is
        # in image size instead, so multiply by the spacing when using
        # GetBoundingBox
        (s1, s2, s3) = [x * s for x, s in zip([s1, s2, s3], image.GetSpacing())]

        volume = statistics.GetNumberOfPixels() * dv

    return volume / (max(s1, s2, s3) ** 3)


def sphericity(image: sitk.Image) -> float:
    r""" Measure the sphericity of a object.

    The sphericity is defined [5]_ [6]_ as the ratio between the surface of a sphere
    with the same volume of the object and the surface of the object itself.

    .. math::
        \frac{\pi^{\frac{1}{3}}(6V)^{\frac{2}{3}}}{A}

    A sphere has sphericity equal to 1, non-spherical objects have sphericity
    strictly lesser than 1.

    References
    ----------
    .. [5] Wadell, Hakon. "Volume, shape, and roundness of quartz particles."
           The Journal of Geology 43.3 (1935): 250-280.
    .. [6] Lehmann, Gaëthan. "Label object representation and manipulation with ITK"
           Insight Journal, July-December 2007

    Parameters
    ----------
    image : sitk.Image
        Input binary (sitkUInt8) image.

    Returns
    -------
    float
        A floating point value of sphericity in the interval [0, 1].
    """

    if image.GetPixelID() != sitk.sitkUInt8:
        raise Exception('Unsupported %s pixel type' % image.GetPixelIDTypeAsString())

    lssif = sitk.LabelShapeStatisticsImageFilter()
    lssif.Execute(image)
    return lssif.GetRoundness(1)


def average_jacobian_error(jacobian: sitk.Image, mask: sitk.Image = None) -> float:
    r""" Compute the average error between the input Jacobian and
    a constant unit Jacobian map.

    The average Jacobian error is defined as
    :math:`\frac{1}{|mask|} \sum_{x \in mask} (J(x) - 1)`

    Parameters
    ----------
    jacobian : sitk.Image
        Input Jacobian map.
    mask : sitk.Image
        Optional binary mask, to compute the error only on a ROI.

    Returns
    -------
    float
        The average Jacobian error.
    """

    jac = sitk.GetArrayViewFromImage(jacobian)
    msc = sitk.GetArrayViewFromImage(mask) if mask is not None else np.ones(jac.shape)

    idx = msc > 0.5

    return np.sum(np.abs(jac[idx] - 1.0)) / np.sum(idx)


def minkowski_compactness(image: sitk.Image, norm: Union[int,str] = 1.0) -> float:
    r""" Compute the Minkovski compactness of a binary image.

    Minkowski compactness of an object is defined as the ratio between the volume
    of the object and the volume of a p-ball centred on the object's centre of
    mass, maximised agaist spatial rotations of the ball.

    Parameters
    ----------
    image : sitk.Image
        Input binary image.
    norm : Union[int,str]
        Order of the Minkowski norm ('inf' or 'max' to use the Chebyshev norm).

    Returns
    -------
    float
        Value of the Minkowski norm.
    """

    if image.GetPixelID() != sitk.sitkUInt8:
        raise Exception('Unsupported %s pixel type' % image.GetPixelIDTypeAsString())

    if norm == 'inf' or norm == 'max':
        descriptor = 'cubeness'
    elif norm == 1.0:
        descriptor = 'octahedroness'
    else:
        raise Exception('Unsupported value %s for the norm parameter' % str(norm))

    lssif = sitk.LabelShapeStatisticsImageFilter()
    lssif.Execute(image)

    a = sitk.GetArrayViewFromImage(image) > 0

    return _disptools.shape_descriptor(a, lssif.GetCentroid(1), image.GetSpacing(), descriptor)


def isovolumteric_radius(image: sitk.Image, norm: Union[int, str]) -> float:
    r""" Compute the radius of a ball isovolumetric to the input object

    Parameters
    ----------
    image : sitk.Image
        Input binary image.
    norm : Union[int,str]
        Order of the Minkowski norm ('inf' or 'max' to use the Chebyshev norm).

    Returns
    -------
    float
        Value of the isovolumetric radius.
    """

    volume = np.sum(sitk.GetArrayViewFromImage(image))

    if norm == 2:
        return math.pow(3.0 / (4.0 * np.pi) * volume, 1.0 / 3.0)
    elif norm == 1:
        return math.pow(3.0 / 4.0 * volume, 1.0 / 3.0)
    elif norm == 'inf' or norm == 'max':
        return math.pow(1.0 / 8.0 * volume, 1.0 / 3.0)
    else:
        raise Exception('Unsupported norm %s' % str(norm))


def fitting_index(
        image   : sitk.Image,
        norm    : Union[int, str] = 2.0,
        centre  : Tuple[float, float, float] = None,
        radius  : float = None,
        padding : bool = True
        ) -> float:
    r""" Compute the fitting index of an input object.

    The fitting index of order `p` is defined as the Jaccard coefficient
    computed between the input object and a p-ball centred in the object's
    centre of mass.

    Parameters
    ----------
    image : sitk.Image
        Input binary image of the object.
    norm : Union[int,str]
        Order of the Minkowski norm ('inf' or 'max' to use the Chebyshev norm).
    centre : Tuple[float, float, float]
        Forces the p-ball to be centred in a specific point.
    radius : float
        Force the radius of the p-ball.
    padding : bool
        If `True`, add enough padding to be sure that the ball will entirely
        fit within the volume.

    Returns
    -------
    float
        Value of the fitting index.
    """

    if image.GetPixelID() != sitk.sitkUInt8:
        raise Exception('Unsupported %s pixel type' % image.GetPixelIDTypeAsString())

    if centre is None:
        # Use the centroid as centre
        lssif = sitk.LabelShapeStatisticsImageFilter()
        lssif.Execute(image)
        centre = lssif.GetCentroid(1)

    if padding:
        # Add some padding to be sure that an isovolumetric 1-ball can fit
        # within the same volume of a sphere touching the boundary
        pad = tuple([x // 4 for x in image.GetSize()])
        image = sitk.ConstantPad(image, pad, pad, 0)
        image.SetOrigin((0,0,0))
        centre = tuple([x + y for x, y in zip(centre, pad)])

    if radius is None:
        radius = isovolumteric_radius(image, norm)

    size = image.GetSize()

    sphere = drawing.create_sphere(radius, size=size, centre=centre, norm=norm) > 0

    data = sitk.GetArrayViewFromImage(image) == 1.0
    reference = sitk.GetArrayViewFromImage(sphere) == 1.0

    union = np.logical_or(data, reference)
    intersection = np.logical_and(data, reference)

    return np.sum(intersection) / np.sum(union)
