import sys
import SimpleITK as sitk
import numpy as np
from scipy.spatial import distance
import math
from typing import *
from disptools import *

try:
    import itk
except ImportError as e:
    print("Warning: cannot import 'itk' module. " +
          "Some functionalities depending upon it may be unavailable.")


def resize_image(
        image        : sitk.Image,
        new_size     : Tuple[int, ...],
        interpolator : int = sitk.sitkLinear
        ) -> sitk.Image:
    """ Resample an image in a grid of given size.

    Parameters
    ----------
    image : sitk.Image
        An input image.
    new_size : Tuple[int, ...]
        A tuple of integers expressing the new size.
    interpolator : int
        A SimpleITK interpolator enum value.

    Returns
    -------
    sitk.Image
        Resized image.
    """

    if type(image) is not sitk.SimpleITK.Image:
        raise Exception("unsupported image object type")

    if type(new_size) != type((1,1)):
        raise Exception("new_size must be a tuple of integers")

    size = image.GetSize()
    if len(new_size) != len(size):
        raise Exception("new_size must match the image dimensionality")

    spacing = []
    for s, x, nx in zip(image.GetSpacing(), size, new_size):
        spacing.append(s * x/nx)

    resampler = sitk.ResampleImageFilter()
    resampler.SetSize(new_size)
    resampler.SetOutputSpacing(tuple(spacing))
    resampler.SetInterpolator(interpolator)

    return resampler.Execute(image)


def polar_conversion(
        norm   : Union[int, str],
        centre : Tuple[float, float, float],
        radius : float,
        x      : float,
        y      : float,
        z      : float
        ) -> Tuple[int, int, int]:
    """ Convert Cartesian to normalised spherical coordinates.

    Parameters
    ----------
    norm : Union[int, str]
        Minkovski norm ('max' or 'inf' for the Chebyshev norm).
    centre : Tuple[float, float, float]
        Coordinates of the centre.
    radius : float
        Radius of the ball.
    x : float
        Cartesian coordinate.
    y : float
        Cartesian coordinate.
    z : float
        Cartesian coordinate.

    Returns
    -------
    Tuple[int, int, int]
        A tuple (rho, theta, phi) of spherical coordinates with rho
        normalised w.r.t. the radius of the sphere.
    """

    metric = 'chebyshev' if norm in ['max', 'inf'] else 'minkowski'
    args = {'p': float(norm)} if norm not in ['max', 'inf'] else {}
    centre = np.array(centre).reshape((1, 3))

    shape = x.shape
    data = np.column_stack((x.flatten(), y.flatten(), z.flatten()))
    rho = distance.cdist(data, centre, metric, **args).reshape(shape)

    xn = x - centre[0,0]
    yn = y - centre[0,1]
    zn = z - centre[0,2]

    theta = np.where(rho > 0.0, np.arccos(np.divide(zn, rho, out=np.zeros(rho.shape), where=(rho > 0.0))), 0.0)
    phi = np.arctan2(yn, xn)
    rho = rho / radius # normalise rho to [0,1]

    return rho, theta, phi


def create_sphere(
        radius         : float,
        size           : Union[int, Tuple[int, int, int]] = None,
        centre         : Tuple[float, float, float] = None,
        fg_val         : float = 1.0,
        bg_val         : float = 0.0,
        norm           : Union[int, str] = 2,
        sigma          : float = 0.0,
        neighbourhood  : int = 3,
        value_function : Callable[[float, float, float], float] = None,
        points_file    : str = ''
        ) -> sitk.Image:
    """ Create a volume image containing a spherical object.

    Parameters
    ----------
    radius : float
        Radius of the ball.
    size : Union[int, Tuple[int, int, int]]
        Size of the volume image.
    centre : Tuple[float, float, float]
        Centre of the sphere.
    fg_val : float
        Value in the foreground, i.e. inside the ball;
        for variable values, see value_function.
    bg_val : float
        Value in the background, i.e. outside the ball.
    norm : Union[int, str]
        Minkovski norm used to define the metric ('inf'
        or 'max' to use the Chebyshev norm).
    sigma : float
        Sigma parameter for a Gaussian smoothing of the result.
    neighbourhood : int
        Neighbourhood size for a Gaussian smoothing of the result.
    value_function : Callable[[float, float, float], float]
        Function that associate an intensity to each point within
        the ball; it takes three arguments (normalised spherical
        coordinates of the point) and returns a scalar intensity;
        the input is in vector form, and operations within this
        function should be vectorised (e.g. using numpy functions).
    points_file : str
        Name of the file used to store a set of reference
        points associated to the ball (empty string means to not
        store the reference points).

    Returns
    -------
    sitk.Image
        A volume image containing a ball.
    """

    if size is None:
        size = 2 * radius + 3

    if type(size) is not tuple:
        shape = (size, size, size)
    else:
        shape = tuple(reversed(size))

    if centre is None:
        centre = tuple([x // 2 for x in shape])

    def function(z, y, x):
        rho, theta, phi = polar_conversion(norm, centre, radius, x, y, z)
        if value_function is None:
            return np.where(rho <= 1.0, fg_val, bg_val)
        else:
            return np.where(rho <= 1.0, value_function(rho, theta, phi), bg_val)

    a = np.fromfunction(function, shape).astype(np_float_type)

    sphere = sitk.GetImageFromArray(a)

    if sigma > 0.0:
        sphere = sitk.DiscreteGaussian(sphere, sigma, neighbourhood)

    if points_file != '':

        def within(point: Tuple[int, int, int]) -> bool:
            """ Check whether a point is within the volume.
            """

            point_x = point[0] >= 0.0 and point[0] < shape[2]
            point_y = point[1] >= 0.0 and point[1] < shape[1]
            point_z = point[2] >= 0.0 and point[2] < shape[0]
            return point_x and point_y and point_z

        with open(points_file, 'w') as f:
            c = np.array(centre)
            x = np.array([1,0,0])
            y = np.array([0,1,0])
            z = np.array([0,0,1])
            r = radius
            points = [
                c,
                c + r*x, c + 0.5 * r*x,
                c - r*x, c - 0.5 * r*x,
                c + r*y, c + 0.5 * r*y,
                c - r*y, c - 0.5 * r*y,
                c + r*z, c + 0.5 * r*z,
                c - r*z, c - 0.5 * r*z,
                c + 0.5 * r*x + 0.5 * r*y,
                c - 0.5 * r*x + 0.5 * r*y,
                c + 0.5 * r*x - 0.5 * r*y,
                c - 0.5 * r*x - 0.5 * r*y,
                c + 0.5 * r*x + 0.5 * r*z,
                c - 0.5 * r*x + 0.5 * r*z,
                c + 0.5 * r*x - 0.5 * r*z,
                c - 0.5 * r*x - 0.5 * r*z,
                c + 0.5 * r*y + 0.5 * r*z,
                c - 0.5 * r*y + 0.5 * r*z,
                c + 0.5 * r*y - 0.5 * r*z,
                c - 0.5 * r*y - 0.5 * r*z,
            ]
            points = [p for p in points if within(p)]
            f.write('point\n%d\n' % len(points))
            for x, y, z in points:
                f.write('%f\t%f\t%f\n' % (x, y, z))

    return sphere


def mask(
        image    : sitk.Image,
        mask     : sitk.Image,
        jacobian : bool = False
        ) -> sitk.Image:
    """ Mask an image (special meaning for Jacobian maps).

    Parameters
    ----------
    image : sitk.Image
        Input image to be masked (possibly float).
    mask : sitk.Image
        Mask (possibly float).
    jacobian : bool
        If true, the background after masking is set to
        one, if false to zero.

    Returns
    -------
    sitk.Image
        The masked image.
    """

    if jacobian:
        background = np.logical_not(sitk.GetArrayViewFromImage(mask)).astype(np_float_type)
        background = sitk.GetImageFromArray(background)
        background = sitk.Cast(background, image.GetPixelID())
        background.CopyInformation(image)

        cast_mask = sitk.Cast(mask, image.GetPixelID())
        cast_mask.CopyInformation(image)

        result = sitk.Multiply(image, cast_mask)
        result = sitk.Add(result, background)

    else:
        cast_mask = sitk.Cast(mask, image.GetPixelID())
        cast_mask.CopyInformation(image)
        result = sitk.Multiply(image, cast_mask)

    return result


def float_dilate(image: sitk.Image, dilate: int) -> sitk.Image:
    """ Dilate a float valued image.

    Parameters
    ----------
    image : sitk.Image
        Image to be dilated.
    dilate : int
        Radius of the structuring element.

    Returns
    -------
    sitk.Image
        Dilated image.
    """

    if type(dilate) is not int or dilate < 1:
        raise Exception('float_dilate: invalid dilation value')

    original_type = image.GetPixelID()
    image = sitk.Cast(image, sitk.sitkUInt8)
    image = sitk.BinaryDilate(image, dilate)
    image = sitk.Cast(image, original_type)

    return image


def sitk_to_itk(image: sitk.Image) -> Any:
    """ Function to convert an image object from SimpleITK to ITK.

    .. note::
        Data is copied to the new object (deep copy).

    Parameters
    ----------
    image : sitk.Image
        Input image.

    Returns
    -------
    any
        Image in ITK format.
    """

    if 'itk' not in sys.modules:
        raise Exception('sitk_to_itk: itk module is required to use this feature.')

    a = sitk.GetArrayViewFromImage(image)

    if len(a.shape) < 4:
        result = itk.GetImageFromArray(a)

    else:
        # NOTE: This workaround is implemented this way since it
        # seems that itk.GetImageFromArray() is not working properly
        # with vector images.

        region = itk.ImageRegion[3]()
        region.SetSize(image.GetSize())
        region.SetIndex((0, 0, 0))

        PixelType = itk.Vector[itk_float_type, 3]
        ImageType = itk.Image[PixelType, 3]

        # Create new image
        result = ImageType.New()
        result.SetRegions(region)
        result.Allocate()
        result.SetSpacing(image.GetSpacing())

        # Copy image data
        b = itk.GetArrayViewFromImage(result)
        b[:] = a[:]

    result.SetSpacing(image.GetSpacing())
    result.SetOrigin(image.GetOrigin())

    return result


def sin_vector_field(n: int) -> sitk.Image:
    """ Generate an example of continuous vector field.

    Parameters
    ----------
    n : int
        Size of the volume image.

    Returns
    -------
    sitk.Image
        A vector-valued volume image containing a continuous
        vector field.
    """

    x = np.linspace(0, 0.5, n)
    xx, yy, zz = np.meshgrid(x, x, x)
    data = np.empty((n, n, n, 3), dtype=np_float_type)
    data[:,:,:,0] = np.sin(xx**2 + yy**2 + zz**2)
    data[:,:,:,1] = data[:,:,:,0]
    data[:,:,:,2] = data[:,:,:,0]

    return sitk.GetImageFromArray(data)


def extract_slice(
        image    : sitk.Image,
        index    : int = None,
        axis     : int = 2,
        rescale  : bool = False,
        colormap : int = sitk.ScalarToRGBColormapImageFilter.Grey,
        window   : Tuple[float, float] = None
        ) -> sitk.Image:
    """ Extract a slice.

    This function takes care of resampling the image on a uniform grid with
    unit spacing, that can be exported to format that do not support
    anisotropic voxel spacing.

    Parameters
    ----------
    image : sitk.Image
        Input image.
    index : int
        Integer index of the slice to extract.
    axis : int
        Integer index of the axis orthogonal to the slice plane.
    rescale : bool
        If True, do min-max rescaling of the image intensity to [0,255].
    colormap : int
        One of the colormaps defined in sitk.ScalarToRGBColormapImageFilter.
        If None, the result will be a grayscale image; if not None, map
        the intensity values to RGB colours through the selected colormap.
    window : (float, float)
        Tuple composed by a couple of minimum and maximum values for
        intensity windowing.

    Returns
    -------
    sitk.Image
        A sitk.sitkUInt8 image containing the desired slice.
    """

    if window is not None:
        image = sitk.IntensityWindowing(image, *window)

    if rescale:
        image = sitk.RescaleIntensity(image)

    image = sitk.Cast(image, sitk.sitkUInt8)

    size = list(image.GetSize())

    if index is None:
        index = size[axis] // 2

    size[axis] = 0
    size = tuple(size)

    indices = tuple([index if i == axis else 0 for i in range(0, len(size))])

    slice_image = sitk.Extract(image, size=size, index=indices)

    new_size = [int(x * dx) for x, dx in zip(slice_image.GetSize(), slice_image.GetSpacing())]
    slice_image = sitk.Resample(slice_image, tuple(new_size))
    slice_image.SetSpacing(tuple([1 for _ in new_size]))

    if colormap is not None:
        slice_image = sitk.ScalarToRGBColormap(slice_image, colormap, useInputImageExtremaForScaling=False)

    return slice_image
