import math
import os
import re
import sys
import threading
import warnings
from functools import reduce
from glob import glob
from typing import *

import numpy as np
import scipy.interpolate
import SimpleITK as sitk

from disptools import *
import disptools.drawing as drawing
import disptools.io as io

import _disptools


def regularise(jacobian: sitk.Image, epsilon: float = 1e-5) -> sitk.Image:
    r""" Regularise the Jacobian, removing singularities.

    Given a 3D scalar image, replace all the entries that are smaller
    than `epsilon` with `epsilon`.

    Parameters
    ----------
    jacobian : sitk.Image
        Input Jacobian map
    epsilon  : float
        Lower threshold for the Jacobian.

    Returns
    -------
    sitk.Image
        Thresholded Jacobian.
    """

    jacobian = sitk.Cast(jacobian, sitk_float_type)

    if (3 != len(jacobian.GetSize())):
        raise Exception("Wrong jacobian dimensionality")

    # Object for the result
    result = sitk.Image(jacobian)

    # Call function from the underlying C library
    _disptools.regularise(sitk.GetArrayViewFromImage(result), epsilon)

    return result


def jacobian(field: sitk.Image) -> sitk.Image:
    r""" Compute the Jacobian map of a vector field.

    Parameters
    ----------
    field : sitk.Image
        Input vector field.

    Returns
    -------
    sitk.Image
        The Jacobian of the given vector field.
    """

    field = sitk.Cast(field, sitk_vector_float_type)

    image = sitk.GetArrayViewFromImage(field)

    if 4 != len(image.shape) or 3 != image.shape[3]:
        raise Exception("Wrong input dimensionality")

    # Convert to the library's memory layout
    shape = image.shape[0:3]
    image_data = np.empty((3, *shape), dtype=np_float_type)
    image_data[0,:,:,:] = image[:,:,:,0]
    image_data[1,:,:,:] = image[:,:,:,1]
    image_data[2,:,:,:] = image[:,:,:,2]

    # Object for the result
    result = np.zeros(shape, dtype=np_float_type)

    # Call function from the underlying C library
    _disptools.jacobian(field.GetSpacing(), image_data, result)

    result = sitk.GetImageFromArray(result)
    result.CopyInformation(field)

    return result


# This function wraps the C call.
def _displacement(
        image         : sitk.Image,
        mask          : sitk.Image = None,
        initial_guess : sitk.Image = None,
        epsilon       : float = 9.99e-4,
        tolerance     : float = 0.2,
        it_max        : int = 50000,
        alpha         : float = 1.2,
        beta          : float = 0.5,
        gamma         : float = .1,
        delta         : float = 1e-3,
        zeta          : float = 100.0,
        theta         : float = 1e-6,
        iota          : float = 1e-9,
        strict        : bool = False,
        eta           : float = 0.1,
        eta_max       : float = 0.4,
        algorithm     : str = 'gradient',
        gpu_id        : int = -1,
        ):
    r""" Compute a displacement field that realises a prescribed Jacobian.

    .. note::
        This function should not be called directly. Please use its wrapper `displacement'.

    Parameters
    ---------
    See the documentation for `displacement`.

    Returns
    -------
    sitk.Image
        A displacement field whose Jacobian matches the input.
    """

    jacobian = sitk.GetArrayViewFromImage(image).astype(np_float_type, copy=True)

    if mask == None:
        # Default mask: whole image
        mask = np.ones(jacobian.shape, dtype=np.bool)
    else:
        # Get mask as numpy array of booleans
        mask = sitk.GetArrayViewFromImage(mask).astype(bool, copy=True)
        if mask.shape != jacobian.shape:
            raise Exception("The shapes of the Jacobian and the mask must agree")

    # Create objects for the result, initialise initial guess
    if initial_guess is not None:
        data = sitk.GetArrayViewFromImage(initial_guess)
        if data.shape[0:3] != jacobian.shape:
            raise Exception("The shapes of the Jacobian and the initial guess must agree")
        field_tmp = np.empty((3, *jacobian.shape), dtype=np_float_type)
        field_tmp[0,:,:,:] = data[:,:,:,0]
        field_tmp[1,:,:,:] = data[:,:,:,1]
        field_tmp[2,:,:,:] = data[:,:,:,2]
    else:
        field_tmp = np.zeros((3, *jacobian.shape), dtype=np_float_type)

    # Arguments
    args = [
        image.GetSpacing(),
        jacobian,
        mask,
        epsilon,
        tolerance,
        eta,
        eta_max,
        alpha,
        beta,
        gamma,
        delta,
        zeta,
        theta,
        iota,
        strict,
        it_max,
        field_tmp,
        algorithm,
        gpu_id,
    ]

    # Call the function in a separate thread
    #
    # The thread cannot actually be killed from the REPL while busy, but
    # at least this way the REPL can handle SIGINT. Using a process
    # instead of a thread would allow to kill the call, but it does not
    # work in the REPL, so there is no point.
    #
    # The real solution would be to put a signal handler inside the C
    # routines. There are resources to be freed before interrupting, and
    # they are shared betweeen OpenMP threads, so it requires to use a
    # termination flag inside the OpenMP sections.
    t = threading.Thread(target=_disptools.displacement, args=args, daemon=True)
    t.start()
    t.join()

    # Convert to ITK's image memory layout
    field = np.empty((*jacobian.shape, 3), dtype=np_float_type)
    field[:,:,:,0] = field_tmp[0,:,:,:]
    field[:,:,:,1] = field_tmp[1,:,:,:]
    field[:,:,:,2] = field_tmp[2,:,:,:]

    # Convert the result to match the input type
    field = sitk.GetImageFromArray(field)
    field.CopyInformation(image)

    return field


def redistribute_volume_change(image: sitk.Image, mask: sitk.Image) -> sitk.Image:
    r""" Redistribute the volume change over the image.

    Redistribute the change of volume within the body on the background,
    and enforce the total volume change over the entire image to be zero.

    Parameters
    ----------
    image : sitk.Image
        Input Jacobian.
    mask : sitk.Image
        Binary mask of the body volume.

    Returns
    -------
    sitk.Image
        A new Jacobian map with redistributed volume changes.
    """

    data = sitk.GetArrayFromImage(image)
    index = sitk.GetArrayViewFromImage(mask)

    correction = -(np.sum(data) - data.size) / (data.size - np.count_nonzero(index))
    data[index == 0.0] += correction

    result = sitk.GetImageFromArray(data)
    result.CopyInformation(image)

    return result


# This function is a wrapper of `_displacement()` that adds padding/cropping
# and handles the multi-resolution pyramid.
def displacement(
        jacobian     : sitk.Image,
        *,
        levels       : int = 1,
        pad          : int = 0,
        redistribute : bool = False,
        mask          : sitk.Image = None,
        initial_guess : sitk.Image = None,
        epsilon       : float = 9.99e-4,
        tolerance     : float = 0.2,
        it_max        : Union[int, List[int]] = 50000,
        alpha         : float = 1.2,
        beta          : float = 0.5,
        gamma         : float = 0.1,
        delta         : float = 1e-3,
        zeta          : float = 10.0,
        theta         : float = 1e-6,
        iota          : float = 1e-9,
        strict        : bool = False,
        eta           : float = 0.1,
        eta_max       : float = 0.4,
        algorithm     : str = 'gradient',
        gpu_id        : int = -1
        ) -> sitk.Image:
    r""" Generate a displacement field that realises a given Jacobian.

    Given a 3D scalar image encoding a Jacobian map, compute a
    3D vector image encoding a vector field whose Jacobian map
    matches the input up to a certain tolerance.
    The three algorithms provided are:

    * ``gradient``: a gradient descent method (default).
    * ``greedy``: a greedy search method based on the method proposed in [1]_.
    * ``matching``: a volume matching routine based on gradient descent,
      published in [2]_ and [3]_. The implementation comes from
      the `atrophysim tool`_.

    The initial value of the step length in the gradient descent is
    given by the parameter ``eta``. The ``gradient`` and ``matching``
    algorithms use an `Armijo condition`_ to control the step length,
    in the form

    .. math::
        E(d - \eta \nabla E(d)) - E(d) \le -\gamma \eta ||\nabla E(d)||^2

    where :math:`d` is the displacement, :math:`E` the loss function,
    :math:`\eta` the current step length, and :math:`\gamma \in (0, 1)` a
    parameter of the condition. At each iteration the step length is increased
    multiplying it by ``alpha``, and if the Armijo condition is not met after
    the update, ``eta`` is decreased multiplying it by ``beta`` until the
    truth of the inequality in the Armijo condition is restored. A maximum
    value for ``eta`` can be set through the parameter ``eta_max``.

    The ``gradient`` and ``matching`` algorithms have a regularisation term
    that penalises values of the Jacobian below a certain threshold, given
    by ``delta``. The importance of the regularisation term is controlled
    by the parameter ``zeta`` (set to ``0`` to have no regularisation).

    Termination is controlled by a condition on the improvement on the result
    and one on the step length. If the percentual improvement of the cost
    drops below the value given by ``theta``, the algorithm terminates.
    Termination happens also if the step length becomes smaller than
    ``iota``.

    .. _atrophysim tool: https://www.nitrc.org/projects/atrophysim
    .. _Armijo condition: https://en.wikipedia.org/wiki/Wolfe_conditions

    .. note::
        The displacement is generally not accurate on image boundary voxels.

    .. note::
        The C verbose output is written to `stdout`. If you want to capture
        it from within Python, the `wurlitzer package`_ might be helpful.

    .. warning::
        This function calls a C routine which cannot be interrupted from
        the REPL.

    .. _wurlitzer package: https://github.com/minrk/wurlitzer

    References
    ----------
    .. [1] van Eede, M. C., Scholz, J., Chakravarty, M. M., Henkelman, R. M., and Lerch, J. P.
           "Mapping registration sensitivity in MR mouse brain images." Neuroimage 82 (2013), 226–236.
    .. [2] Karaçali, B., and Davatzikos, C. "Estimating topology preserving and smooth displacement fields."
           IEEE Transactions on Medical Imaging 23, 7 (2004), 868–880.
    .. [3] Karaçali, B., and Davatzikos, C. "Simulation of tissue atrophy using a topology preserving
           transformation model." IEEE transactions on medical imaging 25, 5 (2006), 649–652.

    Parameters
    ----------
    jacobian : sitk.Image
        Input Jacobian.
    levels : int
        Number of levels in the multi-resolution pyramid; the size of
        the image along each direction is halven at each level.
    pad : int
        Thickness of the zero-padding around the volume (0 for
        the mask, 1.0 for the Jacobian) to be used during the
        computation. The padding is removed before returning the result.
    redistribute : bool
        Redistribute the volume change inside the mask to the background.
    mask : sitk.Image
        Binary mask for the region of interest.
    initial_guess : sitk.Image
        Initial estimation of the solution. The default is a null
        displacement field.
    epsilon : float
        A floating point value, representing the tolerance per
        voxel on the Jacobian of the resulting field.
    tolerance : float
        Tolerance on Jacobian outside the mask.
    it_max : Union[int, List[int]]
        Maximum number of iterations allowed. If it is a list, its
        length must match the number of levels in the multi-resolution
        pyramid, and each value is used for a single level, with the
        first element of the list representing the level with lowest
        resolution. If it is a scalar, then the same number of
        iterations is used for all pyramid levels.
    alpha : float
        Coefficient that controls the increase of the step length.
    beta : float
        Coefficient that controls the decrease of the step length.
    gamma : float
        Armijo-Goldstein parameter.
    delta : float
        Lower threshold for Jacobian regularisation.
    zeta : float
        Weight for the regularisation term.
    theta : float
        Terminate if the percentage improvement of the cost per
        iteration drops below this value.
    iota : float
        Terminate if the step length drops below this value.
    strict : bool
        If True, reject iterations that not decrease the maximum
        voxel error.
    eta : float
        Initial step length.
    eta_max : float
        Maximum step length allowed.
    algorithm : str
        Algorithm to generate the field, one of `greedy`, `gradient`, or `matching`.
    gpu_id : int
        Id of the CUDA device used to run the GPU implementation. If
        equal to `-1`, the CPU implementation is used instead. Requires
        a build of disptools with CUDA support enabled.

    Returns
    -------
    sitk.Image
        A displacement field whose Jacobian matches the input.
    """

    # Filter parameters used locally and parameters propagated to the
    # wrapped function
    parameters = locals().copy()
    used = ['jacobian', 'levels', 'pad', 'redistribute', 'mask',
            'initial_guess', 'it_max']
    [parameters.pop(p) for p in used]

    jacobian = sitk.Cast(jacobian, sitk_float_type)

    if mask is None:
        mask = np.ones(tuple(reversed(jacobian.GetSize())), dtype=np_float_type)
        mask = sitk.GetImageFromArray(mask)
        mask.CopyInformation(jacobian)
    else:
        mask = sitk.Cast(mask, sitk_float_type)

    size = jacobian.GetSize()
    origin = jacobian.GetOrigin()
    spacing = jacobian.GetSpacing()
    direction = jacobian.GetDirection()

    # Ensure consistency of the coordinate system
    def make_consistent(img, name, interpolator):
        if img is None:
            return img

        inconsitent = []
        if img.GetSize() != size:
            inconsitent.append('size')
        if img.GetOrigin() != origin:
            inconsitent.append('origin')
        if img.GetSpacing() != spacing:
            inconsitent.append('spacing')
        if img.GetDirection() != direction:
            inconsitent.append('direction')

        if inconsitent != []:
            inconsitent = ' and '.join(inconsitent)
            warnings.warn("%s of '%s' " % (inconsitent, name) +
                          "inconsistent with the Jacobian, " +
                          "resampling to a common coordinate space")
            if interpolator != sitk.sitkNearestNeighbor:
                img = sitk.SmoothingRecursiveGaussian(img, 2.0)
            return sitk.Resample(img, jacobian, sitk.Transform(), interpolator)
        else:
            return img
    mask = make_consistent(mask, 'mask', sitk.sitkNearestNeighbor)
    initial_guess = make_consistent(initial_guess, 'initial_guess', sitk.sitkLinear)

    # Add a voxel of zero-flux padding anyway since the algorithm
    # will not compute the displacement field on boundary voxels
    if pad > 0:
        pad += 1
    else:
        pad = 1
    pad = ((pad, pad, pad), (pad, pad, pad))

    if redistribute:
        jacobian = redistribute_volume_change(jacobian, mask)

    mask = sitk.ConstantPad(mask, *pad, 0)
    jacobian = sitk.ZeroFluxNeumannPad(jacobian, *pad)

    # Create image pyramid
    jacobian_pyramid = [jacobian]
    mask_pyramid = [mask]
    for i in range(1, levels):
        new_size = tuple(map(lambda x: x // 2, jacobian_pyramid[i-1].GetSize()))
        jacobian_pyramid.append(drawing.scale_image(jacobian, new_size))
        mask_pyramid.append(drawing.scale_image(mask, new_size, sitk.sitkNearestNeighbor))

    # Set maximum number of iterations for each pyramid level
    if not isinstance(it_max, list):
        it_max = [it_max for i in range(levels)]
    elif len(it_max) != levels:
        raise ValueError('len(it_max) should equal the value of `levels`')

    # Set initial guess
    field = initial_guess

    # Multi-resolution algorithm
    for jacobian, mask, it in zip(jacobian_pyramid[::-1], mask_pyramid[::-1], it_max):
        size = jacobian.GetSize()
        logging.info('Size %s' % str(size))
        field = drawing.scale_image(field, size) if field is not None else None
        field = _displacement(jacobian, mask, initial_guess=field, it_max=it, **parameters)

    # Remove padding from the result
    field = sitk.Crop(field, *pad)
    field.SetOrigin(origin)
    field.SetSpacing(spacing)
    field.SetDirection(direction)

    return field


def average_jacobian_from_displacements(input_filenames_pattern: str, epsilon: float = 1e-5) -> sitk.Image:
    r""" Compute the average Jacobian of a set of displacement fields.

    This function reads a collection of displacement fields from files
    (``rvf`` or any other format readable by SimpleITK) and computes
    the average Jacobian of the deformation associated to them
    (defined as the geometric mean computed under logarithm).
    It accepts a string argument containing a `glob pattern`_ to the
    input displacement files, and a second optional argument setting
    a minimum threshold for the Jacobian.

    For instance, assuming there is a folder ``/home/user/my_displacements``
    containing a set of displacement fields in ``vtk`` format, the average
    Jacobian can be computed with

    >>> average_jacobian = disptools.displacements.average_jacobian_from_displacements('/home/user/my_jacobians/*.vtk')

    The average Jacobian is defined as the geometric mean computed under
    logarithm.

    .. _glob pattern: https://en.wikipedia.org/wiki/Glob_(programming)

    Parameters
    ----------
    input_filenames_pattern : str
        A glob pattern for the  displacement files in RVF or another
        format readable by SimpleITK.
    epsilon : float
        Minimum threshold for the Jacobian, all values below `epsilon`
        will be replaced with `epsilon`.

    Returns
    -------
    sitk.Image
        The average Jacobian of the given displacements.
    """

    total_jacobian = None

    filenames = glob('%s' % input_filenames_pattern)
    n = len(filenames)
    logging.debug('Files to process: %d' % n)

    logging.debug('Starting')
    i = 1
    for f in filenames:
        logging.debug('Processing file %3d/%d %s' % (i, n, f))

        try:
            file_format = re.match(r'.*\.([^.]+)$', f).group(1)
        except:
            logging.debug('Skipping file %s' %i)
            continue

        if file_format == 'rvf':
            I = io.read_rvf(f)
        else:
            I = sitk.ReadImage(f)

        J = jacobian(I)
        J = regularise(J, epsilon)

        if total_jacobian is None:
            total_jacobian = np.zeros(tuple(reversed(I.GetSize())), dtype=np_float_type)

        total_jacobian += np.log(sitk.GetArrayViewFromImage(J))

        i += 1

    average_jacobian = np.exp(1.0/n * total_jacobian)
    output_jacobian = sitk.GetImageFromArray(average_jacobian)
    output_jacobian.CopyInformation(I)

    return output_jacobian


def average_jacobian(input_filenames_pattern: str, epsilon: float = 1e-5) -> sitk.Image:
    r""" Compute the average of a set of Jacobians.

    This function reads a collection of Jacobian maps from files (any format
    readable by SimpleITK) and computes their average Jacobian (defined as the
    geometric mean computed under logarithm). It accepts a string argument
    containing a `glob pattern`_ to the input files, and a second optional
    argument setting a minimum threshold for the Jacobian.

    For instance, assuming there is a folder ``/home/user/my_jacobians``
    containing a set of Jacobian maps in ``vtk`` format, the average can
    be computed with

    >>> average_jacobian = disptools.displacements.average_jacobian('/home/user/my_jacobians/*.vtk')

    The average Jacobian is defined as the geometric mean computed under
    logarithm.

    .. _glob pattern: https://en.wikipedia.org/wiki/Glob_(programming)

    Parameters
    ----------
    input_filenames_pattern : str
        A glob pattern for the  displacement files in a format
        readable by SimpleITK.
    epsilon : float
        A lower threshold for the Jacobian, all values below `epsilon`
        will be replaced with `epsilon`.

    Returns
    -------
    sitk.Image
        The geometric mean of the input Jacobian maps.
    """

    total_jacobian = None

    filenames = glob('%s' % input_filenames_pattern)
    n = len(filenames)
    logging.debug('Files to process: %d' % n)

    logging.debug('Starting')
    i = 1
    for f in filenames:
        logging.debug('Processing file %3d/%d %s' % (i, n, f))

        image = sitk.Cast(sitk.ReadImage(f), sitk_float_type)
        image = sitk.Threshold(image, lower=epsilon, upper=1e9, outsideValue=epsilon)
        jacobian = sitk.GetArrayFromImage(image)

        if total_jacobian is None:
            total_jacobian = np.zeros(jacobian.shape, dtype=jacobian.dtype)

        total_jacobian += np.log(jacobian)

        i += 1

    average_jacobian = np.exp(1.0/n * total_jacobian)
    output_jacobian = sitk.GetImageFromArray(average_jacobian)
    output_jacobian.CopyInformation(image)

    return output_jacobian


def jacobian_to_volume_change(jacobian: sitk.Image, epsilon: float = 1e-5) -> sitk.Image:
    r""" Convert a Jacobian map to a volume change map.

    A volume change map is defined as

    .. math::
        VC[f](x) =
        \begin{cases}
            1 - \frac{1}{J[f](x)}  \quad &J[f](x) \in (0,1) \\
            J[f](x) - 1            \quad &J[f](x) \ge 1
        \end{cases}

    Parameters
    ----------
    jacobian : sitk.Image
        Input Jacobian map.
    epsilon : float
        Lower threshold for the Jacobian; any value below
        `epsilon` will be replaced with `epsilon`.

    Returns
    -------
    sitk.Image
        Volume change map associated to the input Jacobian.
    """

    data = sitk.GetArrayFromImage(jacobian)
    processed = np.empty(data.shape, dtype=data.dtype)

    ind_expa = data >= 1.0
    ind_comp = data < 1.0
    ind_sing = data <= epsilon

    data[ind_sing] = epsilon

    processed[ind_expa] = data[ind_expa] - 1.0
    processed[ind_comp] = 1.0 - (1.0 / data[ind_comp])

    result = sitk.GetImageFromArray(processed)
    result.CopyInformation(jacobian)
    return result


def volume_change_to_jacobian(volume_change: sitk.Image) -> sitk.Image:
    r""" Convert a volume change map to a Jacobian map.

    A volume change map is defined as

    .. math::
        VC[f](x) =
        \begin{cases}
            1 - \frac{1}{J[f](x)}  \quad &J[f](x) \in (0,1) \\
            J[f](x) - 1            \quad &J[f](x) \ge 1
        \end{cases}

    Parameters
    ----------
    volume_change : sitk.Image
        Input volume change map.

    Returns
    -------
    sitk.Image
        A Jacobian map associated to the input volume changes.
    """

    data = sitk.GetArrayViewFromImage(volume_change)
    processed = np.empty(data.shape, dtype=data.dtype)

    ind_expa = data >= 0.0
    ind_comp = data < 0.0

    processed[ind_expa] = data[ind_expa] + 1.0
    processed[ind_comp] = -1.0 / (data[ind_comp] - 1.0)

    result = sitk.GetImageFromArray(processed)
    result.CopyInformation(volume_change)
    return result


def deformation_to_displacement(deformation: sitk.Image) -> sitk.Image:
    r""" Convert a deformation field to a displacement field.

    A deformation field :math:`D` is given by the sum of the identity
    transform and a displacement field :math:`d`:

    .. math::
        D(x) = x + d(x)

    Parameters
    ----------
    deformation :
        Input deformation field.

    Returns
    -------
    sitk.Image
        Displacement field associated to the deformation.
    """

    a = sitk.GetArrayFromImage(deformation)

    for x, y, z in np.ndindex(deformation.GetSize()):
        a[z,y,x,0] -= x
        a[z,y,x,1] -= y
        a[z,y,x,2] -= z

    displacement = sitk.GetImageFromArray(a)
    displacement.CopyInformation(deformation)

    return displacement


def compose_displacements(*fields: sitk.Image) -> sitk.Image:
    r""" Compose multiple displacement fields.

    Compute the composition pairwise and iteratively. For a couple
    of displacements :math:`d_1` and :math:`d_2`
    associated to the transforms :math:`f_1` and :math:`f_2`, the
    composition

    .. math::
        (f_2 \circ f_1) (x) = f_2(f_1(x))

    is obtained by resampling :math:`d_2` with :math:`d_1` and then
    summing.

    Parameters
    ----------
    fields : sitk.Image
        Variadic list of displacement fields.

    Returns
    -------
    sitk.Image
        The composition of the input displacement fields.
    """

    fields = list(fields)
    total_field = sitk.Image(fields.pop(0))

    for field in fields:
        resampled_field = sitk.Warp(field,
                                    total_field,
                                    outputSize=total_field.GetSize(),
                                    outputSpacing=total_field.GetSpacing(),
                                    outputOrigin=total_field.GetOrigin(),
                                    outputDirection=total_field.GetDirection())
        resampled_field.CopyInformation(total_field)
        total_field = sitk.Add(total_field, resampled_field)

    return total_field


def decompose_displacements(
        field1: sitk.Image,
        field2: sitk.Image
        ) -> sitk.Image:
    r""" Decompose two displacement fields.

    Given two displacement fields :math:`d_1` and :math:`d_2`
    associated to the transforms :math:`f_1` and :math:`f_2`,
    find a third displacement :math:`d_3` associated to the
    transform :math:`f_3`, such that

    .. math::
        f_1 &= f_3 \circ f_2 \\
        d_1(x) &= d_2(x) + d_3(d_2(x))

    Parameters
    ----------
    field1 : sitk.Image
        Total displacement.

    field2 : sitk.Image
        Component to be decomposed from the total displacement.

    Returns
    -------
    sitk.Image
        A vector image representing a displacement field such
        that its composition with the second argument gives
        the first argument.
    """

    field3 = sitk.Warp(field1 - field2,
                       sitk.InvertDisplacementField(field2),
                       outputSize=field1.GetSize(),
                       outputSpacing=field1.GetSpacing(),
                       outputOrigin=field1.GetOrigin(),
                       outputDirection=field1.GetDirection())
    field3.CopyInformation(field1)
    return field3


def field_zero_padding(
        field  : sitk.Image,
        size_x : Tuple[int,int] = (1,1),
        size_y : Tuple[int,int] = (1,1),
        size_z : Tuple[int,int] = (1,1)
        ) -> sitk.Image:
    r""" Add a zero padding to a vector field.

    Set the zero padding manually, since `sitk.ConstantPad()` does not
    support vector images.

    Parameters
    ----------
    field : sitk.Image
        Input vector field.
    size_x : (int, int)
        Amount of padding at the beginning and end of x direction.
    size_y : (int, int)
        Amount of padding at the beginning and end of y direction.
    size_z : (int, int)
        Amount of padding at the beginning and end of z direction.

    Returns
    -------
    sitk.Image
        A padded vector field.
    """

    a = np.pad(sitk.GetArrayViewFromImage(field),
               (size_x, size_y, size_z, (0,0)),
               'constant',
               constant_values=0.0)

    field_pad = sitk.GetImageFromArray(a)
    field_pad.SetSpacing(field.GetSpacing())
    field_pad.SetOrigin(field.GetOrigin())
    field_pad.SetDirection(field.GetDirection())

    return field_pad


def invert_displacement_padded(field: sitk.Image) -> sitk.Image:
    r""" Invert a displacement field using one voxel of padding in the computation.

    Parameters
    ----------
    field : sitk.Image
        Input displacement field.

    Returns
    -------
    sitk.Image
        The inverse of the displacement, computed under padding.
    """

    inverse = sitk.InvertDisplacementField(field_zero_padding(field))
    inverse = sitk.Crop(inverse, (1,1,1), (1,1,1))
    inverse.CopyInformation(field)

    return inverse


def warp_points_by_displacement(
        points       : np.ndarray,
        displacement : sitk.Image
        ) -> np.ndarray:
    r""" Warp a set of Elastix points by a displacement field.

    Parameters
    ----------
    points : np.ndarray
        A :math:`n \times m` array representing :math:`n` points
        with :math:`m` components.
    displacement : sitk.Image
        Displacement field.

    Returns
    -------
    np.ndarray
        An array representing the warped points.
    """

    data = sitk.GetArrayViewFromImage(displacement)

    (ox, oy, oz) = displacement.GetOrigin()
    (nx, ny, nz) = displacement.GetSize()
    (sx, sy, sz) = displacement.GetSpacing()

    x = np.linspace(ox, (nx-1) * sx, nx)
    y = np.linspace(oy, (ny-1) * sy, ny)
    z = np.linspace(oz, (nz-1) * sz, nz)

    # NOTE: ITK images are indexed as [z,y,x]
    f = scipy.interpolate.RegularGridInterpolator((z, y, x),
                                                  data,
                                                  bounds_error=False)

    p = np.empty(points.shape)
    p[:,0] = points[:,2]
    p[:,1] = points[:,1]
    p[:,2] = points[:,0]

    return points + f(p)


def inverse_consistency_error(
        forward  : sitk.Image,
        backward : sitk.Image,
        mask     : sitk.Image = None
        ) -> Tuple[sitk.Image, float, float]:
    r""" Compute the inverse consistency error (ICE).

    Parameters
    ----------
    forward : sitk.Image
        A displacement from the registration (maps from
        reference space to moving image space).
    backward : sitk.Image
        A displacement from the inverse registration (maps
        from moving image space to reference space)
    mask : sitk.Image
        ROI in the reference image space

    Returns
    -------
    sitk.Image
        Map of the average inverse consistency error magnitude.
    float
        Average inverse consistency error.
    float
        Maximum inverse consistency error.
    """

    composition = compose_displacements(forward, backward)

    if mask is not None:
        composition = sitk.Mask(composition, mask > 0)
        n = np.sum(sitk.GetArrayViewFromImage(mask) > 0)
    else:
        n = reduce(int.__mul__, forward.GetSize())

    vme = np.linalg.norm(sitk.GetArrayViewFromImage(composition), axis=3)

    ic = sitk.GetImageFromArray(vme)
    ic.CopyInformation(forward)

    aic = np.sum(vme) / n
    mic = np.max(vme)

    return ic, aic, mic
