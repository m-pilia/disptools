import subprocess
import tempfile

import SimpleITK as sitk
import numpy as np
from disptools import *


def jacobian_to_atrophy_map(image: sitk.Image) -> sitk.Image:
    r""" Convert a Jacobian map to a atrophy map.

    The atrophy rate :math:`a` is defined as the pointwise percentage of
    volume loss, so :math:`a = -(J-1)` (where :math:`J` is the Jacobian).

    Parameters
    ----------
    image : sitk.Image
        Input image.

    Returns
    -------
    sitk.Image
        Corresponding atrophy map.
    """

    return -(image - 1.0)


def mask_to_simulatrophy_mask(
        image: sitk.Image,
        radius: int = None,
        kernel: int = sitk.sitkBall,
        ) -> sitk.Image:
    r""" Convert a binary mask to a Simul\@atrophy mask.

    The mask used by Simul\@atrophy has five labels:
      - 0: skull
      - 1: cerebro-spinal fluid (CSF)
      - 2: gray matter
      - 3: white matter
      - 4: falx cerebri

    This function takes as input a binary mask, and returns another mask
    in the Simul\@atrophy format, where the ROI of the original mask is
    mapped to white matter, a surrounding region of CSF is created
    around it, and the remaining is set to skull.

    Parameters
    ----------
    image : sitk.Image
        Input binary mask.

    radius : int
        Radius for the dilation, determines the amount of CSF
        surrounding the ROI. If `None`, all the volume outside the ROI
        is set to CSF.

    kernel : int
        Kernel used for the dilation, among the values in
        `itk::simple::KernelEnum`_.

        .. _itk::simple::KernelEnum: https://itk.org/SimpleITKDoxygen/html/namespaceitk_1_1simple.html#a38998f2c7b469b1ad8e337a0c6c0697b

    Returns
    -------
    sitk.Image
        A Simul\@atrophy mask constructed from the input mask.
    """

    image = image > 0
    dilated = sitk.BinaryDilate(image, radius, kernel) if radius is not None else 1
    return 2 * image + dilated


def run(
        jacobian   : sitk.Image,
        mask       : sitk.Image,
        scale      : Tuple[float, float, float] = None,
        sigma      : float = 2.0,
        lame       : Tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0),
        origin     : Tuple[int, int, int] = None,
        size       : Tuple[int, int, int] = None,
        executable : str = None,
        ) -> sitk.Image:
    r""" Wrapper around Simul\@atrophy.

    Use Simul\@atrophy [7]_ to generate a displacement that realises the
    volume changes prescribed by the input Jacobian map. Optionally, the
    function can operate on a downsampled version of the image.

    .. note::
        Requires a working installation of Simul\@atrophy.

    References
    ----------
    .. [7] Khanal, Bishesh, Nicholas Ayache, and Xavier Pennec.
           "Simulating Longitudinal Brain MRIs with Known Volume Changes
           and Realistic Variations in Image Intensity." Frontiers in
           neuroscience 11 (2017): 132.

    Parameters
    ----------
    jacobian : sitk.Image
        Jacobian map.

    mask : sitk.Image
        Binary mask marking the ROI whose Jacobian shall be matched.

    scale : Tuple[float, float, float]
        If not ``None``, operate on an image downsampled by dividing
        each size component by the given factors.

    sigma : float
        Amount of smoothing prior to resampling. Only relevant when
        ``scale`` is not ``None``.

    lame : Tuple[float, float, float, float]
        Lamé parameters, in the following order:
          * :math:`\mu` within the ROI
          * :math:`\mu` outside the ROI
          * :math:`\lambda` within the ROI
          * :math:`\lambda` outside the ROI

    origin : Tuple[int, int, int]
        If not `None` then only a region of the image is processed,
        defined by origin and size. Requires to specify ``size``.

    size : Tuple[int, int, int]
        If not `None` then only a region of the image is processed,
        defined by origin and size. Requires to specify ``origin``.

    executable : str
        Path to the Simul\@atrophy executable. If `None`, the executable
        is searched within the system path.

    Returns
    -------
    sitk.Image
        A displacement field realising the volume changes of the given
        Jacobian map.
    """
    if executable is None:
        executable = 'simul_atrophy'

    with tempfile.TemporaryDirectory() as tmpdir:

        atrophy_file = os.path.join(tmpdir, 'atrophy_map.mha')
        mask_file = os.path.join(tmpdir, 'mask.mha')

        args = [
            executable,
            '-parameters', ','.join([str(p) for p in lame]),
            '-boundary_condition', 'dirichlet_at_walls',
            '--relax_ic_in_csf',
            '-atrophyFile', atrophy_file,
            '-maskFile', mask_file,
            '-imageFile', mask_file, # dummy parameter
            '--invert_field_to_warp',
            '-numOfTimeSteps', '1',
            '-resPath', os.path.join(tmpdir, ''),
            '-resultsFilenamesPrefix', 'out_',
        ]

        if origin is not None and size is not None:
            args += [
                '-domainRegion', ' '.join([str(x) for x in list(origin) + list(size)]),
            ]

        mask_s = mask_to_simulatrophy_mask(mask, radius=None)
        atrophy = jacobian_to_atrophy_map(jacobian)
        if scale is not None:
            img_origin = jacobian.GetOrigin()
            img_size = [x // s for x, s in zip(jacobian.GetSize(), scale)]
            img_spacing = [x * s for x, s in zip(jacobian.GetSpacing(), scale)]
            mask_s = sitk.Resample(mask_s, img_size, sitk.Transform(), sitk.sitkNearestNeighbor, img_origin, img_spacing)
            atrophy = sitk.SmoothingRecursiveGaussian(atrophy, sigma)
            atrophy = sitk.Resample(atrophy, img_size, sitk.Transform(), sitk.sitkBSpline, img_origin, img_spacing)
        sitk.WriteImage(mask_s, mask_file)
        sitk.WriteImage(atrophy, atrophy_file)
        del mask_s
        del atrophy

        proc = subprocess.Popen(args, cwd=tmpdir)
        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError('Execution failed with status {}'.format(proc.returncode))

        displacement = sitk.ReadImage(os.path.join(tmpdir, 'out_T1vel.nii.gz'))
        if scale is not None:
            displacement = sitk.SmoothingRecursiveGaussian(displacement, sigma)
            displacement = sitk.Resample(displacement, jacobian, sitk.Transform(), sitk.sitkBSpline)

        return displacement

