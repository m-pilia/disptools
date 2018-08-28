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
        radius: int = 8,
        kernel: int = sitk.sitkBall,
        ) -> sitk.Image:
    r""" Convert a binary mask to a Simul@trophy mask.

    The mask used by Simul@trophy has five labels:
      - 0: skull
      - 1: cerebro-spinal fluid (CSF)
      - 2: gray matter
      - 3: white matter
      - 4: falx cerebri

    This function takes as input a binary mask, and returns another mask
    in the Simul@atrophy format, where the ROI of the original mask is
    mapped to white matter, a surrounding region of CSF is created
    around it, and the remaining is set to skull.

    Parameters
    ----------
    image : sitk.Image
        Input binary mask.

    radius : int
        Radius for the dilation, determines the amount of CSF
        surrounding the ROI.

    kernel : int
        Kernel used for the dilation, among the values in
        `itk::simple::KernelEnum`.

    Returns
    -------
    sitk.Image
        A Simul@atrophy mask constructed from the input mask.
    """

    image = image > 0
    dilated = sitk.BinaryDilate(image, radius, kernel)
    return 2 * image + dilated

