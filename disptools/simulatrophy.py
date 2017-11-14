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
