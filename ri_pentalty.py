import math

"""
Two functions used in mass spectrometry / chromatographic data analysis

Both functions calculate a penalty factor meant to be applied to match factors.

Approach is to incorporate retention index information to supplement spectral
data.

Mostly used in isomer differentiaition where spectra may be similar, but
published retention indices may be different
"""

def gaussian_penalty(ri_sample: int, ri_ref: int, shape: float):
    """
    small penalties for small shifts but has exponential scaling
    
    Parameters
    ----------
    ri_sample : int
        retention index of sample
    ri_ref : int
        retention index of reference
    shape : float
        parameter used to control penalty "steepness" or shape of gaussian

    Returns
    -------
    penalty_factor : float
        factor calculated based on difference in retention indices

    """
    delta = abs(ri_sample - ri_ref)
    penalty_factor = math.exp(-1 * shape * delta)
    return penalty_factor

def trapezoidal_penalty(ri_sample: int,
                        ri_ref: int,
                        slope: float,
                        tolerance: float):
    """
    tolerance with penalty-free zone with linear scaling penalties once past
    threshold
    
    Parameters
    ----------
    ri_sample : int
        retention index of sample
    ri_ref : int
        retention index of reference
    slope : float
        parameter used to control penalty "steepness" or slope of line
    tolerance : float
        parameter used to define a penalty-free window with no penalty

    Returns
    -------
    penalty_factor : float
        factor calculated based on difference in retention indices

    """
    delta = abs(ri_sample - ri_ref)
    if delta < tolerance:
        return 1
    else:
        penalty_factor = 1 - (slope * (delta - tolerance))
        return penalty_factor
    
    