# -*- coding: utf-8 -*-
import numpy as np
from scipy.fftpack import dct

def costf(f, fac=True):
    """
    This routine transform an input array from real to Chebyshev space

    :param f: the input array
    :type f: numpy.ndarray
    :param fac: normalisation factor is used
    :type f: bool
    :returns: a transformed array
    :rtype: numpy.ndarray
    """
    if fac:
        norm = np.sqrt(0.5/(f.shape[-1]-1))
    else:
        norm = 1.

    fhat = norm*dct(f, type=1, axis=-1)

    return fhat


def get_dr(f):
    """
    This routine calculates the first radial derivative of a input array using
    Chebyshev recurrence relation.

    :param f: the input array
    :type f: numpy.ndarray
    :returns: the radial derivative of f
    :rtype: numpy.ndarray
    """
    Nr = f.shape[-1]
    fhat = costf(f)

    df = np.zeros_like(fhat)
    df[..., -1] = 0.
    df[..., -2] = (Nr-1)*fhat[..., -1]

    for i in range(Nr-3, -1, -1):
        df[..., i] = df[..., i+2]+2.*(i+1)*fhat[..., i+1]

    df[..., :] = 2.*df[..., :]

    df = costf(df)

    return df
