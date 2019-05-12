"""
Correlate_template function borrowed from obspys 
forthcoming (1.1.2) release.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from future.utils import native_str

from bisect import bisect_left
from copy import copy
import ctypes as C  # NOQA
from distutils.version import LooseVersion
import warnings

import numpy as np
import scipy

from obspy import Stream, Trace
from obspy.core.util.misc import MatplotlibBackend
from obspy.signal.headers import clibsignal
from obspy.signal.invsim import cosine_taper


def _pad_zeros(a, num, num2=None):
    """Pad num zeros at both sides of array a"""
    if num2 is None:
        num2 = num
    hstack = [np.zeros(num, dtype=a.dtype), a, np.zeros(num2, dtype=a.dtype)]
    return np.hstack(hstack)


def _call_scipy_correlate(a, b, mode, method):
    """
    Call the correct correlate function depending on Scipy version and method.
    """
    if LooseVersion(scipy.__version__) >= LooseVersion('0.19'):
        cc = scipy.signal.correlate(a, b, mode=mode, method=method)
    elif method in ('fft', 'auto'):
        cc = scipy.signal.fftconvolve(a, b[::-1], mode=mode)
    elif method == 'direct':
        cc = scipy.signal.correlate(a, b, mode=mode)
    else:
        msg = "method keyword has to be one of ('auto', 'fft', 'direct')"
        raise ValueError(msg)
    return cc


def _window_sum(data, window_len):
    """Rolling sum of data."""
    window_sum = np.cumsum(data)
    # in-place equivalent of
    # window_sum = window_sum[window_len:] - window_sum[:-window_len]
    # return window_sum
    np.subtract(window_sum[window_len:], window_sum[:-window_len],
                out=window_sum[:-window_len])
    return window_sum[:-window_len]


def correlate_template(data, template, mode='valid', normalize='full',
                       demean=True, method='auto'):
    """
    Normalized cross-correlation of two signals with specified mode.

    If you are interested only in a part of the cross-correlation function
    around zero shift consider using function
    :func:`~obspy.signal.cross_correlation.correlate` which allows to
    explicetly specify the maximum shift.

    :type data: :class:`~numpy.ndarray`, :class:`~obspy.core.trace.Trace`
    :param data: first signal
    :type template: :class:`~numpy.ndarray`, :class:`~obspy.core.trace.Trace`
    :param template: second signal to correlate with first signal.
        Its length must be smaller or equal to the length of ``data``.
    :param str mode: correlation mode to use.
        It is passed to the used correlation function.
        See :func:`scipy.signal.correlate` for possible options.
        The parameter determines the length of the correlation function.
    :param normalize:
        One of ``'naive'``, ``'full'`` or ``None``.
        ``'full'`` normalizes every correlation properly,
        whereas ``'naive'`` normalizes by the overall standard deviations.
        ``None`` does not normalize.
    :param demean: Demean data beforehand. For ``normalize='full'`` data is
        demeaned in different windows for each correlation value.
    :param str method: Method to use to calculate the correlation.
         ``'direct'``: The correlation is determined directly from sums,
         the definition of correlation.
         ``'fft'`` The Fast Fourier Transform is used to perform the
         correlation more quickly.
         ``'auto'`` Automatically chooses direct or Fourier method based on an
         estimate of which is faster. (Only availlable for SciPy versions >=
         0.19. For older Scipy version method defaults to ``'fft'``.)

    :return: cross-correlation function.

    .. note::
        Calling the function with ``demean=True, normalize='full'`` (default)
        returns the zero-normalized cross-correlation function.
        Calling the function with ``demean=False, normalize='full'``
        returns the normalized cross-correlation function.

    .. rubric:: Example

    >>> from obspy import read
    >>> data = read()[0]
    >>> template = data[450:550]
    >>> cc = correlate_template(data, template)
    >>> index = np.argmax(cc)
    >>> index
    450
    >>> round(cc[index], 9)
    1.0
    """
    # if we get Trace objects, use their data arrays
    if isinstance(data, Trace):
        data = data.data
    if isinstance(template, Trace):
        template = template.data
    data = np.asarray(data)
    template = np.asarray(template)
    lent = len(template)
    if len(data) < lent:
        raise ValueError('Data must not be shorter than template.')
    if demean:
        template = template - np.mean(template)
        if normalize != 'full':
            data = data - np.mean(data)
    cc = _call_scipy_correlate(data, template, mode, method)
    if normalize is not None:
        tnorm = np.sum(template ** 2)
        if normalize == 'naive':
            norm = (tnorm * np.sum(data ** 2)) ** 0.5
            if norm <= np.finfo(float).eps:
                cc[:] = 0
            elif cc.dtype == float:
                cc /= norm
            else:
                cc = cc / norm
        elif normalize == 'full':
            pad = len(cc) - len(data) + lent
            if mode == 'same':
                pad1, pad2 = (pad + 2) // 2, (pad - 1) // 2
            else:
                pad1, pad2 = (pad + 1) // 2, pad // 2
            data = _pad_zeros(data, pad1, pad2)
            # in-place equivalent of
            # if demean:
            #     norm = ((_window_sum(data ** 2, lent) -
            #              _window_sum(data, lent) ** 2 / lent) * tnorm) ** 0.5
            # else:
            #      norm = (_window_sum(data ** 2, lent) * tnorm) ** 0.5
            # cc = cc / norm
            if demean:
                norm = _window_sum(data, lent) ** 2
                if norm.dtype == float:
                    norm /= lent
                else:
                    norm = norm / lent
                np.subtract(_window_sum(data ** 2, lent), norm, out=norm)
            else:
                norm = _window_sum(data ** 2, lent)
            norm *= tnorm
            if norm.dtype == float:
                np.sqrt(norm, out=norm)
            else:
                norm = np.sqrt(norm)
            mask = norm <= np.finfo(float).eps
            if cc.dtype == float:
                cc[~mask] /= norm[~mask]
            else:
                cc = cc / norm
            cc[mask] = 0
        else:
            msg = "normalize has to be one of (None, 'naive', 'full')"
            raise ValueError(msg)
    return cc