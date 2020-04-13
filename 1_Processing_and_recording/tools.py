"""
Useful functions for the Processing and Recording module of the GPHS445
notebooks.

Author: Calum J Chamberlain
Date:   13/04/2020
"""

import numpy as np

from obspy import Trace


def plot_fft(
    x: np.array, 
    y: np.array, 
    strict_length: bool = True,
    reconstruct: bool = True,
    log_y: bool = True,
):
    """
    Calculate the FFT of a discretely sampled array and plot the spectra.
    
    Parameters
    ----------
    x:
        Array of sample points - must be regularly sampled, assumed to be in
        units seconds, giving frequencies in Hz.
    y:
        Array of amplitudes - x[n] should be the sample point for y[n].
    strict_length:
        Whether to compute an FFT of length len(y), or the nearest fast 
        length. The default (True) computes an FFT of length len(y)
    reconstruct:
        Whether to plot the reconstructed time-series or not.
    log_y:
        Whether to use a logathirmic scale on the y (amplitude) axis of the
        amplitude spectra or not.

    Returns
    -------
    A figure of the original time-series, the amplitude spectrum, and the
    time-series computed from inverse FFT of the amplitude spectrum.
    """
    import matplotlib.pyplot as plt
    from scipy import fftpack
    
    N = len(y)
    if strict_length:
        fft_len = N
    else:
        # Find the next fast length for the FFT
        fft_len = fftpack.next_fast_len(N)
    
    # Work out the sampling interval
    dt = x[1] - x[0]
    # Check that the data are regularly sampled - note that due to 
    # floating-point rounding errors we need to check that they are close, 
    # rather than exactly the same.
    assert np.allclose(x[1:] - x[0:-1], dt, atol=1e-20), "Data are not regularly sampled"

    # Compute the FFT
    yf = fftpack.fft(y, n=fft_len)
    # Make an array of frequencies that the FFT has been computed for
    xf = np.linspace(0.0, 1.0 / (2. * dt), int(N / 2))
    # Compute the inverse FFT to check that we recover the data.
    yr = fftpack.ifft(yf)
    # Get the positive and real component of the FFT - this is the amplitude spectra
    amplitude_spectra = np.abs(yf[:N//2])
    # Multiply to normalise amplitude spectra to 1.
    amplitude_spectra *= 2./N

    # Everything from here is just code used to make the plot.
    if reconstruct:
        nrows = 3
    else:
        nrows = 2
    fig, ax = plt.subplots(nrows=nrows, figsize=(15, 8))
    
    # Plot the original time-series
    ax[0].plot(x, y, label="Time series")
    ax[0].set_xlabel("Time (s)")
    ax[0].set_ylabel("Amplitude")
    ax[0].autoscale(enable=True, axis='both', tight=True)
    ax[0].legend()
    
    # Plot the amplitude spectrum
    if log_y:
        ax[1].loglog(xf, amplitude_spectra, label="Amplitude spectra")
    else:
        ax[1].semilogx(xf, amplitude_spectra, label="Amplitude spectra")
    ax[1].set_xlabel("Frequency (Hz)")
    ax[1].set_ylabel("Normalised amplitude")
    ax[1].autoscale(enable=True, axis='both', tight=True)
    ax[1].legend()
    
    if reconstruct:
        # Plot the reconstructed time-series
        ax[2].plot(x, np.real(yr)[0:len(x)], label="Reconstructed Time-series")
        ax[2].set_xlabel("Time (s)")
        ax[2].set_ylabel("Amplitude")
        ax[2].autoscale(enable=True, axis='both', tight=True)
        ax[2].legend()

    return fig


def resample_and_plot(tr: Trace, sampling_rate: float):
    """
    Simplified (made less general) from obspy Trace.resample.
    
    Uses a frequency domain method to resample, and a hanning window.
    
    Parameters
    ----------
    tr:
        Obspy Trace to resample. Works on a copy of the trace. 
        Your data are safe here.
    sampling_rate:
        Desired sampling rate in Hz.

    Returns
    -------
    Resampled Trace, Figure.
    """ 
    from scipy.signal import get_window
    from scipy.fftpack import rfft, irfft
    import matplotlib.pyplot as plt

    factor = tr.stats.sampling_rate / float(sampling_rate)
    
    # Copy the trace and work on this copy
    tr_out = tr.copy()

    # Copy things for plotting
    data_in = tr.data
    max_time = tr.stats.npts * tr.stats.delta
    dt = tr.stats.delta
    N = tr.stats.npts
    # resample in the frequency domain. Make sure the byteorder is native.
    x = rfft(tr_out.data.newbyteorder("="))
    # Cast the value to be inserted to the same dtype as the array to avoid
    # issues with numpy rule 'safe'.
    x = np.insert(x, 1, x.dtype.type(0))
    if tr_out.stats.npts % 2 == 0:
        x = np.append(x, [0])
    x_r = x[::2]
    x_i = x[1::2]

    # Multiply by a hanning window to stabilise the interpolation
    large_w = np.fft.ifftshift(
        get_window('hanning', tr_out.stats.npts))
    x_r *= large_w[:tr_out.stats.npts // 2 + 1]
    x_i *= large_w[:tr_out.stats.npts // 2 + 1]

    # interpolate
    num = int(tr_out.stats.npts / factor)
    df = 1.0 / (tr_out.stats.npts * tr_out.stats.delta)
    d_large_f = 1.0 / num * sampling_rate
    f = df * np.arange(0, tr_out.stats.npts // 2 + 1, dtype=np.int32)
    n_large_f = num // 2 + 1
    large_f = d_large_f * np.arange(0, n_large_f, dtype=np.int32)
    large_y = np.zeros((2 * n_large_f))
    large_y[::2] = np.interp(large_f, f, x_r)
    large_y[1::2] = np.interp(large_f, f, x_i)

    large_y = np.delete(large_y, 1)
    if num % 2 == 0:
        large_y = np.delete(large_y, -1)
    tr_out.data = irfft(large_y) * (float(num) / float(tr_out.stats.npts))
    tr_out.stats.sampling_rate = sampling_rate
    
    # Plot
    fig, axes = plt.subplots(ncols=2, figsize=(14, 6))
    axes[0].plot(np.arange(0, max_time, dt), data_in)
    axes[1].semilogx(
        np.linspace(0.0, 1.0 / (2. * dt), int(N / 2)),
        2./N * np.abs(x[:N//2]), label="Original")
    axes[1].semilogx(
        np.linspace(dt, dt + 1.0 / (2. * tr_out.stats.delta), num // 2),
        2./N * np.abs(large_y[:num//2]), label="Resampled")
    axes[0].plot(np.arange(0, max_time, 1 / sampling_rate), tr_out.data)
    axes[1].legend()
    axes[0].set_xlabel("Time (s)")
    axes[1].set_xlabel("Frequency (Hz)")
    return tr_out, fig