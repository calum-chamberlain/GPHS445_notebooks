"""
Useful functions for the Processing and Recording module of the GPHS445
notebooks.

Author: Calum J Chamberlain
Date:   13/04/2020
"""

import numpy as np

from obspy import Trace


AMP_LIMIT = 1e-5  # Limit phase to only plot amplitudes above this.


def plot_fft(
    x: np.array, 
    y: np.array, 
    strict_length: bool = True,
    reconstruct: bool = True,
    log_y: bool = True,
    log_x: bool = True,
    plot_phase: bool = False,
    fft_len: int = None
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
        Whether to use a logathirmic scale on the y (amplitude)
        axis of the amplitude spectra or not.
    log_x:
        Whether to use a logaritmic scale on the x (frequency) axes of the
        spectra or not.
    plot_phase:
        Whether to plot the phase spectrum as well as the amplitude spectrum.
    fft_len:
        Length of fft to use - ignored if strict_len = True.

    Returns
    -------
    A figure of the original time-series, the amplitude spectrum, and the
    time-series computed from inverse FFT of the amplitude spectrum.
    """
    import matplotlib.pyplot as plt
    from scipy import fftpack
    import numpy.fft as fftlib
    
    N = len(y)
    if strict_length:
        fft_len = N
    elif not fft_len:
        # Find the next fast length for the FFT
        fft_len = fftpack.next_fast_len(N)
    
    # Work out the sampling interval
    dt = x[1] - x[0]
    # Check that the data are regularly sampled - note that due to 
    # floating-point rounding errors we need to check that they are close, 
    # rather than exactly the same.
    assert np.allclose(x[1:] - x[0:-1], dt, atol=1e-20), "Data are not regularly sampled"

    # Compute the FFT
    yf = fftlib.fft(y, n=fft_len)
    # Make an array of frequencies that the FFT has been computed for
    # xf = np.linspace(0.0, 1.0 / (2. * dt), int(N / 2))
    xf = np.fft.fftfreq(fft_len, dt)[:fft_len//2]
    # Compute the inverse FFT to check that we recover the data.
    yr = fftlib.ifft(yf)
    # Get the positive magnitude of the FFT - this is the amplitude spectrum
    amplitude_spectrum = np.abs(yf[:fft_len//2])
    # Get phase spectrum - this is the argument of the complex sequence
    phase_spectrum = np.angle(yf[:fft_len//2], deg=True)
    # Wrap phase 0 -> 2 pi
    # phase_spectrum = phase_spectrum % (2 * np.pi)
    # Multiply to normalise amplitude spectra to 1.
    amplitude_spectrum *= 2./fft_len
    # Mask phase spectrum to only show values when there is power in the amplitude spectrum
    # masked_phase = np.ones_like(phase_spectrum) * np.nan
    # masked_phase[amplitude_spectrum > AMP_LIMIT] = phase_spectrum[amplitude_spectrum > AMP_LIMIT]

    # Everything from here is just code used to make the plot.
    nrows = 2
    if reconstruct:
        nrows += 1
    if plot_phase:
        nrows += 1
    
    ts_row, amp_row, phase_row, rs_row = 0, 1, -2, -1
    
    fig, ax = plt.subplots(nrows=nrows)
    
    # Plot the original time-series
    ax[ts_row].plot(x, y, label="Time series")
    ax[ts_row].set_xlabel("Time (s)")
    ax[ts_row].set_ylabel("Amplitude")
    ax[ts_row].autoscale(enable=True, axis='both', tight=True)
    ax[ts_row].legend()
    
    # Plot the amplitude spectrum
    if log_y and log_x:
        ax[amp_row].loglog(xf, amplitude_spectrum, label="Amplitude spectra")
    elif log_x:
        ax[amp_row].semilogx(xf, amplitude_spectrum, label="Amplitude spectra")
    else:
        ax[amp_row].plot(xf, amplitude_spectrum, label="Amplitude spectra")
    ax[amp_row].set_xlabel("Frequency (Hz)")
    ax[amp_row].set_ylabel("Normalised \namplitude")
    ax[amp_row].autoscale(enable=True, axis='both', tight=True)
    ax[amp_row].legend()
    ax[amp_row].grid("on")
    
    if plot_phase:
        # ax[amp_row].get_shared_x_axes().join(ax[amp_row], ax[phase_row])
        ax[amp_row].sharex(ax[phase_row])
        if log_x:
            ax[phase_row].semilogx(xf, phase_spectrum, label="Phase spectra")
        else:
            ax[phase_row].plot(xf, phase_spectrum, label="Phase spectra")
        ax[phase_row].set_xlabel("Frequency (Hz)")
        ax[phase_row].set_ylabel("Phase \n(degrees)")
        ax[phase_row].autoscale(enable=True, axis='x', tight=True)
        ax[phase_row].set_ylim(-180, 180)
        ax[phase_row].legend()
        ax[phase_row].grid("on")
    
    if reconstruct:
        # ax[ts_row].get_shared_x_axes().join(ax[ts_row], ax[rs_row])
        ax[ts_row].sharex(ax[rs_row])
        # Plot the reconstructed time-series
        ax[rs_row].plot(x, np.real(yr)[0:len(x)], label="Reconstructed Time-series")
        ax[ts_row].set_xlabel("Time (s)")
        ax[rs_row].set_ylabel("Amplitude")
        ax[rs_row].autoscale(enable=True, axis='both', tight=True)
        ax[rs_row].legend()

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
        get_window('hann', tr_out.stats.npts))
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
    fig, axes = plt.subplots(ncols=2)
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


def design_filter(sample_rate, data_length, corners, order=4, window_type='butter',
                  filter_type='bandpass', ripple=None, attenuation=None):
    """
    Design a frequency-domain filter.
    
    :type sample_rate: float
    :param sample_rate: Sampling-rate in Hz
    :type data_length: int
    :param data_length: 
        Length of data to apply to - will use the next-fast fft length from this.
    :type corners: array-like
    :param corners: list of corners for filter in order, in Hz
    :type order: int
    :param order: Filter order
    :type window_type: str
    :param window_type: 
        Type of window to use, must be one of:
        'butter' : Butterworth
        'cheby1': Chebyshev I
        'cheby2': Chenyshev II
        'ellip': Cauer/elliptic
        'bessel': Bessel/Thomson
    :type filter_type: str
    :param filter_type: 
        Type of band to use, must be one of: 
        'bandpass', 'lowpass', 'highpass', 'bandstop'
    """
    from scipy import signal, fftpack
    
    nyquist = .5 * sample_rate
    # Check that highpass is usefully less than the nyquist
    if max(corners) > (nyquist * .98):
        raise NotImplementedError(
        "Highcut {0} is higher than Nyquist {1}.".format(
            max(corners), nyquist))
    fft_len = fftpack.next_fast_len(data_length)
    # N is the "order" of the filter, Wn is the filter window.
    b, a = signal.iirfilter(
        N=order, Wn=corners, btype=filter_type, analog=False, ftype=window_type,
        output='ba', rp=ripple, rs=attenuation, fs=sample_rate)
    _, filt = signal.freqz(b, a, worN=fft_len, fs=sample_rate)
    return filt


def filter_and_plot(data, dt, filt):
    """
    Filter data using a simple filter and plot the data, the transfer 
    function and the filtered data. Plots will be in both time and
    frequency domain.
    
    Note: Filter is maximum-phase
    
    :type data: `numpy.ndarray`
    :param data: Data to be filtered
    :type dt: float
    :param dt: Sample-interval (assumed to be in seconds).
    :type filt: `numpy.ndarray`
    :param filt: Frequency-domain representation of filter.
    """ 
    import matplotlib.pyplot as plt
    from scipy import signal, fftpack
    
    N = len(data)
    filt_time = fftpack.ifft(filt)  # Generate a time-domain representation of the filter
    x_time = np.arange(0, N * dt, dt)
    
    fft_len = fftpack.next_fast_len(N)  # Find the next fast length for the FFT
    data_freq = fftpack.fft(data, n=fft_len)
    filtered_freq = data_freq * filt  # Filtering is multiplication in frequency domain
    filtered = fftpack.ifft(filtered_freq)
   
    x_freq = np.linspace(0.0, 1.0 / (2. * dt), int(N / 2))
    
    fig, axes = plt.subplots(nrows=2, ncols=3)
    
    axes[0][0].plot(x_time, data)
    axes[0][0].set_title("Input data")
    axes[0][1].plot(np.arange(0, len(filt_time) * dt, dt), np.real(filt_time))
    axes[0][1].set_title("Filter")
    axes[0][2].plot(x_time, np.real(filtered))
    axes[0][2].set_title("Filtered")
    
    axes[1][0].semilogx(x_freq, 2./N * np.abs(data_freq[:N//2]))    
    axes[1][1].semilogx(x_freq, 2./N * np.abs(filt[:N//2]))
    axes[1][2].semilogx(x_freq, 2./N * np.abs(filtered_freq[:N//2]))
    
    for ax in axes[0]:
        ax.set_xlabel("Time (s)")
        ax.autoscale(enable=True, axis='both', tight=True)
    for ax in axes[1]:
        ax.set_xlabel("Frequency (Hz)")
        ax.autoscale(enable=True, axis='both', tight=True)
    return fig
