"""
Code to compute spectral ratios.

Author: Calum J Chamberlain
Date:   13/04/2020
"""

import numpy as np

from obspy import Trace


class Spectra():
    def __init__(self, frequencies, amplitudes):
        assert len(frequencies) == len(amplitudes), "Shapes do not match"
        frequency_interval = frequencies[1] - frequencies[0]
        assert np.allclose(frequencies[1:] - frequencies[0:-1],
                          frequency_interval), "Frequencies is not regularly sampled."
        self._frequency_interval = frequency_interval
        self.frequencies = frequencies
        self.amplitudes = amplitudes

    def __truediv__(self, other):
        assert np.all(self.frequencies == other.frequencies), "Cannot divide with difference frequency sampling"
        out = self.copy()
        out.amplitudes /= other.amplitudes
        return out

    def copy(self):
        from copy import deepcopy
        return deepcopy(self)
        
    def plot(self, log_y: bool = True, log_x: bool = True):
        """
        Plot the spectra

        Returns a Figure.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        if log_y and log_x:
            ax.loglog(self.frequencies, self.amplitudes)
        elif log_x:
            ax.semilogx(self.frequencies, self.amplitudes)
        elif log_y:
            ax.semilogy(self.frequencies, self.amplitudes)
        else:
            ax.plot(self.frequencies, self.amplitudes)
        ax.set_ylabel("Amplitude")
        ax.set_xlabel("Frequency")
        return fig

    def decimate(self, factor: int):
        """
        Smooth the spectra - works in place.

        Parameters
        ----------
        factor:
            Decimation factor.
        """
        self.frequencies = self.frequencies[0::factor]
        self.amplitudes = self.amplitudes[0::factor]
        self._frequency_interval = self.frequencies[1] - self.frequencies[0]


def amplitude_spectra(
    tr: Trace, 
    strict_length: bool = False
) -> Spectra:
    """
    Compute the amplitude spectra for a trace

    Parameters
    ----------
    tr:
        The trace to compute for
    strict_length:
        Whether to compute the FFT for the length of data points (True),
        or the nearest fast length (False).

    Returns
    -------
    The amplitude spectra.
    """
    from scipy import fftpack
    
    N = tr.stats.npts
    if strict_length:
        fft_len = N
    else:
        # Find the next fast length for the FFT
        fft_len = fftpack.next_fast_len(N)

    # Compute the FFT
    yf = fftpack.fft(tr.data, n=fft_len)
    # Make an array of frequencies that the FFT has been computed for
    xf = np.linspace(0.0, 1.0 / (2. * tr.stats.delta), int(N / 2))
    # Compute the inverse FFT to check that we recover the data.
    yr = fftpack.ifft(yf)
    # Get the positive and real component of the FFT - this is the amplitude spectra
    amplitude_spectra = np.abs(yf[:N//2])
    # Multiply to normalise amplitude spectra to 1.
    amplitude_spectra *= 2./N
    return Spectra(frequencies=xf, amplitudes=amplitude_spectra)


def spectral_ratio(
    tr1: Trace, 
    tr2: Trace, 
    decimate: int = None,
    plot: bool = False,
    log_y: bool = False,
    log_x: bool = False,
) -> Spectra:
    """
    Compute the spectral ratio between two traces.

    Parameters
    ----------
    tr1:
        First trace, this will be the numerator
    tr2:
        Second trace, this will be the donominator
    decimate:
        Factor to decimate spectra by prior to computing the ratio.
    plot:
        Whether to plot the resulting ratio or not.
    log_x:
        Whether to use a log axis for the frequency (x) axis.
    log_y:
        Whether to use a log axis for the Y axis.

    Returns
    -------
    Spectra, [Figure]

    Returns a figure of plot = True.
    """
    assert tr1.stats.sampling_rate == tr2.stats.sampling_rate, "Traces must have the same sampling-rate"
    tr1 = tr1.detrend()
    tr2 = tr2.detrend()

    # Both traces need to be the same length to generate the same frequency sampling
    length = max(tr1.stats.npts, tr2.stats.npts)
    if tr1.stats.npts != length:
        tr1.data = np.append([tr1.data, np.zeros(length - tr1.stats.npts)])
    if tr2.stats.npts != length:
        tr2.data = np.append([tr2.data, np.zeros(length - tr2.stats.npts)])
    
    spectra1 = amplitude_spectra(tr1)
    spectra2 = amplitude_spectra(tr2)

    if decimate:
        spectra1.decimate(decimate)
        spectra2.decimate(decimate)
    ratio = spectra1 / spectra2

    if plot:
        fig = ratio.plot(log_x=log_x, log_y=log_y)
        ax = fig.gca()
        ax.set_ylabel("Spectral Ratio")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_title(f"{tr1.id} / {tr2.id}")
        ax.grid()
        return ratio, fig
    return ratio
