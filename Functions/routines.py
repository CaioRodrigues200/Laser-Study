import numpy as np
import matplotlib.mlab as mlab
from scipy.constants import c

def get_spectrum(x, Fs, Fc, xunits = 'm', yunits = 'dBm', window=mlab.window_none, sides="twosided"):
    spectrum, frequency = mlab.magnitude_spectrum(
        x, Fs=Fs, window=window, sides=sides
    )
    frequency = c/(frequency + Fc) if (xunits=='m') else frequency + Fc
    spectrum = spectrum*np.conj(spectrum)
    if (yunits=='dBm'):
        spectrum = 10*np.log10(1e3*spectrum)

    return frequency, spectrum

def signal_power(x):
    return np.mean(np.abs(x) ** 2)