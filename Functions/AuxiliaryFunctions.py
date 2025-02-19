import numpy as np

def rc_filter(R, C, Fs, filter_type):
    """Implements a digital filter based on a RC circuit, converting an analog RC filter characteristics into a digital RC filter using the bilinear transformation (A).    
       The outputs of this functions are coefficients and should be used with **signal.lfilter_zi(b, a)** and **signal.lfilter(b, a, sigTx)**   

        Parameters
        ----------
        R : float
            Resistance value, in Ohms.
        C : float
            Capacitance value, in F.
        Fs : float
            Sampling frequency in Samp/s.
        filter_type : {'low', 'high'}
            Filter type.

    """

    # Analog cut frequency cannot be greater than Fs/2, by Nyquist criteria.
    if (1/(2*np.pi*R*C) > Fs/2):
        ('This analog filter cannot be realized with this sample rate')    
    # Default to allpass if invalid type is selected

    # Filter coefficients
    b = [1, 0]  # Numerator coefficients
    a = [1, 0]  # Denominator coefficients
    # Constants
    RC = R * C
    T  = 1 / Fs
    # Analog Cutoff Fc
    w = 1 / (RC)
    # Prewarped coefficient for Bilinear transform
    A = 1 / (np.tan((w*T) / 2))
    if(filter_type=='high'):
        b[0] = (A)     / (1 + A)
        b[1] = -b[0]
        a[1] = (1 - A) / (1 + A)
    if(filter_type=='low'):
        b[0] = (1)     / (1 + A)
        b[1] = b[0]
        a[1] = (1 - A) / (1 + A)
    return b,a
