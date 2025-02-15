# Python Functions Used

## Scipy
- from scipy.ndimage import **shift**(input, shift, output=None, order=3, mode='constant', cval=0.0, prefilter=True)

    * Shift an array.
    * The array is shifted using spline interpolation of the requested order. Points outside the boundaries of the input are filled according to the given mode.

## Matplotlib
- from matplotlib.mlab import **magnitude_spectrum**(x, Fs=None, window=None, pad_to=None, sides=None)

    * Compute the magnitude (absolute value) of the frequency spectrum of x. Data is padded to a length of *pad_to* and the windowing function *window* is applied to the signal.
    * The OptiCommPy function **get_spectrum** uses this function and adapts it to complex input values (double sided) and converts to show power in y-axis