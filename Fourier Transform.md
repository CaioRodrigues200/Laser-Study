---
tags:
  - ITQuadrant
  - Faculdade
  - CommunicationSystems
  - TCC
aliases:
  - Transformada de Fourier
  - DFT
  - FFT
  - Fourier Transform
title-class: 0➡️
note-length-class: 0➡️
inlink-class: 5➡️
outlink-class: 3➡️
progressive-summarization-maturity: 0➡️
note-maturity: 2➡️
---

# Table of contents

- [[#Definition]]
- [[#Basic Properties]]
- [[#Transform table]]
- [[#DFT and FFT]]
	- [[#DFT]]
	- [[#FFT]]
		- [[#Where FFT is ill-suited for?]]
- [[#Spectrum Normalizations]]
	- [[#Types of DFTs (or FFTs)]]
	- [[#Consideration about sampling time and frequency bins]]
	- [[#Why to use window functions]]
	- [[#Spectrum Scaling]]
- [[#Sources]]


# Definition
- The Fourier Transform is defined as: $$\mathfrak{F}[f(x)] = F(f) = \int^{\infty}_{-\infty}f(x)~e^{-j2\pi f x}~dx$$

# Basic Properties

- ### Linearity property
	- The linearity property states that, for the signals $x(t)$ and $y(t)$, and its fourier transform $X(f)$ and $Y(f)$: $$\mathfrak{F}[~ax(t)+by(t)~] = aX(f)+bY(f)$$
- ### Duality property
- ### Time shift property
	- The time shift property states that, for a signal $x(t)$ and its fourier transform $X(f)$: $$\mathfrak{F}[~x(t-t_0)~] = e^{-j2\pi f t_0}X(f)$$

- ### Frequency shift property
	- The frequency shift property states that, for a signal $x(t)$ and its fourier transform $X(f)$: $$\mathfrak{F}[~e^{j2\pi f_0 t}x(t)~] = X(f-f_0)$$
	  This is known as the principle of **Modulation**.


- ### Scaling property
	- The time scaling property states that, for a signal $x(t)$ and its fourier transform $X(f)$: $$\mathfrak{F}[x(at)]=\frac{1}{|a|}X\left(\frac{f}{a}\right)$$

- ### Convolution property
- ### [[Autocorrelation]] property
	- The autocorrelation property states that, for a signal $x(t)$ and its fourier transform $X(f)$: $$\mathfrak{F}[R_x(\tau)] = |X(f)|^2 = \int^{\infty}_{-\infty}R_x(\tau)~e^{-j2\pi f \tau}~d\tau$$


# Transform table
**Proofs and deductions:** [[Fourier Transform Table Deductions]]

![[media_c7a_c7afa2a9-47ad-4363-a0c8-49f6b107d1ff_phpZzh3jE 1.png]]

- ### 21. Why is the Fourier TF of a Sum of Deltas also a Sum of Deltas?
	- [Why is the Fourier TF of a Sum of Deltas also a Sum of Deltas? (youtube.com)](https://www.youtube.com/watch?v=ry171Hgvm-8)



# DFT and FFT

## DFT

- The DFT (Discrete Fourier Transform) transforms a sequence of $N$ complex numbers ($x[n]$) into another ($f[k]$). Each $f[k]$ coefficient can be described matematically as: $$f[k] = \sum_{i=0}^{N-1}x[i]~e^{-j2\pi ik/N}$$
- Before any transformation can be made, its necessary to specifiy the **fundamental frequency**. This frequency often is obtained by taking the inverse of the signal period $T$:  $f_0 = 1/T$. The subsequent frequencies are the multiple of it ($2f_0, 3f_0, 4f_0, \text{etc.}$). 
- If its desired an lower frequency than the provided by $1/T$, one can extend the signal using method like [[Zero Padding]], that increases the frequency resolution (but creates some signal artifacts in spectrum, like **lobesides or bumpings**).
- One way to visualize the DFT is by observing a unitary circle:

	- First, we divide de circle in $N$ points, $N$ is the number of sample points of the discrete signal to be transformed in the DFT.
	  ![[Pasted image 20240818235652.png]]
	  
	- Each coefficient $f[k]$ associates to a multiple of the fundamental frequency given by $k \cdot \frac{2\pi}{N}$. Each $f[k]$ is calculated by summing up all multiplications of $x[i]$ and $e^{-j2\pi ik/N}$ complex numbers 
	- The result is the complex number which its magnitudes represents the senoid's magnitude and its phase represents the senoid's phase.
	  ![[Pasted image 20240819001702.png]]
	- Note that a different approach can be taken if we consider: $$e^{-j2\pi ik/N} = cos(2\pi kt/N)-j sin(2\pi kt/N)$$
		- This can generate two coefficients to be calculated, one for cosine and another for sine.

## FFT

### Where FFT is ill-suited for?

- Since the FFT must wait for a large number of samples to begin the transformation, it will be less suited for **low-latency applications**.
- Microcontrollers and low power devices, like noise cancelling headphones, may not have the computational power or memory needed to perform a FFT.


# Spectrum Normalizations

Although the algorithm is matematically well described, on the application field it demands more specifications of what enters ($x[i]$) and what it gives ($y[k]$). In short, a DFT is an operation that takes $N$ real numbers (signal sampling with sampling frequency of $f_s$) and generates $N$ complex numbers. Although, no consideration or limitations was made to the algorithm to care about normalizations and meaning on the results, therefore some considerations and adjustments must be made.

## Types of DFTs (or FFTs)

There are mainly three types of DFTs that differ only by a factor. This factor is applied by convenience depending on the problem, and is very important to know beforehand which type of DFT the signal data package in use calculates: $$\begin{cases} y_k^{(1)} = \sum_{i=0}^{N-1}x[i]~\text{exp}(-j2\pi ik/N)~,~~~k=0~...~N-1 \\ y_k^{(2)} = y_k^{(1)}/\sqrt{N} \\ y_k^{(3)} = y_k^{(1)}/N \end{cases}$$
- The case 1 is the most common case and demands the inverse DFT being divided by the factor $N$.
- The case 2 is called “symmetrical case” because applying successively DFT and inverse DFT reproduces the same signal.
- The case 3 is the simplest form for computation, and demands no division at the inverse DFT.


## Consideration about sampling time and frequency bins

The nyquist frequency determines that the maximum useful frequency is $f_s/2$, and above this there will be non-conclusive frequency bins. If there is $N$ frequency bins, there must be a equally spaced distance between them, called frequency resolution ($f_{res}$), and by the DFT expression its defined by $f_{res} = f_s/N$.

As a consequence of the fact that the input is always real, the output will obey the following relationship:
$$ y_{N-k} = y^*_k $$

This is called the **hermitian symmetry**. Then, considering $N$ even, the DFT will produce only $N/2 + 1$ distinct (independent) complex numbers. These corresponds to the frequencies:
$$ f_k = k \cdot f_{res} = k \cdot \frac{f_s}{N} ~~,~~~~~~\text{for }k=0~...~N/2 $$
The $k = 0$ term corresponds to the DC average of the signal, and it has no imaginary part. The last term $k = N/2$ corresponds to the Nyquist frequency and also has no imaginary part because of complex-conjugate symmetry of the output array ( $y[N/2] = y^*[N/2]$, so the only solution would be a real number).


## Why to use window functions

By taking a time series sampling of length $N$, containing a sinusoidal signal, one can expect that after the DFT operation, in the spectrum will be shown up a single sharp peak in one frequency bin (the frequency of the sinusoidal signal). In reality is more probable that the spectrum will be something like this:
- ![[Pasted image 20250406113935.png | Spectral response of a rectangular window (also a spectrum shape that will probably be reproduced by taking the DFT of a sinusoid without any window function)]]
Apart from the central peak (the desired one), there are many more peaks at regular intervals, these are called *sidelobes* and are undesirable. The reason why it happens is because the DFT implicitly assumes that the input signal is periodic, when in fact there is a big chance it is not. If the frequency of the sinusoid is not a multiple of $f_{res}$, then a mismatch will occur between the first and last sample. This causes a spreading in power all across the spectrum. 

To avoid this problem, is useful to apply a **window function** before appying the DFT on the signal, in a manner that the mismatch between the first and last sample diminishes. This window is basically a sequence of real numbers $\{w_j\}$ where $j = 0~...~N-1$, in which multiplied by the signal sampling will produce $x'_j = x_j \cdot w_j$ as the new input for the DFT. Most windows have a symmetry relationship that implies $w_j = w_{N-j}$. Although, because $j$ maximum value is $N-1$, these windows will not be exactly symmetric.

For example, the Hamming window function is defined as: $$w_j = \frac{1}{2}\left[1-cos\left(\frac{2 \pi \cdot j}{N}\right) \right];~~~~j=0~...~N-1$$ 
- ![[Captura de tela 2025-04-06 120325.png | Left image descbribes the shape of a Hamming window function with a continuous interpolation. Right image shows the Hamming window applied on a signal sampling with $N=8$, and shows that its not exactly symmetric by the fact that there is no $w_8$.]]

## Spectrum Scaling

Until now, the discussed DFT algorithm generates what we call the **linear spectrum**. Lets assume a signal input measured in volts ($V$), the common spectral measure should be given in $V$. By squaring it, another spectrum will be obtained: the **power spectrum**, given in $V^2$.

- The amplitude spectrum given in $V$, will tell the amplitude of each frequency component
- The power spectrum given in $V^2$, will tell the power of each frequency component

Beside these two measures, there are measures that consider only the amplitude/power density in a given equivalent band noise (ENBW). This bandwidth considers the noise influence on the input signal and try to normalize it by taking a kind of average of the window bins. This should be explained better if we imagine a white noise dominating the input signal. This noise should affect only directly all frequency bins across the spectrum and just summing up all the spectrum amplitude by a constant, however, because the sidelobes of the window function, this noise will influence not only the given frequency bin, but all adjacent bins in a given level. 

**Dividing the result by the ENBW should correct this phenomenon**, this band and its normalized version (NENBW) can be estimated by: $$\text{ENBW} = f_s \frac{\sum_{j=0}^{N-1} w_j^2}{\left( \sum_{j=0}^{N-1} w_j \right)^2}~~[\text{samples}/s]~~~,~~~~~~~\text{NENBW} = N \frac{\sum_{j=0}^{N-1} w_j^2}{\left( \sum_{j=0}^{N-1} w_j \right)^2}~~[\text{bins}]$$
If no window is used, is assumed a rectangular window where $w_j = 1$ for every $j$. In this case, $\frac{\sum_{j=0}^{N-1} w_j^2}{\left( \sum_{j=0}^{N-1} w_j \right)^2} = \frac{N}{N^2} = 1/N$. Therefore, in this scenario: $$\text{ENBW} = f_s/N = f_{res}~~~~\text{and}~~~~\text{NENBW} = 1$$
In summary, taking this correction into account, the spectral density measures would be:
- The amplitude spectral density, which is given in $V/\sqrt{Hz}$
- The power spectral density, which is given in $V^2/Hz$

![[Pasted image 20250405201206.png | Naming convention for DFT outputs]]

# Sources

- [Applied DSP No. 9: The z-Domain and Parametric Filter Design](https://youtu.be/xIN5Mnj_MAk?si=48TtIltcvdPx0Ru0)
- [Fourier transform - Wikipedia](https://en.wikipedia.org/wiki/Fourier_transform#Properties_of_the_Fourier_transform)