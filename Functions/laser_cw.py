import numpy as np
from routines import get_spectrum
from scipy.constants import c

class parameters():
    pass

class laser():
    def __init__(self, parameters, numberOfSamples:int, Fs:float) -> None:
        # Laser parameters from user
        self.wavelength = getattr(parameters, "wavelength", 1550e-9) #m
        self.power  = getattr(parameters, "power", 1e-3) #W
        self.linewith = getattr(parameters, "linewith", 100e6) #hz
        self.initialPhase = getattr(parameters, "initialPhase", 0) #rad
        self.polarization_alfa = getattr(parameters, "polarization_alfa", 0) #rad
        self.polarization_ellp = getattr(parameters, "polarization_ellp", 0) #rad
        self.Ns = numberOfSamples
        self.Fs = Fs
        # Laser parameters calculated from others parameters
        self.k = self.set_k()
        self.theta = self.set_theta()
        self.laser_phase = self.set_laser_phase()

    def set_k(self) -> float:
        a = np.tan(2*self.polarization_alfa)
        b = np.sin(2*self.polarization_ellp)
        if (self.polarization_alfa >= 0):
            return 0.5 - 0.5*np.sqrt((1-b**2)/(1+a**2))
        return 0.5 + 0.5*np.sqrt((1-b**2)/(1+a**2))
    
    def set_theta(self) -> float:
        b = np.sin(2*self.polarization_ellp)
        if self.k in [1, 0]:
            return 0
        else:
            return np.arcsin(b/(2*np.sqrt(self.k*(1-self.k))))
        
    def set_laser_phase(self) -> np.ndarray:
        mean_dphase = 0
        variance_dphase = (2*np.pi*self.linewith*1/self.Fs)
        laser_dphase = mean_dphase + np.random.randn(self.Ns, 1) * np.sqrt(variance_dphase)
        return self.initialPhase + np.cumsum(laser_dphase)
    
    @property
    def e_x(self) -> np.ndarray:
        if self.k == 1:
            return np.zeros(self.Ns)
        temp = np.ones(self.Ns)*np.sqrt(self.power)*np.sqrt(1-self.k)
        return temp * np.exp(1j*self.laser_phase)

    @property
    def e_y(self) -> np.ndarray:
        if self.k == 0:
            return np.zeros(self.Ns)
        temp = np.ones(self.Ns)*np.sqrt(self.power)*np.sqrt(self.k)*np.exp(1j*self.theta)
        return temp * np.exp(1j*self.laser_phase)

    def plot(self, ax, xunits:str = 'm', yunits:str = 'dBm') -> None:
        # X plot
        freq, spectrum = get_spectrum(self.e_x, Fs=self.Fs, Fc=c/self.wavelength, xunits=xunits, yunits=yunits)
        ax.plot(freq, spectrum, 'k', linewidth = 2, label = 'X pol.')
        # Y plot
        freq, spectrum = get_spectrum(self.e_y, Fs=self.Fs, Fc=c/self.wavelength, xunits=xunits, yunits=yunits)
        ax.plot(freq, spectrum, 'r', linewidth = 2, label = 'Y pol.', alpha=0.5)
        ax.set_xlabel(f'Frequency [{xunits}]')
        ax.set_ylabel(f'Power [{yunits}]')
        ax.legend()
        ax.grid()

'''
def laser_rate_equations(self, t, x):
    # Carrier density
    n = carrier_density(self, x)
    # Foton density
    s = photon_density(self, x)
    # Optical phase
    phi = optical_phase(self, x)
    return np.asarray([n,s,phi])

def carrier_density(self,):
    F_N = np.randn*np.sqrt(2*self.v*self.n_sd/self.tau_n*(self.beta*self.v*self.s_sd+1))
    return I_t/(q*self.v)-y(2)/self.tau_n-self.vg*self.a0*y(1)*(y(2)-self.n_t)/(1+self.epsilon*y(1)) + F_N

def photon_density(self,):
    F_S = np.randn*np.sqrt(2*self.beta*self.v*self.n_sd*(self.v*self.s_sd+1)^3/self.tau_n)    
    return (self.gamma * self.a0 * self.vg * ((y(2)-self.n_t)/(1+self.epsilon*y(1)))-1/self.tau_p)*y(1)+self.beta*self.gamma*y(2)/self.tau_n + F_S
    
def optical_phase(self, ):
    F_phi = np.randn*np.sqrt(self.gamma*self.vg*self.sigma*(y(2)-self.n_t)/y(1))
    return self.alpha/2 * (self.gamma*self.vg*self.a0*(y(2)-self.n_t)-1/self.tau_p) + F_phi

def plot(self):
    P = 0.5 * Y(:,1) * laserP.Vactv * laserP.Mune * h * laserP.v ./ (laserP.Gamma .* laserP.Tphot);
    C = 1/(2.*pi).*(laserP.Alpha./2).*(laserP.Gamma.*laserP.vg.*laserP.a0.*(Y(:,2)-laserP.Nnull)-1/laserP.Tphot)+randn(length(Y(:,2)),1).*sqrt(laserP.Gamma*laserP.Vg*(Y(:,2)-laserP.Nnull) ./ Y(:,1));


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import random

    N = 10    
    i_t = random.randint(0,2**N - 1)
    laser = laser_dfb()

def __init__:
    I_bias = laserP.Bias;
    Im     = laserP.CPeak;
    T  = 417e-12;
    tr = 100e-12;
    if t < 2*T
        if t < 0
            Ip = 0;
        elseif t < T
            Ip = Im * (1 - exp(-(t./tr).^2));
        else
            Ip = Im * exp(-(t-T).^2./tr.^2);
        end
    else
        temp = t - 2*T;
        if temp < 0
            Ip = 0;
        elseif temp < T
            Ip = Im * (1 - exp(-(temp./tr).^2));
        else
            Ip = Im * exp(-(temp-T).^2./tr.^2);
        end
    end
    I_t = I_bias + Ip;
    plot(t/T, I_t*1000,'r-o');
    hold on;
'''