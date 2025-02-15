# References:
# [1] - J. C. Cartledge and R. C. Srinivasan, "Extraction of DFB laser rate equation parameters for system simulation purposes," in Journal of Lightwave Technology, vol. 15, no. 5, pp. 852-860, 1997.
# [2] - I. Fatadin, D. Ives and M. Wicks, "Numerical simulation of intensity and phase noise from extracted parameters for CW DFB lasers," in IEEE Journal of Quantum Electronics, vol. 42, no. 9, pp. 934-941, 2006.
#
from numba import jit
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
from scipy.integrate import solve_ivp,odeint
from scipy.constants import c, h, e

# Routines for DFB laser simulation.
class laser_dfb():
    def __init__(self, parameters):
        ## add noise terms
        self.noise_terms = getattr(parameters, "noise_terms", False)
        # active layer volume
        self.v = getattr(parameters, "v", 1.5e-10)
        # electron lifetime
        self.tau_n = getattr(parameters, "tau_n", 1.0e-9)
        # active layer gain coefficient
        self.a0 = getattr(parameters, "a0", 2.5e-16)
        # group velocity
        self.vg = getattr(parameters, "vg", 8.5e+9)
        # carrier density at transparency
        self.n_t = getattr(parameters, "n_t", 1.0e+18)
        # gain compression factor
        self.epsilon = getattr(parameters, "epsilon", 1.0e-17)
        # mode confinement factor
        self.gamma = getattr(parameters, "gamma", 0.4)
        # photon lifetime
        self.tau_p = getattr(parameters, "tau_p", 3.0e-12)
        # fraction of spontaneous emission coupling
        self.beta = getattr(parameters, "beta", 3.0e-5)
        # linewidth enchancement factor
        self.alpha = getattr(parameters, "alpha", 5)
        # gain cross section
        self.sigma = getattr(parameters, "sigma", 2e-20)
        # i_bias
        self.i_bias = getattr(parameters, 'i_bias', 0.078)
        # i_max
        self.i_max  = getattr(parameters, 'i_max', 0.188)
        # total differential quantum efficiency
        self.eta_0 = getattr(parameters, 'eta_0', 0.4)
        # wavelength
        self.lmbd = getattr(parameters, 'lmbd', 1300e-9)
        # frequency
        self.freq_0  = c/self.lmbd
        # threshold current
        self.ith = e/self.tau_n*self.v*(self.n_t + 1/(self.gamma * self.tau_p * self.vg * self.a0))
        
    def get_current(self, t):
        return np.real(self.current[np.argwhere(t >= self.t)[-1]])
    
    def set_current(self, t, current):
        # set variables
        self.current = self.i_max * current + self.i_bias
        self.t = t
        self.t_step = t[1]-t[0]
        
    def set_laser_dynamic_response(self, y):
        # get variables
        self.N = y[0,:]
        self.S = y[1,: ]
        self.phase = y[2,:]
        # calculate others parameters
        self.power = (self.S/2) * (self.v * h * self.freq_0 * self.eta_0)/(self.gamma * self.tau_p)
        self.chirp = 1/(2*np.pi) * (self.alpha / 2) * (self.gamma * self.vg * self.a0 * (self.N - self.n_t) - 1/self.tau_p) # np.diff(y[2,:],prepend=y[2,0])/self.t_step
        self.e_out = np.sqrt(np.real(self.power)) * np.exp(1j * self.phase)

    def get_initial_conditions(self):
        y0 = np.zeros([3])        
        y0[1] = self.gamma * self.tau_p/(self.v * e) * (self.get_current(0)-self.ith)
        y0[0] = self.n_t + (1 + self.epsilon * y0[1]) / (self.vg * self.a0 * self.gamma * self.tau_p)
        y0[2] = (self.alpha / 2) * (self.gamma * self.vg * self.a0 * (y0[0] - self.n_t) - 1 / self.tau_p)
        return y0

    def solve(self, t, current):
        self.set_current(t, current)
        # Initial conditions        
        sol = solve_ivp(
            self.laser_rate_equations,
            t_span=[t[0], t[-1]],
            y0=self.get_initial_conditions(),
            method='RK45',
            t_eval = t,
            dense_output=True,
            #rtol = 1e-4,
            #atol = 1e-6,
        )
        self.set_laser_dynamic_response(sol['y'])
        return sol

    def laser_rate_equations(self, t, y):
        dy = np.zeros([3])        
        # Carrier density
        dy[0] = self.carrier_density(t,y)
        # Foton density
        dy[1] = self.photon_density(y)
        # Optical phase
        dy[2] = self.optical_phase(y)
        return self.add_noise_rate_equations(y,dy)

    def add_noise_rate_equations(self,y,dy):
        # noise terms
        dn = np.zeros(3)
        if self.noise_terms:
            dn = self.laser_noise_sources(y)
        return (dy + dn)

    def laser_noise_sources(self, y):
        # noise term 
        # diffusion coefficient - see ref. [2]
        dss = (self.beta*y[0]*y[1]/self.tau_n)
        #dnn = (y[0]/self.tau_n)*(1+self.beta*y[1])
        dpp = (self.beta*y[0])/(4*self.tau_n*y[1])
        dzz = y[0]/self.tau_n
        # noise sources
        fs = np.random.randn()*np.sqrt(2*dss/self.t_step)
        fz = np.random.randn()*np.sqrt(2*dzz/self.t_step)
        fp = np.random.randn()*np.sqrt(2*dpp/self.t_step)
        fn = fz - fs
        return [fn, fs, fp]

    def carrier_density(self, t, y):        
        return self.get_current(t)/(e*self.v)-y[0]/self.tau_n-self.vg*self.a0*y[1]*(y[0]-self.n_t)/(1+self.epsilon*y[1])
    
    def photon_density(self, y):
        return (self.gamma * self.a0 * self.vg * ((y[0]-self.n_t)/(1+self.epsilon*y[1]))-1/self.tau_p)*y[1]+self.beta*self.gamma*y[0]/self.tau_n
    
    def optical_phase(self, y):
        return self.alpha/2 * (self.gamma*self.vg*self.a0*(y[0]-self.n_t)-1/self.tau_p)

    def get_im_response(self, f, type='exact'):
        Y, Z = self.im_response_yz()
        if type=='exact':
            return Y, Z, self.im_response_hf(f,Y,Z)
        elif type=='approx.':
            return Y, Z, self.im_response_haf(f,Y,Z)
        else:
            print('Invalid type of IM response.')
            return -1

    def im_response_hf(self,f,Y,Z):
        return Z/((1j*2*np.pi*f)**2+1j*2*np.pi*f*Y+Z)

    def im_response_haf(self,f,Y,Z):
        fr = 1/(2*np.pi)*np.sqrt(Z)
        return fr**2/((fr**2-f**2)+1j*f*Y/(2*np.pi))

    def im_response_yz(self):
        # Intensity modulation frequency response - see ref. [1]
        Y = self.vg*self.a0*self.S[-1]/(1+self.epsilon*self.S[-1]) + 1/self.tau_n + 1/self.tau_p
        Y = Y - self.gamma*self.a0*self.vg*(self.N[-1]-self.n_t)/(1+self.epsilon*self.S[-1])**2
        Z = self.vg*self.a0*self.S[-1]/(1+self.epsilon*self.S[-1]) * 1/self.tau_p + 1/(self.tau_p*self.tau_n)
        Z = Z + (self.beta-1)*self.gamma*self.a0*self.vg/self.tau_n*(self.N[-1]-self.n_t)/(1+self.epsilon*self.S[-1])**2
        return Y,Z
    
    def plot(self):
        fig,ax=plt.subplots(2,2,figsize=(12,6))
        ax[0,0].plot(1e9*self.t, 1e3*self.power)
        self._extracted_from_plot_4(ax, 0, 0, 'Optical Power [mW]')
        ax[1,0].plot(1e9*self.t, 1e-9*np.real(self.chirp))
        self._extracted_from_plot_4(ax, 1, 0, 'Chirp [GHz]')
        ax[0,1].plot(1e9*self.t, np.real(self.N))
        self._extracted_from_plot_4(ax, 0, 1, r'Carrier density $N(t)$ [cm$^{-3}$]')
        ax[1,1].plot(1e9*self.t, np.real(self.S))
        self._extracted_from_plot_4(ax, 1, 1, r'Photon Density $S(t)$ [cm$^{-3}$]')
        plt.tight_layout()
        return ax

    def _extracted_from_plot_4(self, ax, arg1, arg2, arg3):
        ax[arg1, arg2].set_xlabel('Time [ns]')
        ax[arg1, arg2].set_ylabel(arg3)
        ax[arg1, arg2].set_xlim(1e9*np.array([self.t.min(),self.t.max()]))
        ax[arg1, arg2].grid(True)
        
if __name__ == "__main__":
    from commpy.utilities  import upsample
    from optic.metrics import signal_power
    from optic.dsp import firFilter, pulseShape
    from optic.core import parameters

    import matplotlib.pyplot as plt
    import random

    from Functions.laser_dfb import laser_dfb

    # simulation parameters
    SpS = 64
    Rs = 2.5e9 # Symbol rate (for OOK case Rs = Rb)
    Tsymb = 1/Rs # Symbol period in seconds
    Fs = 1/(Tsymb/SpS) # Signal sampling frequency (samples/second)
    Ts = 1/Fs # Sampling period

    laser_dfb_parameters = parameters()
    laser = laser_dfb(laser_dfb_parameters)

    # generate pseudo-random bit sequence
    bitsTx = [0, 1, 0, 1, 0] #np.random.randint(2, size=10)
    #n = np.arange(0, bitsTx.size)

    # upsampling
    symbolsUp = upsample(bitsTx, SpS)

    # typical NRZ pulse
    pulse = pulseShape('nrz', SpS)
    pulse = pulse/max(abs(pulse))

    # pulse formatting
    sigTx = firFilter(pulse, symbolsUp)

    t = np.arange(0,SpS*len(bitsTx))*Ts
    plt.plot(1e9*t, sigTx, label = 'RF binary signal', linewidth=2)
    plt.plot(1e9*t[::SpS], bitsTx, 'ro')
    plt.ylabel('Amplitude (a.u.)')
    plt.xlabel('Time (ns)')
    plt.xlim([1e9*t.min(), 1e9*t.max()])
    plt.legend(loc='upper left')
    plt.grid()

    sol = laser.solve(t, sigTx)
    laser.plot()
    plt.figure()
    plt.plot(laser.t, laser.current)