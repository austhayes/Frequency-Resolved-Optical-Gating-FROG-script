""" 
Welcome to the FROG package! Here I will display quite a bit of my code utilized throughout my project. 

Imports below:

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import scipy.constants as const
import StackGP as sgp
from StackGP.Wrapper import StackGPRegressor
from scipy.signal import find_peaks, peak_widths

"""
Beginning by a defining a class where the functions load in the spectrometer files (function: load_files()), organize the dataframe correctly (function: cleanup()), then plot the heatmap that should be the same as that seen in there spectrometer software (function: heat_map). NOTE: Pay close attention to the inputs! The inputs get slightly dense especially for the heat_map() function, so it is best to review inputs thoroughly.

"""

class loading_mapping():
    """
    
    For starters we will need to load in the files. The 3 files one should take from the spectrometer in their frequency resolved optical gating (FROG) setup, should consist of a data file (intensities), a wavelengths file, and a time file. Utilizing these, I will load in said files via a function, called loading_files.
    
    """
    
    def load_files(self, filename = 'data_1400_afterBS.dat', wave = 'wavelengths_1400_afterBS.dat', time = 'times_1400_afterBS.dat'):
        spec = pd.read_csv(filename, header = None, delimiter =',')
    
        lam = pd.read_csv(wave, header = None, delimiter =',')
    
        time = pd.read_csv(time, header = None, delimiter =',')
    
        return spec, lam, time
    
    
    """ 
    
    The times correspond to the column labels, the rows correspond to the wavelengths. So below the cleanup function converts the inputted time and lam files into lists, then sets them as the row and column labels respectively.
    
    """
    def cleanup(self, spec, lam, time):
        lam = lam[0].to_list() 
        time = time[0].to_list()
        
        fixed_df = pd.DataFrame(data = spec.values, index = lam, columns = time)
    
        return fixed_df, lam, time
    
    def heat_map(self, df, y, x, frequency_units, time_units, x_boundary, y_boundary, colormap, title_text): 
        #Example inputs: heat_map(df = fixed_df, y = wave_clean, x = times_clean,.. etc.) all obtained from cleanup function. 
    
        """
        
        Heat map plot, showcases frequency (over wavelength) vs. time, with the colors correspconding to intensity.
        
        """
        
        times, wavelength = np.meshgrid(x, y)
        freq = []
        for i in range(len(y)):
            freq.append((const.c *10**9)/(y[i]*10**12)) #Converting from wavelength (nm) to frequency (nm to THz).

        #Normalizing the intensity data utilizing the min-max method to make the maximum intensity 1 and the minimum 0.

        df_  = (df - df.min()) / (df.max() - df.min())
    
        plt.pcolormesh(times, freq, df_, shading = 'nearest', cmap=f'{colormap}')
        plt.colorbar(label='Intensity')
        plt.xlim(x_boundary[0], x_boundary[1])
        plt.ylim(y_boundary[0], y_boundary[1])
        plt.ylabel(f'Frequency ({frequency_units})')
        plt.xlabel(f'Time Delay ({time_units})')
        plt.title(f'{title_text}')
        plt.show()



    def freq(self, wave):
        freq = []
        for i in range(len(wave)):
            freq.append((const.c *10**9)/(wave[i]*10**12)) #Converting from wavelength (nm) to frequency (nm to THz).

        return freq

"""
Utilizing StackGP in order to fit the data, this will return the equation that best fits the data, we will then begin error checking after working with said outputted equation.
    
"""

class fitting():
    def fitted_data(self, x, y, z, gens = 100):
        """
Not correct yet, function will not run properly. Still in progress. 

        """
        X, Y = np.meshgrid(x, y)
        X_ = np.vstack([X.ravel(), Y.ravel()])
        model = StackGPRegressor(generations=gens, ops="allOps" )
        model.fit(X_, z)
        return model.get_model_string()


class obtain_tau():
    def get_tau(self, delta_x):
        """

        If users only have the steps taken by the stage, we can convert that to time delay time using this: "get_tau"    function within the obtain_tau class. This is not mandatory as it is also done in the get_E_sig_SHG function below but if users would like to have the values for the time delay in a list, this is a great function to use. 
        
        """
        tau = []
        for i in range(len(delta_x)):
            tau.append(2*delta_x[i]/const.c)
        return tau

class obtain_field():
    def field(self, freq,intensity):
        """

The goal within this function is to setup the field, the projection algorithm will be cleaned up, as well as the docstrings to utilize the obtain field function and continuously tweak it. 

        """
        N = len(freq)
        dt = 0.1
        t = np.arange(-N/2, N/2)*dt
        
        E_sig_guess = np.sqrt(2*math.log(2)*np.exp(-(t[i])**2))

        E_spectrum_t = np.fft.fftshift(np.fft.fft(E_sig_guess, axis=0), axes=0)
        I_t = np.abs(E_spectrum_t) ** 2

        updated_field = np.sqrt(intensity) * np.exp(1j * np.angle(E_spectrum_t))
        

        
        
        return E_spectrum_t, I_t, updated_field

    def frog_projection_new(self, intensity_2d, max_iterations=1000, threshold_RMS=0.05):

       """

    The Projection algorithm is very mathematical, the goal here is to find the correct electric field equation that describes the radiation being      observed within the system. I utilized a variety of resources to help critique the function, google AI was utilized to help with the projection     algorithm. URLS: https://link.springer.com/article/10.1186/s43020-024-00139-4, https://arxiv.org/html/1708.03560v1.

        """
        
        I_meas = intensity_2d.to_numpy() if hasattr(intensity_2d, 'to_numpy') else intensity_2d
        N_f, N_t = I_meas.shape  # N_f = 2048, N_t = 133
        
        # 1: E_t must match N_f (2048) because it maps directly to your frequencies
        E_t = np.exp(-0.5 * (np.linspace(-3, 3, N_f))**2) * np.exp(1j * np.random.uniform(0, 2*np.pi, N_f))
        
        FROG_mag = np.sqrt(np.maximum(I_meas, 0))
        rms_errors = []
    
        # Map the 133 delay steps to index shifts on the 2048 time grid
        # This centers the shifts around 0 delay
        delay_indices = np.linspace(-N_t // 2, N_t // 2, N_t).astype(int)
    
        for iteration in range(max_iterations):
            E_sig_t = np.zeros((N_f, N_t), dtype=complex)
            
            # 1. Construct 2D Signal Field using the mapped delays
            for j in range(N_t):
                shift_amount = delay_indices[j]
                E_shifted = np.roll(E_t, shift_amount) 
                
                # CRITICAL FIX 2: Both sides are now length 2048
                E_sig_t[:, j] = E_t * E_shifted
            
            # 2. Convert to Frequency Domain along the rows (axis 0)
            E_sig_f = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(E_sig_t, axes=0), axis=0), axes=0)
            
            # 3. Calculate Error Metric
            I_calc = np.abs(E_sig_f)**2
            I_calc_norm = I_calc / (np.max(I_calc) + 1e-12)
            I_meas_norm = I_meas / (np.max(I_meas) + 1e-12)
            rms_error = np.sqrt(np.mean((I_calc_norm - I_meas_norm)**2))
            rms_errors.append(rms_error)
            
            if rms_error < threshold_RMS:
                print(f"Converged at iteration {iteration} with RMS: {rms_error:.4f}")
                break
                
            # 4. Intensity Constraint Substitution
            calculated_phase = np.angle(E_sig_f)
            E_sig_f_updated = FROG_mag * np.exp(1j * calculated_phase)
            
            # 5. Transform back to Time Domain
            E_sig_t_updated = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(E_sig_f_updated, axes=0), axis=0), axes=0)
            
            # 6. Extract New 1D E(t) field 
            # Grabbing the central delay column (zero-delay profile)
            E_t_new = E_sig_t_updated[:, N_t // 2] 
            E_t = E_t_new / (np.max(np.abs(E_t_new)) + 1e-12)

            #7. Extract temporal and spectral phase from the Electric field.

            temporal_phase = np.unwrap(np.angle(E_t))

            #E_f = np.fft(E_t)
            E_f = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(E_t)))
            spectral_phase = np.unwrap(np.angle(E_f))
            
            
        return E_t, rms_errors, I_calc, temporal_phase, spectral_phase

    def calc_time(self, wavelengths, E_t):
        bw = const.c/(np.min(wavelengths)*1e-9) - const.c/(np.max(wavelengths)*1e-9)
        dt = (1/bw) * 1e15 #units of fs
        
        N = len(E_t) #Number of steps, should be obtained from the frog_projection_new function

        time_range = np.linspace( -N//2, N//2 - 1, N) * dt
        I_t = np.abs(E_t)**2
        I_norm = I_t/I_t.max()

        half_max = np.where(I_norm >= 1/2)[0]

        if len(half_max) >= 2:
            start = half_max[0]
            end = half_max[-1]

            fwhm = time_range[end] - time_range[start]
            print('Time step:', dt,'fs')
            print('Pulse duration (fwhm):', fwhm, 'fs')

        else:
            print('Could not calculate time.')
        return fwhm, dt
        
