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

    def projection(self, intensity, freq, threshold_RMS=0.01, max_iterations=500,):

        """

    The Projection algorithm is very mathematical, the goal here is to find the correct electric field equation that describes the radiation being      observed within the system. I utilized a variety of resources to help critique the function, google AI was utilized to help with the projection     algorithm. URLS: https://link.springer.com/article/10.1186/s43020-024-00139-4, https://arxiv.org/html/1708.03560v1.

        """
        N = len(intensity.columns)

        t = np.linspace(-7, 7, N)

        E_t = np.sqrt(2*math.log(2))*np.exp(-t**2) * np.exp(1j * np.random.uniform(0, 2*np.pi, N))
        
        rms_errors = []
        print(max_iterations)

        for i in range(0, max_iterations):
            # Transform to the frequency domain, utilizing a fast fourier transform.
            E_f = np.fft.fft(E_t)
            I_calc = np.abs(E_f) ** 2
            
            # Root mean squared error between calculated and measured intensity.
            rms_error = np.sqrt(np.mean((I_calc - intensity) ** 2))
            rms_errors.append(rms_error)
            
            # Convergence condition
            if rms_error < threshold_RMS:
                print(f"Converged at iteration {i} with RMS error: {rms_error:.4f}")
                break
                
            # Applying projection
            calculated_phase = np.angle(E_f)
            E_f_updated = np.sqrt(intensity) * np.exp(1j * calculated_phase)
            
            # Utilizing an inverse fourier transform to transform back into the time domain. 
            E_t_updated = np.fft.ifft(E_f_updated)
            
            # Projection with a mask related to assumed pulse length to filter out noise. 
            E_t = E_t_updated  

            #Need to find FWHM to find time:
            maximum_val = intensity.max()
            half_max = maximum_val/2
                    
            I_t = np.abs(E_t) ** 2

            # Finding main peak
            #peaks, _ = find_peaks(I_t, height=np.max(I_t) * 0.5)
            
            # Calculating FWHM, 0.5 corresponds to half of the peak height
            #results_half = peak_widths(I_t, peaks, rel_height=0.5)
            
            # Convert width from index count to units of time
           ## dt = t[1] - t[0]
            ##fwhm_samples = results_half[0][0]
            ##fwhm_time = fwhm_samples * dt
                
            ##time_window_mask = (fwhm_time >= -100) & (fwhm_time <= 100)
            ##E_t[~time_window_mask] = 0  # Force everything outside the window to be 0
        else:
            print(f'Maxed iterations reached. Max iterations: {max_iterations}')

        
        E_t = new_arr = np.nan_to_num(E_t)
            
        return E_t, rms_errors
        

"""
Analysis

"""
def Analysis():
    print("This function will output the analysis of the FROG data")
    """
    Mathematically, the first step, Fourier transform to convert data from the frequency domain to the time domain, and vice versa=. This is a critical next step that will be relied on heavily throughout the rest of the project. Other than this, there is a lot of sort of splitting up of the data that needs to be done. The matlab code does a lot of extra work with the data files that I am not positive is really all that necesssary so I have been spending a bit of time trying to figure out what I need, and how to do this in a clear and concise manner to make sure it is easy to understand, as it is currently, not at all. 

    """
