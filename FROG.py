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

class obtain_time():
    def time(self, freq):
        """

The goal within this function is to obtain the measurement time to later be used within the electric field calculations. The current time we have that is outputted from the spectrometer is the time delay between the two pulses in reference to the crystal.

        """
        t = np.fft.ifft(freq)
        return t

"""

The E field class will work on computing the electric field, via a variety of ways. The first, we will utilize a less specific approach where we simply find the magnitude of the electric field for the radiation (light pulse). Utilizing the equation: I(t) = 1/2*c*n*eps_0*|E(t)|^2, we can solve for the (modulus) of E(t) by solving for E(t) => = +-|E(t)| = (2*I(t)/(c*n*eps_0))^(1/2).

Variables:

    c = speed of light (299,792,458 m/s)
    eps_0 = permitivity of free space (8.854*10^(-12) F/m (Farads/meter)
    n = Refractive index of medium

"""

class E_field():
    
    def get_E_mod_df(self, df, n):
        E_t = np.sqrt(2*df/(const.c * n * const.epsilon_0))
        return E_t

    def get_E_sig_SHG(self, df, n, tau, delta_x = None):   
        """

        Since second harmonic generation (SHG, aka frequency double) is of key importance in our lab, we will be focusing on E fields of said form,         E_signal(t, tau) = E(t)*E(t-tau). Individuals can calculate the time delay (tau) from the difference in path lengths the two beams within           the FROG take if needed. Users should input either a delta_x or a tau value. 
        
        """

        if delta_x != None:
            tau = []
            for i in range(len(delta_x)):
                tau.append(2*delta_x[i]/const.c)
        else:
            tau = tau

        
        E_t_approx = []
        E_t_tau_approx = []
        
        E_sig = []
        
        for i in range(len(t)):
            E_field_guess_t = np.sqrt(2*math.log(2)*np.exp(-(t[i])**2))
            E_t_approx.append(E_field_guess_t)
            
            E_field_guess_t_tau = np.sqrt(2*math.log(2)*np.exp(-(t[i] - tau[i])**2))
            E_t_tau_approx.append(E_field_guess_t_tau)

            E_signal = (E_field_guess_t * E_field_guess_t_tau)
            E_sig.append(E_signal)
        
        
        return  E_t_approx, E_t_tau_approx, E_sig


"""
Analysis

"""
def Analysis():
    print("This function will output the analysis of the FROG data")
    """
    Mathematically, the first step, Fourier transform to convert data from the frequency domain to the time domain, and vice versa=. This is a critical next step that will be relied on heavily throughout the rest of the project. Other than this, there is a lot of sort of splitting up of the data that needs to be done. The matlab code does a lot of extra work with the data files that I am not positive is really all that necesssary so I have been spending a bit of time trying to figure out what I need, and how to do this in a clear and concise manner to make sure it is easy to understand, as it is currently, not at all. 

    """

"""
SCRATCH WORK:

Attempt at fitting the data:


def fitted_data(x, y, z, gens):
    inputted_data = np.concatenate([(x.to_numpy(), y.to_numpy(), z.to_numpy())])
    inputted_data = np.array(inputted_data)
    response = inputted_data[0]/(1.3-inputted_data[1]/inputted_data[2])
    models=sgp.evolve(inputted_data,response,generations=gens,tracking=True,ops=sgp.allOps())
    return sgp.printGPModel(models[0]), models[0]
"""

    
    

