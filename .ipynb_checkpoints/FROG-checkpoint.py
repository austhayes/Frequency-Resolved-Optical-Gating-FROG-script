""" 
Welcome to the FROG package! Here I will display quite a bit of my code utilized throughout my project. 

Imports below:

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import plotly.express as px
import scipy.constants as const
import StackGP as sgp
from StackGP.Wrapper import StackGPRegressor


"""

For starters we will need to load in the files. The 3 files one should take from the spectrometer in their frequency resolved optical gating (FROG) setup, should consist of a data file (intensities), a wavelengths file, and a time file. Utilizing these, I will load in said files via a function, called loading_files.

"""

def load_files(filename = 'data_1400_afterBS.dat', wave = 'wavelengths_1400_afterBS.dat', time = 'times_1400_afterBS.dat'):
    spec = pd.read_csv(filename, header = None, delimiter =',')

    lam = pd.read_csv(wave, header = None, delimiter =',')

    time = pd.read_csv(time, header = None, delimiter =',')

    return spec, lam, time


""" 

The times correspond to the column labels, the rows correspond to the wavelengths. So below the cleanup function converts the inputted time and lam files into lists, then sets them as the row and column labels respectively.

"""
def cleanup(spec, lam, time):
    lam = lam[0].to_list() 
    time = time[0].to_list()
    
    fixed_df = pd.DataFrame(data = spec.values, index = lam, columns = time)

    return fixed_df, lam, time

"""

Heat map plot, showcases frequency (over wavelength) vs. time, with the colors correspconding to intensity.

"""

def heat_map(df, y, x, frequency_units, time_units, x_boundary, y_boundary, colormap, title_text): #df = fixed_df, y = wave_clean, x = times_clean, all from cleanup. 
    times, wavelength = np.meshgrid(x, y)
    freq = []
    for i in range(len(y)):
        freq.append(const.c/y[i]) #Converting from wavelength to frequency.

    plt.pcolormesh(times, freq, df, shading = 'nearest', cmap=f'{colormap}')
    plt.colorbar(label='Intensity')
    plt.xlim(x_boundary[0], x_boundary[1])
    plt.ylim(y_boundary[0], y_boundary[1])
    plt.ylabel(f'Frequency ({frequency_units})')
    plt.xlabel(f'Time ({time_units})')
    plt.title(f'{title_text}')
    plt.show()

"""
Utilizing StackGP in order to fit the data, this will return the equation that best fits the data, we will then begin error checking after working with said outputted equation.

"""
"""
def fitted_data(x, y, z, gens):
    inputted_data = np.concatenate([(x.to_numpy(), y.to_numpy(), z.to_numpy())])
    inputted_data = np.array(inputted_data)
    response = inputted_data[0]/(1.3-inputted_data[1]/inputted_data[2])
    models=sgp.evolve(inputted_data,response,generations=gens,tracking=True,ops=sgp.allOps())
    return sgp.printGPModel(models[0]), models[0]
"""
"""
def fitted_data(x, y, z, gens = 100):
    X, Y = np.meshgrid(x, y)
    X_ = np.vstack([X.ravel(), Y.ravel()])
    model = StackGPRegressor(generations=gens, ops="allOps" )
    model.fit(X_, z)
    return model.get_model_string()
"""
"""
Analysis

"""
def Analysis():
    print("This function will output the analysis of the FROG data")
    """
    Mathematically, the first step, Fourier transform to convert data from the frequency domain to the time domain, and vice versa=. This is a critical next step that will be relied on heavily throughout the rest of the project. Other than this, there is a lot of sort of splitting up of the data that needs to be done. The matlab code does a lot of extra work with the data files that I am not positive is really all that necesssary so I have been spending a bit of time trying to figure out what I need, and how to do this in a clear and concise manner to make sure it is easy to understand, as it is currently, not at all. 

    """

    
    

