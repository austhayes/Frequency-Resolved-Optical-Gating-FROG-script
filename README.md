# FROG README.md
# By Austin Hayes

# Summary

The main goal of this project is to rewrite the current Frequnecy Resolved Optical Gating (FROG) MATLAB script(s) in Python. In doing so, this will allow my lab much more usability, and the potential for much more usability of the files obtained from the spectrometer in the FROG setup. A chirped light pulse is a pulse that makes certain colors travel faster than others, ie. light of different frequencies arrive at different times relative to some detection device. The importance of this is after sending a collimated beam through a nonlinear crystal relative to the beams focus after traveling through a lens, we can chirp the light pulse(s), which spreads the frequencies out, ie. increases the bandwidth. Utilizing the FROG setup, we can then detect the bandwidth relative to the original by splitting the pulse into two pulses, and then utilizing a stage with mirrors on it to change the time delay between the two pulses before the head through another crystal. From here, we determine how much the beams bandwidth increased, which in turn lets us know the minimum length in time this light pulse can be. Overall, the FROG setup is a device that measures how much we've increased the bandwidth which then tells us how short we can make laser pulses in time. In Terahertz-Scanning Tunneling Microscopy (THz-STM), we rely on short pulses of light (In the THz frequency regime) to drive electrons across a tunnel junction (For us that includes the tungsten (W) STM tip (The probe), the gap separating the tip and the sample, and the sample itself). The shorter we can make these pulses, the better we can resolve the atomic structure of these samples in time. Overall, after the beam goes through the FROG setup and hits the spectrometer, we are given a few .M files with the main input being intensity vs. delay time and wavelength of the beam. From here we get some very important results, with the main two being intensity vs. time and wavelength, as well as the relative phase vs. wavelength which can tell us a ton about our limitations when it comes to the individual pulses. 

The inputs expected for this project are the intensity vs. wavelength and delay time. I will include functions that can do conversions to the frequnecy domain (useful for obtaining the E field and time), as well as converting step size (dx) and number of steps (N) to delay times as well. The overall workflow doesn't have a ton of parts but there are 4 scripts that need to be rewritten: NewloadFrogTrace.m, FROGAnalysis.m, FROGAnalyzer.m, Quick_Graphing.m. The workflow will follow the route of: reading in the files into my version of NewloadFrogTrace.m (1_NewLoadFrogTrace), from here I will work on extracting the data from said files and utilizing them within my versions of FROGAnalysis.m and FROGAnalyzer.m (2_FROGAnalysis) This is where majority of my time will be spent because I will need to thorughly go through interpolations, fourier transforms, and other mathematical (as well as estimations) that will result in accurate results (Majority of this work done within FROG.py). From there I will just need to write a plotting script from the outcomes as well as showcase given import values in optics such as the group delay dispersion (GDD) (All results will be within 3_FROGResults).

# Workflow

First two weeks: Continue working on understanding the matlab script, read in the files from the spectrometer via python, ideally begin on rewriting the analysis files. 

Essentially the rest of the project: Work through and organize everything in the analysis python code but make sure all of the equations are worked out correctly so that there isn't a domino effect of errors in calculations. Also, cross check the equations with the matlab code, ie. run sections of the matlab code so that every answer in the python code is double checked with the original code. Rather than continuing to try to go block by block mimicking the matlab script, instead I will clean and create a new structure that is more user friendly, as found out by the first few weeks of this project, trying to even start on FROGAnalysis the way the matlab code is not the way to go about this. 

Final week to two weeks: Rewrite the plotting software in python, and continuously test the completed files to ensure they are correct. 

As for the overall workflow of the project, you should start in 1_NewLoadFrogTrace, which will generate your initial wavelength (or frequency) vs. time with the heatmap representing intensity. From there, you should take this data into 2_FROGAnalysis file where you should be able to put in your file(s) and get back, intensity vs. time, phase vs. time, intensity vs. frequency,  phase vs. frequency, FROG Error, Group delay dispersion (GDD), GDD error, and how the modeled result differs from the normal result. Lastly, we will take all results and analyze them in comparison to the initial FROG MATLAB script within 3_FROGResults.

# Points to focus on

Continuously push everything to github, backup files, work in an organized manner, and double check everything with the results we get in matlab. The original software starts with an initial guess, then iterates back and forth between the simulated data and the real data, from here it continues iterating back and forth until error is minimized (denoted by the threshold_RMS variable within the frog_projection_new() function), the path taken throughout this project will be slightly different with the same goals in mind, create a program that minimizes the error within the simulated waveform, ideally improving not only usability of the code, but lessening the overall error of the results. This will be cross checked throughout the project by double checking results against the matlab program. 

# Getting started with FROG software

Within the file: 1_NewLoadFrogTrace.ipynb, the main focus will be to access the environment and then load in the main files obtained from a spectrometer, those files include: intensities, wavelengths, and times. We will then plot out the newly loaded in data after constructing a proper data frame in a few steps. Follow the steps below:

1. We will begin by making sure users are within the environment: FROG_ENV. Access the environment by running the following code in the command line:

     A. make init
     
     
     B. conda activate FROG_ENV

 
     C. conda install ipykernel
     

     D. python -m ipykernel install --user --name=FROG_ENV --display-name "Python (FROG_ENV)"
     

     E. If you use jupyer lab write: jupyter lab, if you use jupyter notebook write: jupyter notebook

2. Declare an instance of the class "loading_mapping()" by setting a variable equal to frog.loading_mapping()

3. Load in the files utilizing the load_files() function from the FROG.py package.

4. Cleanup the dataframe using the cleanup() function and setting the time to be the x axis labels and the wavelengths to be the y axis labels within said function.

5. The frequency of the pulses is a very important variable to keep track of, we do the conversion twice, once in the freq() function, and another time within heat_map function. The code may be improved to just utilize the other function but for simplicity, the heat_map function will continue doing the conversion as well from the wavelengths file.

6. Store the new_df, new_wave, new_time, and freq data to be utilized in later notebooks.

After this step will be utilizing this newly formed dataset to create a heatmap of the wavelength vs. time with the intensity showcased as a value along a colorbar.

7. Plot the heatmap utilizing the heat_map() function, pay close attention to the inputs.

From here, the result should showcase a heat map like that of which is seen on your spectrometer.

8. Once confirmed that the heat map above matches that seen on the users spectrometer software:
    Move onto notebook 2, file: 2_FROGAnalysis.ipynb
   