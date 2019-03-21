# VUW GPHS445: Observational Earthquake Seismology
## Jupyter notebooks for the GPHS445 course taught at Victoria University of Wellington

[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CircleCI](https://circleci.com/gh/calum-chamberlain/GPHS445_notebooks.svg?style=svg&circle-token=1a58ff80d0b826d24bf8b733e205775e545d9bee)](https://circleci.com/gh/calum-chamberlain/GPHS445_notebooks)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/calum-chamberlain/GPHS445_notebooks/master)

This is the devlopment branch!

This repository contains relevant notebooks and links to notebooks to help students
through practical elements of observational seismology.

The course covers:
1. [Recording and processing digital seismic data](1_Processing_and_recording/README.md)
    - [Introduction to Python](1_Processing_and_recording/0_Python.ipynb)
    - [How seismographs work](1_Processing_and_recording/1_Seismographs.ipynb)
    - [How digitizers work](1_Processing_and_recording/2_Digitizers.ipynb)
    - [Fourier transforms](1_Processing_and_recording/3_Fourier_Transforms.ipynb)
    - [Filtering and resampling](1_Processing_and_recording/4_Filtering_Resampling.ipynb)
    - [Introduction to Obspy](1_Processing_and_recording/5_Intro_To_Obspy.ipynb) seismic data with [Obspy](http://docs.obspy.org/)
2. [Introduction to seismic waves](2_Seismic_waves/README.md)
    - [Main seismic wave types and fundamental properties](2_Seismic_waves/1_Seismic_Waves.ipynb)
    - [Reflection, refraction and conversion of phases](2_Seismic_waves/2_Transmission.ipynb)
    - [Dispersion and attenuation](2_Seismic_waves/3_Dispersion_and_Attenuation.ipynb)
3. [Building an earthquake catalog](3_Building_a_catalog/README.md)
    - [Earthquake detection](3_Building_a_catalog/1_Earthquake_Detection.ipynb)
    - [Phase picking and the impact of filters](3_Building_a_catalog/2_Phase_Picking.ipynb)
    - [Earthquake location](3_Building_a_catalog/3_Earthquake_Location.ipynb)
    - [Magnitude calculation](3_Building_a_catalog/4_Magnitudes.ipynb)
    - [Focal mechanism determination](3_Building_a_catalog/5_Focal_Mechanisms.ipynb)
    - [Source-time functions](3_Building_a_catalog/6_Source_Time_Functions.ipynb)
4. [Seismic hazard](4_Hazard/README.md)
    - [Magnitude-frequency relations](4_Hazard/1_Magnitude_Frequency.ipynb)
    - [Earthquake clustering](4_Hazard/2_Earthquake_Clustering.ipynb)
    - [Aftershock forecasting](4_Hazard/3_Aftershock_Forecasting.ipynb)
    
    
## How to use these notebooks
1. If using github:
    - Fork this repository;
    - Clone (`git clone ...` where ... is the url given by the *big green button*) to your local machine
1. If just downloading:
    - Download using the *big green button*
2. Change into the newly created directory;
3. Install the requirements (recommended to use [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#id2)):
    ```bash
    conda env create -f environment.yml  # Create an environment called gphs445
    source activate gphs445  # Activate the environment - You will need to do this everytime you use the notebooks
    ```
4. Start jupyter (run `jupyter notebook`) in your repository directory and navigate to the notebook you 
   want to work on;
5. Save any changes you make in the notebook;
6. When you want to back-up your changes, or when you are happy with them, commit the
   changes and push to your upstream repository 
   (check out the [github cheatsheet](https://services.github.com/on-demand/downloads/github-git-cheat-sheet.pdf) for more commands):
   ```bash
   git add <notebook you have been working on> # Replace with the filename you were working on
   git commit -m "Some memorable commit message"
   git push origin master
   ```
