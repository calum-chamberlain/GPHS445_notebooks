# VUW GPHS445: Observational Earthquake Seismology
## Jupyter notebooks for the GPHS445 course taught at Victoria University of Wellington

[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This repository contains relevant notebooks and links to notebooks to help students
through practical elements of observational seismology.

The course covers:
1. Recording and processing digital seismic data
    - How seismographs work
    - How digitizers work
    - Fourier transforms
    - Filtering and resampling
    - Introduction to processing seismic data with [Obspy](http://docs.obspy.org/)
2. Introduction to seismic waves
    - Main seismic wave types and fundamental properties
    - Reflection, refraction and conversion of phases
    - Dispersion and attenuation
3. Building an earthquake catalog
    - Earthquake detection
    - Phase picking and the impact of filters
    - Earthquake location
    - Magnitude calculation
    - Focal mechanism determination
    - Source-time functions
4. Seismic hazard
    - Magnitude-frequency relations
    - Earthquake clustering
    - Aftershock forecasting
    
    
## How to use these notebooks
1. Fork this repository
2. Install the requirements (recommended to use [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#id2)):
    ```bash
    conda create -n gphs445  # Create an environment called gphs445
    source activate gphs445  # Activate the environment - You will need to do this everytime you use the notebooks
    # Install the requirements
    while read requirement; do conda install --yes $requirement; done < requirements.txt
    ```
3. Start jupyter in your repository directory and navigate to the notebook you 
   want to work on.
4. Save any changes you make in the notebook.
5. When you want to back-up your changes, or when you are happy with them, commit the
   changes and push to your upstream repository 
   (check out the [github cheatsheet](https://services.github.com/on-demand/downloads/github-git-cheat-sheet.pdf) for more commands):
   ```bash
   git add <notebook you have been working on> # Replace with the filename you were working on
   git commit -m "Some memorable commit message"
   git push origin master
   ```