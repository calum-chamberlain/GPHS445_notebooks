# VUW GPHS445: Observational Earthquake Seismology
## Jupyter notebooks for the GPHS445 course taught at Victoria University of Wellington

[![Python 3.8](https://img.shields.io/badge/python-3.8-blue.svg)](https://www.python.org/downloads/release/python-380/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![test](https://github.com/calum-chamberlain/GPHS445_notebooks/workflows/test/badge.svg)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/calum-chamberlain/GPHS445_notebooks/master)

This repository contains relevant notebooks and links to notebooks to help students
through practical elements of observational seismology. This is **NOT** an exhaustive
textbook for observational seismology.

These notebooks cover:
1. [Introduction to handling seismic data with ObsPy](1_Intro_to_processing.ipynb)
2. [Fundamentals of Fourier analysis for Seismic data](2_Fourier_analysis.ipynb)
3. [Earthquake detection and seismic phase identification](3_Earthquake_detection_and_phase_analysis.ipynb)
4. [Earthquake location](4_Earthquake_Location.ipynb)
5. [Magnitude and focal mechanisms](5_Magnitudes_and_focal_mechanisms.ipynb)
6. [Earthquake statistics](6_Earthquake_statistics.ipynb)

    
## How to use these notebooks
1. If using github:
    - Fork this repository (forking will allow you to keep your local changes seperate from changes here);
    - Clone (`git clone ...` where ... is the url given by the *big green button*) to your local machine
1. If just downloading:
    - Download using the *big green button*
2. Change into the newly created directory;
3. Install the requirements (recommended to use [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#id2)):
    ```bash
    conda env create -f environment.yml  # Create an environment called gphs445
    source activate gphs445  # Activate the environment - You will need to do this everytime you use the notebooks
    ```
4. Start jupyter-lab (run `jupyter-lab`) in your repository directory and navigate to the notebook you 
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
