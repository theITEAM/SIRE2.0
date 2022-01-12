# Susceptibility Infectivity and Recoverability Estimation (SIRE 2.0)

## Introduction

In the era of rapid expansion of the human population with increasing demands on food security, effective solutions that reduce the incidence and impact of infectious diseases in plants and livestock are urgently needed. Even within a species hosts differ widely in their response to infection and therefore also in their relative contribution to the spread of infection within and across populations. Three key epidemiological host traits affect infectious disease spread: susceptibility (propensity to acquire infection), infectivity (propensity to pass on infection to others) and recoverability (propensity to recover quickly). Disease control strategies aimed at reducing disease spread may, in principle, target improvement in any one of these three traits.

SIRE 2.0 allows for simultaneous estimation of individual additive genetic effects (along with the corresponding covariance matrix), single nucleotide polymorphism (SNP) and treatment effects on these host traits (so identifying potential pleiotropic effects). SIRE implements a Bayesian algorithm which makes use of temporal data (consisting of any combination of recorded infection times, recovery times or disease status measurements) from multiple epidemics whose dynamics can be represented by the susceptible-infectious-recovered (SIR) or susceptible-infectious (SI) model. 


## Download

The following files can be downloaded:

* Windows: [SIREv2.0.zip](https://github.com/theITEAM/SIRE2.0/releases/download/v2.0/SIREv2.0.zip)

* Linux: Not Currently Available

* Mac: Not Currently Available

## Documentation

Information about the software can be obtained from the [SIRE 2.0 manual](https://github.com/theITEAM/SIRE2.0/raw/master/SIRE%202.0%20Manual.pdf) or a [paper on BioRxiv](https://www.biorxiv.org/content/10.1101/618363v3.full))

## Build

To edit and rebuild this software the following instructions must be followed:

* The files from this repository are first downloaded onto your own computer.

* The C++ code in the "Execute" directory must be compiled on your platform of choice (Windows / Linux / Mac). For example "g++ sire.cc tinyxml2.cc -o a.exe -O3" can be used. The resulting "a.exe" executable file is placed into the "Execute" directory.

* This software relies of NW.js to run the graphical user interface. This can be downloaded [here](https://github.com/nwjs/nw.js) for your platform of choice.  

* All the files and folders downloaded from this repository are copies directly into the NW folder. 

* SIRE is run by clicking on "NW.exe".
