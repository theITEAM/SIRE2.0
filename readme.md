# Susceptibility Infectivity and Recoverability Estimation (SIRE 2.0)

## Introduction

In the era of rapid expansion of the human population with increasing demands on food security, effective solutions that reduce the incidence and impact of infectious diseases in plants and livestock are urgently needed. Even within a species hosts differ widely in their response to infection and therefore also in their relative contribution to the spread of infection within and across populations. Three key epidemiological host traits affect infectious disease spread: susceptibility (propensity to acquire infection), infectivity (propensity to pass on infection to others) and recoverability (propensity to recover quickly). Disease control strategies aimed at reducing disease spread may, in principle, target improvement in any one of these three traits.

SIRE 2.0 allows for simultaneous estimation of individual additive genetic effects (along with the corresponding covariance matrix), single nucleotide polymorphism (SNP) and treatment effects on these host traits (so identifying potential pleiotropic effects). SIRE implements a Bayesian algorithm which makes use of temporal data from multiple contact groups whose dynamics can be represented by the susceptible-infectious-recovered (SIR) or susceptible-infectious (SI) model. 

SIRE 2.0 uses an efficient Markov chain Monte Carlo (MCMC) algorithm to draw samples from the posterior destribution. A graphical interface takes as input any combination of information about individuals’ infection and recovery times, disease status measurements and disease diagnostic test results. The genomic/pedigree relationship matrix can be specified (if included in the model), along with the genotypes of SNPs or other fixed effects (if included), as well as details of which individuals belong in which contact group. A range of prior specifications can be made on different model parameters. The outputs consists of posterior trace plots for model parameters, distributions, visualisation of infection and recovery times, dynamic population estimates and summary statistics (means and credible intervals) as well as MCMC diagnostic statistics. Files containing posterior samples of parameters and events can be exported for further analysis using other tools. 

## Download Binaries

The following files can be downloaded:

* Windows: Download the file [SIRE_v2.0_windows.zip](https://github.com/theITEAM/SIRE2.0/releases/download/v2.02/SIRE.v2.0_windows.zip) and unzip. SIRE is run by clicking on the “SIRE.exe” icon.

* Mac: Download the file [SIRE_v2.0_Mac.zip](https://github.com/theITEAM/SIRE2.0/releases/download/v2.02/SIRE_v2.0_Mac.zip). SIRE is run by clicking on the “SIRE.app” icon (if the error message “SIRE can’t be opened because it is from an unidentified developer…” appears, right clicking on “SIRE.app” and selecting “Open” will allow the option to run).

* Linux: Download the file    [SIRE_v2.0_linux.tar.gz](https://github.com/theITEAM/SIRE2.0/releases/download/v2.01/SIRE_v2.0_linux.tar.gz). This can then be extracted by using the terminal command “tar -zxvf SIRE_v2.0_linux.tar.gz”. The code is executed using “./SIRE”. Note, the linux version is slightly older, and some features have not been implemented yet. 

## Future developments

Efforts are currently underway to incorporate the Bayesian methodology presented in this paper into a more widely applicable software tool called  [BICI](https://github.com/theITEAM/BICI) (Bayesian Individual-based Compartmental Inference), which allows for greater flexibility in model definition as well as supporting a richer array of data types. BICI will also include an R-based interface for various (e.g.  ASReml) styles of scripts for data input and model specification.

## Documentation

Information about the software can be obtained from the [SIRE 2.0 manual](https://github.com/theITEAM/SIRE2.0/raw/master/SIRE%202.0%20Manual.pdf) or a [paper on BioRxiv](https://www.biorxiv.org/content/10.1101/618363v3.full)

## Build

To edit and rebuild this software the following instructions must be followed:

* The files from this repository are first downloaded onto your own computer.

* The C++ code in the "Execute" directory must be compiled on your platform of choice (Windows / Linux / Mac). For example "g++ sire.cc tinyxml2.cc -o a.exe -O3" can be used. The resulting "a.exe" executable file is placed into the "Execute" directory.

* This software relies of NW.js to run the graphical user interface. This can be downloaded [here](https://github.com/nwjs/nw.js) for your platform of choice.  

* All the files and folders downloaded from this repository are copied directly into the NW folder. 

* SIRE is run by clicking on "NW.exe".
