![Logo](https://github.com/Chia-vie/mockingpy/blob/main/logo.png) 
# mockingpy 

Create galaxy mock spectra (and more) using Python.

This package is currently under development and will be updated continuously. 

Use at your own risk.

## Attribution

Written by Christine Ackerl, based on work by Alina Boecker

## How to use mockingpy:

1) Clone this repository into a local working directory

2) Install the package using `python setup.py install` in the same directory

3) Download MILES data (see below)

4) In python `import mockingpy` 

   This will create a `config.yml` file in your working directory. 

5) Edit this config file according to your needs and specify the paths to your MILES and particle data (see below). 

6) To create the mock spectra run `mockingpy.MakeMock('config.yml')` or change the config filename accordingly

The resulting spectra will be stored in `.fits` files.

## MILES data:

### 1) MILES SSP models   
Download the MILES model directories e.g. MILES_BASTI_CH_baseFe/*fits from here:

http://research.iac.es/proyecto/miles/pages/webtools/tune-ssp-models.php

You might need to change your browser in case you face any difficulties.

Specify the path to all your models (e.g. milesdata/models/) in the config file. 

### 2) MILES mass to light ratios
Download the M/L files ('out_phot_*.txt') from here (one for each isochrone and IMF-type): 

http://research.iac.es/proyecto/miles/pages/predicted-masses-and-photometric-observables-based-on-photometric-libraries.php

ATTENTION: Missing headers. Put this line on the very top of each file:

`IMF_type slope [M/H] Age U B V R I J H K U-V B-V V-R V-I V-J V-H V-K (M/L)U (M/L)B (M/L)V (M/L)R (M/L)I (M/L)J (M/L)H (M/L)K F439W F555W F675W F814W F439W-F555W F555W-F675W F555W-F814W`

Specify the path to all M/L files (e.g. milesdata/masstolight/) in the config file.

### 3) Flux tables
The code will create intermediate output files for you, such as flux tables. 

Please specify the path where those files should be stored (e.g. milesdata/output/) in the config file. 

## Particle data:
Please provide a header in your particle data files. 

Specify the names of the age, metallicity, and mass columns in the config file.

Specify the path to your particle data (e.g. particles/simulation1/) in the config file. 

