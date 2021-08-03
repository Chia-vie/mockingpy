# mockingpy
Create galaxy mock spectra (and more) using Python.

This package is currently under development and will be updated continuously. 

Use at your own risk (see LICENCE for details). 

## How to use mockingpy:

Clone this repository 

Download MILES data (see below): 

Install the package using `python setup.py install` 

In python `import mockingpy` 

This will create a `config.yml` file in your current working directory. 

Edit this file according to your needs and specify the paths to your MILES and particle data (see below). 

## MILES data:

### 1) MILES SSP models   
Download the MILES model directories e.g. MILES_BASTI_CH_baseFe/*fits from here:

http://research.iac.es/proyecto/miles/pages/webtools/tune-ssp-models.php

You might need to change your browser in case you face any difficulties.

Specify the path to all your models e.g. milesdata/models/ in the config.yml file. 

### 2) MILES mass to light ratios
Download the M/L files ('out_phot_*.txt') from here (one for each isochrone and IMF-type): 

http://research.iac.es/proyecto/miles/pages/predicted-masses-and-photometric-observables-based-on-photometric-libraries.php

ATTENTION: Missing headers. Put this line on the very top of each file:

`IMF_type slope [M/H] Age U B V R I J H K U-V B-V V-R V-I V-J V-H V-K (M/L)U (M/L)B (M/L)V (M/L)R (M/L)I (M/L)J (M/L)H (M/L)K F439W F555W F675W F814W F439W-F555W F555W-F675W F555W-F814W`

## Particle data:
Please provide a header in your particle data files. 

Edit the column names of your the age, metalicity, and mass column accordingly in the config.yml file.


