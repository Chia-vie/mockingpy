'''
This code creates mock spectra for simulated galaxy particle data.
See documentation for required input files and formats.

author: Christine Ackerl, based on work by Alina Boecker
date: July 2021
'''
import yaml
import numpy as np
from scipy import stats
from scipy import interpolate
from astropy.io import fits
import pandas as pd
import pickle
from tqdm import tqdm

print(
'''
                       _     _                             
  _  _  _   ___   __  | |__ (_)  _ _    __ _   _ __   _  _ 
 / '   ' \ / _ \ / _| | / / | | | ' \  / _` | | '_ \ | || |
 |_||_||_| \___/ \__| |_\_\ |_| |_||_| \__, | | .__/  \_, |
                                       |___/  |_|     |__/
     ,~~~~~~                      ,******>
  ,,~*       \                ,***       >
 <       ◘    \       ,,,*****     ,,,,#
  *++,,        \,,,,,/       ,**
       \                     /
       |                    /
       /                   /
       \                  / 
        \                /
         `-.            /
            **\      /**
               \~~~~/
               //,||
            ,,//--*/`
           /--/`  *
''')

class MakeMock():
    def __init__(self, config):
        with open(config, "r") as ymlfile:
            self.config = yaml.load(ymlfile, Loader=yaml.FullLoader)
        self.read_config()
        # loop through all possibilities
        for name in self.name_particles:
            self.read_particles(name)
            age_models, met_models, age_grid, met_grid, age_bin, met_bin = self.set_up_grid()
            for t in self.imf_type:
                for s in self.imf_slope:
                    for a in self.alpha:
                        self.flux_table(t, s, a)
                        self.finer_model_grid(t, s, a, age_grid, met_grid, age_models, met_models)
                        self.make_spectrum(name, t, s, a, age_grid, met_grid, age_bin, met_bin)
                        if self.noise:
                            SNR = np.array([20, 50, 100])
                            self.add_noise(name, self.isochrone, t, s, a, SNR)

    def read_config(self):
        '''
        see config.yaml file for explanation
        '''
        print('Reading config file')
        self.path_particles = self.config['particles']['path']
        self.prefix_particles = self.config['particles']['prefix']
        self.name_particles = self.config['particles']['names']
        self.agecol = self.config['particles']['columns']['age']
        self.metcol = self.config['particles']['columns']['metallicity']
        self.masscol = self.config['particles']['columns']['mass']
        self.path_models = self.config['models']['path_models']
        self.path_mass_to_light = self.config['models']['path_mass_to_light']
        self.path_output = self.config['models']['path_output']
        self.isochrone = self.config['models']['isochrone']
        self.imf_type = self.config['models']['imftype']
        self.imf_slope = self.config['models']['imfslope']
        self.alpha = self.config['models']['alpha']
        self.noise = self.config['noise']

    def read_particles(self,name):
        '''
        get age, metallicities and mass from eagle simulation
        for determining min and max values of age and met
        :param name: Name of particle file
        :return:
        '''
        print(f'Reading particle data: {name}')
        particles = pd.read_table('%s%s%s.dat' % (self.path_particles, self.prefix_particles, name), sep="\s+")
        particles.replace([np.inf, -np.inf], np.nan, inplace=True)
        particles.dropna(subset=[self.agecol, self.metcol, self.masscol], inplace=True)
        self.age_tag = particles[self.agecol].to_numpy()
        self.met_tag = particles[self.metcol].to_numpy()
        self.mass_tag = particles[self.masscol].to_numpy()

    def set_up_grid(self):
        '''
        this sets up the grid and bins for the finer model grid (age
        grid stays the same as models, metallicity grid in increments of 0.01 dex)
        note: data in df_grid corresponds to the actual grid points,
        bins correspond to boundaries (they lie exactly between
        the actual grid points);
        particles data that falls into a bin will be assigned the
        corresponding grid point
        '''
        print('Setting up age-metallicity grid')
        # get all metallicities and ages from BASTI/MILES-Models
        if self.isochrone == 'Padova00':
            # These are the available metallicities and ages for Padova00 (see http://research.iac.es/proyecto/miles/pages/ssp-models.php)
            met_models = np.array([-2.32, -1.71, -1.31, -0.71, -0.40, 0.0, 0.22])
            age_models = np.array([0.063, 0.071, 0.079, 0.089, 0.1, 0.11, 0.13,
                                   0.14, 0.16, 0.18, 0.20, 0.22, 0.25, 0.28, 0.32,
                                   0.35, 0.4, 0.45, 0.5, 0.56, 0.63, 0.71, 0.79,
                                   0.89, 1.0, 1.12, 1.26, 1.41, 1.58, 1.78, 2.0,
                                   2.24, 2.51, 2.82, 3.16, 3.55, 3.98, 4.47, 5.01,
                                   5.62, 6.31, 7.08, 7.94, 8.91, 10.0, 11.22, 12.59,
                                   14.13, 15.85, 17.78])
        else:
            # These are the available metallicities and ages for BaSTI (see http://research.iac.es/proyecto/miles/pages/ssp-models.php)
            met_models = np.array([-2.27, -1.79, -1.49, -1.26, -0.96,
                                   -0.66, -0.35, -0.25, 0.06, 0.15, 0.26, 0.40])
            age_models = np.array([0.03, 0.04, 0.05, 0.06, 0.07,
                                   0.08, 0.09, 0.1, 0.15, 0.2, 0.3, 0.35, 0.4,
                                   0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5,
                                   1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75,
                                   4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5,
                                   9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5,
                                   13.0, 13.5, 14.0])

        # make bins
        age1 = np.arange(4.25, 14.5, 0.5)
        age2 = np.arange(1.125, 4.125, 0.25)
        age3 = np.arange(0.55, 1.05, 0.1)
        age4 = np.arange(0.125, 0.525, 0.05)
        age5 = np.arange(0.035, 0.105, 0.01)
        age6 = min(self.age_tag)

        if age6 < age_models[0]:
            age_bin = np.hstack((age6, age5, age4, age3, age2, age1))
        else:
            age_bin = np.hstack((age5, age4, age3, age2, age1))

        met1 = np.arange(0.005, 0.405, 0.01)
        met2 = np.arange(-2.265, 0.005, 0.01)
        met3 = max(self.met_tag)
        met4 = min(self.met_tag)

        if (met3 > met_models[-1] and met4 < met_models[0]):
            met_bin = np.hstack((met4, met2, met1, met3))
        elif (met3 > met_models[-1] and met4 >= met_models[0]):
            met_bin = np.hstack((met2, met1, met3))
        elif (met3 <= met_models[-1] and met4 < met_models[0]):
            met_bin = np.hstack((met4, met2, met1))
        else:
            met_bin = np.hstack((met2, met1))

        # make grid
        age_grid = age_models
        met_grid = np.hstack((np.arange(-2.27, 0.0, 0.01), np.arange(0.0, 0.41, 0.01)))

        return age_models, met_models, age_grid, met_grid, age_bin, met_bin


    def flux_table(self, imf_type, imf_slope, alpha):
        '''
        this functions makes a flux table for all age/met combinations
        of the specified MILES library and stores it in a txt- and pickle-file
        '''
        print(f'Computing fluxes for IMF-type: {imf_type}, IMF-slope: {imf_slope}, alpha-abundance: {alpha}')
        #get all metallicities and ages from MILES-Models for chosen IMF
        #df = pd.read_csv('milesdata/%s_%s_baseFe_agemetML.txt'%(isochrone,imf_type),header=0, delim_whitespace=True)
        df = pd.read_csv('%sout_phot_%s_%s.txt'%(self.path_mass_to_light,imf_type,self.isochrone.upper()),header=0, delim_whitespace=True)
        df_slope = df[df['slope']==imf_slope]
        df_models = df_slope[['[M/H]','Age']].copy()
        df_models.loc[:,'[M/H]'] = np.around(df_models.loc[:,'[M/H]'], decimals=2)
        df_models = df_models.reset_index(drop=True)
        age_column = df_models['Age'].values
        met_column = df_models['[M/H]'].values

        if alpha == 'baseFe':
            temp = 'iTp0.00'
            alpha_temp = alpha

        elif alpha == 'p0.00':
            temp = 'iTp0.00'
            alpha_temp = 'Ep0.00'
        else:
            temp = 'iTp0.40'
            alpha_temp = 'Ep0.40'

        # make a flux-array in the order of going through all ages with first met, second met etc.
        column_range = np.arange(0, len(met_column), 1)
        for i in column_range:
            if met_column[i] < 0:
                met_string = 'm%4.2f' %(abs(met_column[i]))
            else:
                met_string = 'p%4.2f' %(met_column[i])
            name = '%sMILES_%s_%s_%s/M%s%4.2fZ%sT%07.4f_%s_%s.fits'%(self.path_models,
                self.isochrone, imf_type,alpha_temp, imf_type, imf_slope,
                met_string, age_column[i], temp, alpha_temp)
            hdu = fits.open(name)
            flux = hdu[0].data
            if i == 0:
                models_column = flux
            else:
                models_column = np.vstack((models_column, flux))

        # save in text-file
        np.savetxt('%s%s_%s%4.2f_%s_flux.txt'%(self.path_output,
            self.isochrone, imf_type, imf_slope, alpha), models_column)

        # save in pickle-file
        pickle.dump(models_column, open('%s%s_%s%4.2f_%s_flux.p'%(self.path_output,
            self.isochrone, imf_type, imf_slope, alpha), 'wb'))

    def finer_model_grid(self, imf_type, imf_slope, alpha, age_grid, met_grid, age_models, met_models):
        '''
        this function makes a finer model grid if the original model
        grid is too coarse for your purpose it stores all new age/met
        combinations and the corresponding flux in a txt- and pickle-file
        '''
        print(f'Refining fluxes for IMF-type: {imf_type}, IMF-slope: {imf_slope}, alpha-abundance: {alpha}')
        #dataframe for MILES models
        #all metallicity/age combinations stored in one array
        df = pd.read_csv('%sout_phot_%s_%s.txt' %(self.path_mass_to_light,imf_type, self.isochrone.upper()), header=0, delim_whitespace=True)
        df_slope = df[df['slope']==imf_slope]
        df_models = df_slope[['[M/H]','Age']].copy()
        df_models.loc[:,'[M/H]'] = np.around(df_models.loc[:,'[M/H]'], decimals=2)
        df_models = df_models.reset_index(drop=True)

        #corresponding fluxes stored in one array (made with flux_table)
        models_column = pickle.load(open('%s%s_%s%4.2f_%s_flux.p' %(self.path_output,self.isochrone , imf_type , imf_slope , alpha), 'rb'))

        #set-up defaults
        index = np.arange(0,4300,1)
        grid = np.zeros((len(met_grid),4300))

        '''
        interpolation to new metallicity grid
        attention: this sorts the fluxes like this: first all 
        fluxes for one age & all metallicites, 
        then all fluxes for second age & all metallicities and so on
        '''
        for i, age in enumerate(age_grid):  # extract all models with same age
            x = df_models[df_models['Age'] == age]  # extract corresponding fluxes
            temp = models_column[x.index.values]  # interpolate for every j-th wavelength
            for j in index:
                y = temp[:,j]  # extract flux at j-th wavelength
                f = interpolate.UnivariateSpline(met_models,y)
                # interpolate metallicity with full model grid
                grid[:,j] = f(met_grid) # extract interpolated flux at j-th wavelength for every metallicity in new grid
            # stack fluxes
            if i == 0:
                full_grid = grid
            else:
                full_grid = np.vstack((full_grid, grid))

        #make new dataframe for finer grid and fill with new values
        df_finer = pd.DataFrame(np.nan, index=range(len(met_grid)*len(age_grid)), columns=df_models.columns.values.tolist())
        met_column = met_grid
        age_column = np.zeros(len(met_grid))+age_grid[0]
        for a in age_grid[1:]:
            temp = np.zeros(len(met_grid))+a
            age_column = np.hstack((age_column,temp))
            met_column = np.hstack((met_column,met_grid))
        df_finer.loc[:,'[M/H]'] = met_column
        df_finer.loc[:,'Age'] = age_column

        # save in text-file
        df_finer.to_csv('%s%s_%s%4.2f_%s_agemet_finer.txt'%(self.path_output,
            self.isochrone, imf_type, imf_slope, alpha), index=False, sep=' ')
        np.savetxt('%s%s_%s%4.2f_%s_flux_finer.txt' %(self.path_output,
            self.isochrone, imf_type, imf_slope, alpha), full_grid)

        # save in pickle-file
        df_finer.to_pickle('%s%s_%s%4.2f_%s_agemet_finer.pkl' %(self.path_output,
            self.isochrone, imf_type, imf_slope, alpha))
        pickle.dump(full_grid, open('%s%s_%s%4.2f_%s_flux_finer.p' %(self.path_output,
            self.isochrone, imf_type, imf_slope, alpha), 'wb'))

    def make_spectrum(self, name_particles, imf_type, imf_slope, alpha, age_grid, met_grid, age_bin, met_bin):
        '''
        this function composes the actual eagle spectrum and stores it in a fits_file
        '''
        particles_range = np.arange(0, len(self.age_tag), 1)
        print(f'Computing spectra for {len(particles_range)} particles')

        # compute mass weighted meand and standard deviation of age & met for later
        average_age = np.sum(self.mass_tag * self.age_tag) / np.sum(self.mass_tag)
        average_met = np.sum(self.mass_tag * self.met_tag) / np.sum(self.mass_tag)

        spread_age = np.sqrt((np.sum(self.mass_tag)*np.sum(self.mass_tag*(self.age_tag-average_age)**2))/
                             ((np.sum(self.mass_tag))**2-np.sum(self.mass_tag**2)))
        spread_met = np.sqrt((np.sum(self.mass_tag)*np.sum(self.mass_tag*(self.met_tag-average_met)**2))/
                             ((np.sum(self.mass_tag))**2-np.sum(self.mass_tag**2)))

        '''
        dataframe for model grid note: data in df_grid corresponds to the actual grid
        points, bins correspond to boundaries(they lie exactly between the actual grid points);
        grid points and bins are specified by the mode you’ve chosen
        '''

        df_grid = pd.read_pickle('%s%s_%s%4.2f_%s_agemet_finer.pkl'%(self.path_output,
            self.isochrone, imf_type, imf_slope, alpha))
        flux_column = pickle.load(open('%s%s_%s%4.2f_%s_flux_finer.p'%(self.path_output,
            self.isochrone, imf_type, imf_slope, alpha),'rb'))

        '''    
        binning: eagle data that falls into a bin will be assigned the corresponding
        grid point binnumber gives you the number of the bin in which the value
        falls (has same length as input data and values from 1 to length of bin)
        '''

        age_statistics, age_bin_edges, age_binnumber = stats.binned_statistic(
            self.age_tag, None, 'count', bins=age_bin)

        met_statistics, met_bin_edges, met_binnumber = stats.binned_statistic(
            self.met_tag, None, 'count', bins=met_bin)

        # defaults for eagle spectrum
        spectrum_mass = np.zeros(4300)

        # empty array for masses
        masses = np.zeros([len(age_grid),len(met_grid)])

        # create spectrum and save masses
        for i in tqdm(particles_range):
            # assign the corresponding grid value depending on in which bin the actual eagle data landed
            age = age_grid[age_binnumber[i] - 1]
            met = met_grid[met_binnumber[i] - 1]
            # search for right age/met combination to extract flux
            x = df_grid[(df_grid['Age'] == age) & (df_grid['[M/H]'] == met)]
            flux = flux_column[x.index[0]]
            spectrum_mass = spectrum_mass + self.mass_tag[i] * flux
            # search for right age/met combination to extract mass
            x_mass = np.where(age_grid == age)
            y_mass = np.where(met_grid == met)
            masses[x_mass, y_mass] = self.mass_tag[i]

        # write spectrum to fits-file
        hdu = fits.PrimaryHDU(spectrum_mass)
        hdu.header.append(('CRVAL1', 3540.5, 'starting wavelength'), end=True)
        hdu.header.append(('CDELT1', 0.9, 'angstrom/pixel scale'), end=True)
        hdu.header.append(('NAXIS1', 4300, 'number of pixels'), end=True)
        hdu.header.append(('age', average_age, 'mass-weighted average age'), end=True)
        hdu.header.append(('sig_age', spread_age, 'mass-weighted SD of age'), end=True)
        hdu.header.append(('met', average_met, 'mass-weighted average metallicity'), end=True)
        hdu.header.append(('sig_met', spread_met, 'mass-weighted SD of metallicity'), end=True)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto('result_%s_%s_%s%4.2f_%s.fits' % (
            name_particles, self.isochrone, imf_type, imf_slope, alpha))

        # write masses to csv file
        massdata = pd.DataFrame(masses.T, columns=age_grid.round(3), index=met_grid.round(3))
        massdata.to_csv('masses_agemet_%s_%s_%s%4.2f_%s.csv' % (
            name_particles, self.isochrone, imf_type, imf_slope, alpha))

    '''
    if you want to add noise to any spectrum, 
    this function will do it for you it creates a new fits-file with the noisy flux
    '''

    def add_noise(self, name_particles, imf_type, imf_slope, alpha, SNR):
        hdu = fits.open('result_%s_%s_%s%4.2f_%s.fits' %(
            name_particles , self.isochrone , imf_type , imf_slope , alpha))
        flux = hdu[0].data
        header = hdu[0].header
        for i in SNR:
            noise = np.random.standard_normal(len(flux))*flux*(1./i)
            noiseyflux = flux + noise
            hdu0 = fits.PrimaryHDU(noiseyflux)
            hdu0.header.append(('CRVAL1', 3540.5, 'starting wavelength'), end=True)
            hdu0.header.append(('CDELT1', 0.9, 'angstrom/pixel scale'), end=True)
            hdu0.header.append(('NAXIS1', 4300, 'number of pixels'), end=True)
            hdu0.header.append(('age', header['age'], 'mass-weighted average age'), end=True)
            hdu0.header.append(('sig_age', header['sig_age'], 'mass-weighted SD of age'), end=True)
            hdu0.header.append(('met', header['met'], 'mass-weighted average metallicity'), end=True)
            hdu0.header.append(('sig_met', header['sig_met'], 'mass-weighted SD of metallicity'), end=True)
            hdulist = fits.HDUList([hdu0])
            hdulist.writeto('result_%s_%s_%s%4.2f_%s_SNR%2d.fits' %(
                name_particles, self.isochrone, imf_type, imf_slope, alpha, i))
            hdu.close()

if __name__ == '__main__':
    MakeMock('default_config.yml')