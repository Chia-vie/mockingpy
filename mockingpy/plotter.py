import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import glob
import os
import pandas as pd
from matplotlib.colors import LogNorm
import vaex

class Spectra():
    def __init__(self,files):
        self.path = os.getcwd()
        self.filenames = [filename for filename in glob.glob(self.path+'/'+files)]
        #for filename in self.filenames:
        self.spectra = [fits.open(filename)[0].data for filename in self.filenames]
        self.names = [filename.replace(self.path+'/','') for filename in self.filenames]
        self.plot_spec()

    def colors(self, palette):
        #purpleblue = ['#642e7c','#7251b7','#8984d6','#93bae1','#8ce2ee']*10
        #icey = ['#003279','#4386bb','#8bbedc','#cccce4','#fbb9cd','#896eb6','#6a1364','#49175f','#020316']*10
        icey2 = ['#0f142d','#112a60','#2a2e62','#205094','#50356f','#3a79b9','#8c7fb2','#84acd5','#babed9','#f0c6d6','#ead8e1']*10
        simple = ['purple','royalblue','black']*10
        if palette == 'icey2':
            for color in icey2:
                yield color
        else:
            for color in simple:
                yield color

    def labels(self):
        for name in self.names:
            yield name

    def plot_spec(self, color='icey2'):
        x = np.arange(3540.5,3540.5+(4300*0.9),0.9)
        fig, ax = plt.subplots(1,1,figsize=(20,5))
        colors = self.colors(color)
        labels = self.labels()
        for spectrum in self.spectra:
            #self.data = fits.open(filename)[0].data
            ax.plot(x, spectrum, color=next(colors), label=next(labels))
        ax.set_xlabel(r'wavelength [$\AA$] ')
        ax.set_ylabel('flux')
        ax.set_title('Mock spectra')
        leg=ax.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()
        for line in leg.get_lines():
            line.set_linewidth(10)
        plt.show()

    def add_subset(self, filenum1, filenum2, operation, name):
        if operation.lower() == 'add':
            self.spectra.append(self.spectra[filenum1] + self.spectra[filenum2])
        elif operation.lower() in ['sub','subtract']:
            self.spectra.append(self.spectra[filenum1] - self.spectra[filenum2])
        self.names.append(name)
        print(f'New list of subsets: {self.names}')

    def compare(self, labels, filenums=None, total=0, palette='simple'):
        '''
        Args:
            labels:
            total: the list index in filenames of the file containing the total flux. Defaults to 0.

        Returns:
        '''
        if filenums == None:
            filenums = [i for i,file in enumerate(self.spectra)]
        x = np.arange(3540.5, 3540.5 + (4300 * 0.9), 0.9)
        fig = plt.figure(figsize=[20, 15])
        fig.suptitle('Mock spectra', size=20)
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        #fig, ax = plt.subplots(1, 1, figsize=(20, 5))
        colors = self.colors(palette)
        all = self.spectra[total]
        for i in filenums:
            sub = self.spectra[i]
            label = labels[i]
            #residual = all - sub
            relative = sub/all
            median = np.median(relative)
            comparison = relative/median
            currentcolor = next(colors)
            ax1.plot(x, comparison, color=currentcolor, label=label)
            ax2.plot(x, sub, color=currentcolor, label=label)
            ax3.plot(x, relative, color=currentcolor, label=label)
        ax3.set_xlabel(r'wavelength [$\AA$] ')
        ax1.set_ylabel('Fraction of total flux devided by median')
        ax2.set_ylabel('Flux')
        ax3.set_ylabel('Fraction of total flux')
        leg = ax1.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()
        for line in leg.get_lines():
            line.set_linewidth(10)
        plt.show()

    def residuals_from_files(self, labels, total=0):
        '''
        Args:
            labels:
            total: the list index in filenames of the file containing the total flux. Defaults to 0.

        Returns:
        '''
        x = np.arange(3540.5, 3540.5 + (4300 * 0.9), 0.9)
        fig, ax = plt.subplots(1, 1, figsize=(20, 5))
        colors = self.colors()
        all = fits.open(self.filenames[total])[0].data
        for i,filename in enumerate(self.filenames[1:]):
            sub = fits.open(filename)[0].data
            label = labels[i]
            res = (all-sub)/all
            ax.plot(x, res, color=next(colors), label=label)
        ax.set_xlabel(r'wavelength [$\AA$] ')
        ax.set_ylabel('Fraction of total flux')
        ax.set_title('Mock spectra')
        leg = ax.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()
        for line in leg.get_lines():
            line.set_linewidth(10)
        plt.show()

    def residuals_from_data(self,datasets, labels, colors='simple'):
        x = np.arange(3540.5, 3540.5 + (4300 * 0.9), 0.9)
        fig, ax = plt.subplots(1, 1, figsize=(20, 5))
        colors = self.colors(colors)
        all = datasets[0]
        for i,dataset in enumerate(datasets[1:]):
            label = labels[i]
            res = (all-dataset)/all
            ax.plot(x, res, color=next(colors), label=label)
        ax.set_xlabel(r'wavelength [$\AA$] ')
        ax.set_ylabel('Fraction of total flux')
        ax.set_title('Mock spectra')
        leg = ax.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()
        for line in leg.get_lines():
            line.set_linewidth(10)
        plt.show()

class ViewAgeMet():
    def __init__(self, file):
        self.file = file
        self.path = os.getcwd()
        self.data = pd.read_csv(self.path + '/' + self.file, header=0, index_col=0)
        self.plot_age_met()

    def plot_age_met(self):
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        map = ax.pcolor(self.data.columns.astype(float), self.data.index, self.data, cmap='magma', norm=LogNorm())
        ax.set_title('Age-Metallicity Distribution')
        ax.set_xlabel('Age')
        ax.set_ylabel('Metallicity')
        plt.colorbar(map, label='Mass')
        plt.show()

class ViewParticles():
    def __init__(self, file):
        self.file = file
        self.path = os.getcwd()
        self.data = vaex.from_ascii(self.path + '/' + self.file)
        self.m_sum = self.data.sum(self.data.mass)

    def map(self, valcol, selection=False, limits=[[-200, 200], [-200, 200]], shape=(50, 50)):
        '''
        :param valcol: The column name of the values to map
        :param selection: The expression by which to select, e.g.: 'mass < 10**9'
        :return: None, plots heatmap
        '''
        self.valcol = valcol
        self.values = self.data.mean(self.data[self.valcol], binby=[self.data.x, self.data.y], limits=limits, shape=shape, selection=selection)
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        map = plt.imshow(self.values.T, origin='lower', extent=[-200,200,-200,200])

        if selection:
            title = selection
        else:
            title = 'All'

        ax.set_title(title)
        ax.set_xlabel('x [kpc]')
        ax.set_ylabel('y [kpc]')
        plt.colorbar(map, label=self.valcol)
        plt.show()

    def mass_fractions(self, groupeddata, binnum=20, log_spacing=True):
        '''
        :param groupeddata: tuple of vaex dataframes grouped by individual satellites
        :param binnum: number of bins
        :return:
        '''
        if log_spacing == True:
            colors=['purple','royalblue','black']
            for i, data in enumerate(groupeddata):
                if i == 0:
                    counts, logbins = np.histogram(np.log10(data.mass.values), bins=binnum)
                    label = 'all'
                else:
                    label = f'subset {i}'
                logmassbins = []
                for edge in logbins[1:]:
                    binmembers = data[np.log10(data.mass) <= edge]
                    binmass = binmembers.sum(binmembers.mass)
                    logmassbins.append(binmass)
                norm_logmassbins = logmassbins/self.m_sum
                x = logbins[1:]
                y = norm_logmassbins
                plt.bar(x, y, align='edge', color=colors[i], label=label, width=-.5)
            plt.xlim(logbins[4], logbins[-1])
            plt.legend()
            plt.xlabel('Satellite log(m_stellar)')
            plt.ylabel('Fraction of halo m_stellar')
            plt.show()
        else:
            pass

    def subset(self, selcol, expression, limit):
        if expression == '<':
            sel = self.data[self.data[selcol] < limit]
        elif expression == '>':
            sel = self.data[self.data[selcol] > limit]
        elif expression == '==':
            sel = self.data[self.data[selcol] == limit]
        else:
            sel = self.data[self.data[selcol] != limit]
        out = sel.groupby(by='nsat').agg({'mass': 'sum',
                                  'feh': ['mean', 'std'], 'age': ['mean','std']})
        out['mass/halo_mass'] = out.mass / self.m_sum
        return out

    def satellites(self, name):
        self.name = name
        # Group by nsat, raname nsat column, join grouped dataset and old one
        groupeddata = self.data.groupby(by='nsat',
                        agg={'lboltot': vaex.agg.sum('lbol')})
        groupeddata.rename('nsat','groupednsat')
        joined=self.data.join(groupeddata, left_on='nsat', right_on='groupednsat')

        # remove renamed nsat column, is not needed anymore.
        joined.drop('groupednsat')

        # calculate the weights
        joined['lumweight'] = joined['lbol']/joined['lboltot']

        # calculate the radii and weighted radii
        joined['radius'] = np.sqrt(joined['x']**2+joined['y']**2+joined['z']**2)
        joined['weightrad'] = joined['lumweight']*joined['radius']

        #extract the data we want
        result=joined.groupby(by='nsat', agg={'mean_weightrad': vaex.agg.sum('weightrad'),
                             'mean_rad': vaex.agg.mean('radius'),
                             'tac': vaex.agg.mean('tac'),
                             'mtot': vaex.agg.mean('mtot')})
        tac = result.evaluate("tac")
        meanweightrad = result.evaluate("mean_weightrad")
        #meanrad = result.evaluate('mean_rad')
        mtot = result.evaluate('mtot')

        #plot weighted radii against tac and mtot
        fig = plt.figure(figsize = [18,7])
        fig.suptitle(name, size = 20)

        ax1 = fig.add_subplot(121)
        ax1.set_xlabel('$r_{mean}$ [kpc]', size = 20)
        ax1.set_ylabel('$t_{acc}$ [Gyr ago]', size = 20)
        ax1.set_ylim(0,12.8)
        ax1.set_xlim(0,600)

        ax2 = fig.add_subplot(122)
        ax2.set_xlabel('$r_{mean}$ [kpc]', size = 20)
        ax2.set_ylabel('$log(m_{tot})$ [$M_{sun}$]', size = 20)
        ax2.set_yscale('log')
        ax2.set_xlim(0,600)
        ax2.set_ylim(10**8,10**11)

        a = ax1.scatter(meanweightrad, tac, c=np.log10(mtot), s=30, vmin=8, vmax=11, cmap='inferno')
        b = ax2.scatter(meanweightrad, mtot, c=tac, s=30, vmin=0, vmax=12.8,cmap='inferno')
        fig.colorbar(a, ax=ax1)
        fig.colorbar(b, ax=ax2).set_label(label='$t_{acc}$ [Gyr ago]',size=20)
        plt.subplots_adjust(wspace=0.1)
        plt.show()