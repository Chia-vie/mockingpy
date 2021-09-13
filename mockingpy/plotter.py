import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import glob
import os
import pandas as pd
from matplotlib.colors import LogNorm
import vaex

class PlotSpec():
    def __init__(self,files):
        self.path = os.getcwd()
        self.filenames = [filename for filename in glob.glob(self.path+'/'+files)]
        self.colors()
        self.plot_spec()

    def colors(self):
        #purpleblue = ['#642e7c','#7251b7','#8984d6','#93bae1','#8ce2ee']*10
        #icey = ['#003279','#4386bb','#8bbedc','#cccce4','#fbb9cd','#896eb6','#6a1364','#49175f','#020316']*10
        icey2 = ['#0f142d','#112a60','#2a2e62','#205094','#50356f','#3a79b9','#8c7fb2','#84acd5','#babed9','#f0c6d6','#ead8e1']*10
        for color in icey2:
            yield color

    def labels(self):
        labels = [filename.replace(self.path+'/','') for filename in self.filenames]
        for label in labels:
            yield label

    def plot_spec(self):
        x = np.arange(3540.5,3540.5+(4300*0.9),0.9)
        fig, ax = plt.subplots(1,1,figsize=(20,5))
        colors = self.colors()
        labels = self.labels()
        for filename in self.filenames:
            self.data = fits.open(filename)[0].data
            ax.plot(x, self.data, color=next(colors), label=next(labels))
        ax.set_xlabel('wavelength [nm]')
        ax.set_ylabel('flux')
        ax.set_title('Mock spectra')
        leg=ax.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()
        for line in leg.get_lines():
            line.set_linewidth(10)
        plt.show()

class AgeMet():
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

class Maps():
    def __init__(self, file):
        self.file = file
        self.path = os.getcwd()
        self.data = vaex.from_ascii(self.path + '/' + self.file)

    def make(self, valcol, selcol=None, seltype=None, limit=None, exclude=None, limits=[[-200, 200], [-200, 200]], shape=(50, 50)):
        '''
        :param valcol: The column name of the values to map
        :return: None, plots heatmap
        '''
        if not seltype == None:
            if seltype.lower() == 'l':
                selection = self.data[selcol] > limit
                title = f'{selcol}>{limit}'
            elif seltype.lower() == 'u':
                selection = self.data[selcol] < limit
                title = f'{selcol}<{limit}'
            elif seltype.lower() == 'e':
                selection = self.data[selcol] != exclude
                title = f'Excluded: {selcol} = {exclude}'
        else:
            selection = None
            title = 'All'
        self.valcol = valcol
        self.values= self.data.mean(self.data[self.valcol], binby=[self.data.x, self.data.y], limits=limits, shape=shape, selection=selection)
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        map =  plt.imshow(self.values.T, origin='lower', extent=[-200,200,-200,200])
        ax.set_title(title)
        ax.set_xlabel('x [kpc]')
        ax.set_ylabel('y [kpc]')
        plt.colorbar(map, label=self.valcol)
        plt.show()

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








