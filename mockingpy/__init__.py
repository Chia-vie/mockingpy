import os
from shutil import copyfile
from pkg_resources import resource_stream
from .mocking import MakeMock
from .plotter import PlotSpec, AgeMet

default_config = resource_stream('mockingpy', 'default_config.yml').name
path = os.getcwd()
newfile = path+'/config.yml'
if not os.path.exists(newfile):
      copyfile(default_config, newfile)
      print(f'Copied the default config file: config.yml into your current working directory: \n'
            f'{path}. \n'
            f'Please edit the file according to your needs.')