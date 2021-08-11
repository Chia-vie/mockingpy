import setuptools
# read in required packages from requirements.txt
with open("requirements.txt", "r") as fh:
    requirements = fh.read().splitlines()

setuptools.setup(
    name="mockingpy",
    version="0.0.1",
    author='Christine Ackerl, based on work by Alina Boecker',
    author_email='christine.ackerl@univie.ac.at',
    description='Create galaxy mock spectra (and more) using Python.',
    long_description='''This package is currently under development. Stay tuned for updates.''',
    url="https://github.com/Chia-vie/mockingpy",
    license='BSD',
    keywords='astronomy, astrophysics, galaxies, spectra,...',
    packages=setuptools.find_packages(),
    package_data={'mockingpy': ['default_config.yml']},
    classifiers=[
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: BSD License',],
    install_requires=requirements,
    include_package_data=True
)
