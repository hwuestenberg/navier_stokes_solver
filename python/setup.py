from setuptools import setup

# reading long description from file
with open('DESCRIPTION.txt') as file:
    long_description = file.read()

# specify requirements of your package here
REQUIREMENTS = [
    'numpy',
    'matplotlib',
    'pandas',
    'scipy',
]

# some more details
CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Topic :: Internet',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.8',
]

# calling the setup function 
setup(name='navier_stokes_solver',
      version='0.1',
      description='Solver for 2D Navier-Stokes equations',
      long_description=long_description,
      url='https://github.com/hwuestenberg/navier_stokes_solver',
      author='Henrik Wuestenberg',
      author_email='henrik.wuestenberg@hotmail.de',
      license='MIT',
      packages=['nssolver'],
      classifiers=CLASSIFIERS,
      install_requires=REQUIREMENTS,
      keywords='cfd fluids simulation'
      )