# setup.py: The setup file
from setuptools import setup, Extension


setup(name='spkmeans',
      version='1.0',
      author="Firas, Kareen",
      description='spk algorethim',
      ext_modules=[Extension('spkmeans',['spkmeansmodule.c', 'spkmeans.c'],),])
