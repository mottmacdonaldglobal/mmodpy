# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 15:26:36 2019

@author: kevin.stanton
"""

from setuptools import setup, find_packages

setup(name='mmodpy',
      version='1.0.1',
      description='A collection of automation scripts to expedite LS-DYNA analytics.',
      url='https://github.com/mottmacdonaldglobal/mmodpy.git',
      author='Kevin Stanton',
      author_email='kevin.stanton@mottmac.com',
      install_requires=['subprocess'],
      packages=find_packages(),
      etc_files=[('mmodpy/etc',
                   ['mmodpy/etc/Primer-Automesh_Coarse_15b.js',
                    'mmodpy/etc/PR_P500_Foundation_Setup_005.js']),
                  ('mmodpy/sample-inputs',
                   ['mmodpy/sample-inputs/foundations.key',
                    'mmodpy/sample-inputs/Input.xlsx']),
                  ('mmodpy/sample-inputs/th_msec2',
                   ['mmodpy/sample-inputs/th_msec2/TH1.csv',
                    'mmodpy/sample-inputs/th_msec2/TH2.csv',
                    'mmodpy/sample-inputs/th_msec2/TH3.csv',
                    'mmodpy/sample-inputs/th_msec2/TH4.csv'])],
      include_package_etc=True)