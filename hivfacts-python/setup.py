#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from glob import glob
from shutil import copy2

import setuptools
from setuptools.command.build_py import build_py

version = '2020.4'


class CopyHIVFactsDataCmd(build_py):

    def run(self):
        for fname in glob('../data/aapcnt/*.json'):
            copy2(fname, 'hivfacts/data/aapcnt/')
        for fname in glob('../data/codonpcnt/*.json'):
            copy2(fname, 'hivfacts/data/codonpcnt/')
        for fname in glob('../data/apobecs/*.json'):
            copy2(fname, 'hivfacts/data/apobecs/')
        for fname in glob('../data/*.json'):
            copy2(fname, 'hivfacts/data/')
        super(CopyHIVFactsDataCmd, self).run()


setup_params = dict(
    name='hivfacts',
    version=version,
    url='https://github.com/hivdb/hivfacts/tree/master/hivfacts-python',
    author='Philip Tzou',
    author_email='philiptz@stanford.edu',
    description='Amino acid prevalence data of HIV-1 pol',
    packages=['hivfacts', 'hivfacts/data'],
    include_package_data=True,
    cmdclass={
        'copy_data': CopyHIVFactsDataCmd
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: Public Domain',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    zip_safe=True)

if __name__ == '__main__':
    setuptools.setup(**setup_params)
