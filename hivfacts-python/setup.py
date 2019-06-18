#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from glob import glob
from shutil import copy2

import setuptools
from setuptools.command.build_py import build_py

version = '2019.6a1'


class CopyAAPcntDataCmd(build_py):

    def run(self):
        for fname in glob('../data/*.json'):
            copy2(fname, 'hivaapcnt/data/')
        super(CopyAAPcntDataCmd, self).run()


setup_params = dict(
    name='hiv-aapcnt',
    version=version,
    url='https://github.com/hivdb/hiv-aapcnt/tree/master/hiv-aapcnt-python',
    author='Philip Tzou',
    author_email='philiptz@stanford.edu',
    description='Amino acid prevalence data of HIV-1 pol',
    packages=['hivaapcnt', 'hivaapcnt/data'],
    include_package_data=True,
    cmdclass={
        'copy_aapcnt': CopyAAPcntDataCmd
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
