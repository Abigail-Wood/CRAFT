#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst', 'r') as fh:
    long_description = fh.read()

requirements = [
    'pandas>=0.2',
    'numpy>=1.0, <1.15',
    'scipy>=1.1',
    'PyVcf>=0.1'
]

setup(
    name='bio-craft',
    version='0.0.1',
    author='Abigail Wood',
    author_email='acjwd3@gmail.com',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    description='Credible Refinement and Annotation of Functional Targets',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    include_package_data=True,
    install_requires=requirements,
    license='LICENSE.txt',
    keywords='bio-craft',
    packages=find_packages(),
    url='https://github.com/Abigail-Wood/CRAFT',
)
