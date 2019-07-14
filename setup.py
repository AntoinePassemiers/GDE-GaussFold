# -*- coding: utf-8 -*-
# setup.py
# author : Antoine Passemiers

from setuptools import setup


setup(
    name='gaussfold',
    version='1.0.0',
    description='Fast Template-free Protein Modelling',
    url='https://github.com/AntoinePassemiers/GDE-GaussFold',
    author='Antoine Passemiers',
    author_email='apassemi@ulb.ac.be',
    packages=['gaussfold', 'gaussfold.model', 'gaussfold.aa', 'gaussfold.atom', 'gaussfold.chain',
              'gaussfold.constraints'])
