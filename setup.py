#!/usr/bin/env python
# encoding: utf-8
from setuptools import setup, find_packages
import glob

setup(
    name="synthsb",
    version="0.1",
    packages=find_packages(),
    scripts=glob.glob('scripts/*.py'),

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    # install_requires = ['docutils>=0.3'],

    package_data={
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
        # # And include any *.msg files found in the 'hello' package, too:
        # 'hello': ['*.msg'],
    },

    # metadata for upload to PyPI
    author="Jonathan Sick",
    author_email="jonathansick@mac.com",
    description="Estimate surface brightness from resolved star "
                "catalogs in M31",
    license="BSD",
    keywords="astronomy",
    url="http://jonathansick.ca",   # project home page, if any

    # could also include long_description, download_url, classifiers, etc.
)
