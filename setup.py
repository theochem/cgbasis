#!/usr/bin/env python

import os
import glob
import Cython.Build
import numpy as np
from setuptools import setup, Extension


def get_version():
    """Load the version from version.py, without importing it.
    This function assumes that the last line in the file contains a variable defining the
    version string with single quotes.
    """
    with open('gbasis/version.py', 'r') as f:
        return f.read().split('=')[-1].replace('\'', '').strip()


def get_cxxflags():
    """If the CXXFLAGS variable is defined (clang/osx) then get it."""
    return os.environ.get("CXXFLAGS", "").split()


def get_readme():
    with open('README.rst') as f:
        return f.read()


setup(
    name="gbasis",
    version=get_version(),
    description="",
    long_description=get_readme(),
    author="Toon Verstraelen",
    author_email="Toon.Verstraelen@UGent.be",
    url="https://github.com/theochem/gbasis",
    package_dir={"gbasis": "gbasis"},
    packages=["gbasis", "gbasis.test"],
    cmdclass={"build_ext": Cython.Build.build_ext},
    ext_modules=[Extension(
        "gbasis.cext",
        sources=["gbasis/cext.pyx"] + glob.glob("gbasis/*.cpp"),
        depends=glob.glob("gbasis/*.hpp"),
        include_dirs=[np.get_include()],
        libraries=["int2"],
        extra_compile_args=get_cxxflags(),
        language="c++",
        ),
    ],
    include_package_data=True,
    classifiers=[
        "Environment :: Console",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Intended Audience :: Science/Research",
    ],
    zip_safe=False,
    setup_requires=['numpy>=1.0', 'cython>=0.24.1'],
    install_requires=['numpy>=1.0', 'nose>=0.11', 'cython>=0.24.1']
)
