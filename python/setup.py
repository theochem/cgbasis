#!/usr/bin/env python

import os

from setuptools import setup, Extension

from Cython.Build import cythonize
import numpy as np


def get_version():
    """Get the version string set by Travis, else default to version 0.0.0."""
    return os.environ.get("PROJECT_VERSION", "0.0.0")


def get_cxxflags():
    """If the CXXFLAGS variable is defined (clang/osx) then get it."""
    return os.environ.get("CXXFLAGS", "-std=c++11 -Wall").split()


def get_include_path():
    """If the PREFIX variable is defined (conda) then get the conda include prefix."""
    prefix = os.environ.get("PREFIX", "")
    if prefix:
        return [os.path.join(prefix, "include")]
    else:
        return []


ext_modules = [Extension("gbasis.cext_common", ["gbasis/cext_common.pyx"],
                         include_dirs=[np.get_include()] + get_include_path(),
                         libraries=["gbasis"],
                         extra_compile_args=get_cxxflags(),
                         language="c++",
                         ),
               Extension("gbasis.cext1", ["gbasis/cext1.pyx"],
                         include_dirs=[np.get_include()] + get_include_path(),
                         libraries=["gbasis"],
                         extra_compile_args=get_cxxflags(),
                         language="c++",
                         ),
               Extension("gbasis.cext2", ["gbasis/cext2.pyx"],
                         include_dirs=[np.get_include()] + get_include_path(),
                         libraries=["gbasis"],
                         extra_compile_args=get_cxxflags(),
                         language="c++",
                         ),
               Extension("gbasis.cext_grids", ["gbasis/cext_grids.pyx"],
                         include_dirs=[np.get_include()] + get_include_path(),
                         libraries=["gbasis"],
                         extra_compile_args=get_cxxflags(),
                         language="c++",
                         ),
               Extension("gbasis.cext_sparse", ["gbasis/cext_sparse.pyx"],
                         include_dirs=[np.get_include()] + get_include_path(),
                         libraries=["gbasis"],
                         extra_compile_args=get_cxxflags(),
                         language="c++",
                         ),
               Extension("gbasis.test.cext", ["gbasis/test/cext.pyx"],
                         include_dirs=[np.get_include()] + get_include_path(),
                         libraries=["gbasis"],
                         extra_compile_args=get_cxxflags(),
                         language="c++",
                         ),
               ]

setup(
    name="gbasis",
    version=get_version(),
    description="",
    long_description="A Python interface for GBasis(c++).",
    author="Toon Verstraelen",
    author_email="Toon.Verstraelen@UGent.be",
    url="https://github.com/theochem/gbasis",
    package_dir={"gbasis": "gbasis"},
    packages=["gbasis", "gbasis.test", "gbasis.bsets", "gbasis.test.cached"],
    package_data={"gbasis.bsets": ["*"],
                  "gbasis.test.cached": ["*/*"]},
    ext_modules=cythonize(ext_modules, include_path=["gbasis/pxds", "gbasis", "gbasis/test/pxds"]),
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
