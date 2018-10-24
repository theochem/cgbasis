#!/usr/bin/env python

import os
import glob
import Cython.Build
import numpy as np
from setuptools import setup, Extension


def get_version():
    """Get the version string set by Travis, else default to version 0.0.0"""
    return os.environ.get("PROJECT_VERSION", "0.0.0")


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
    packages=["gbasis", "gbasis.test", "gbasis.bsets", "gbasis.test.cached"],
    package_data={"gbasis.bsets": ["*"],
                  "gbasis.test.cached": ["*/*"]},
    cmdclass={"build_ext": Cython.Build.build_ext},
    ext_modules=[Extension(
            "gbasis.cext",
            sources=["gbasis/cext.pyx"] + glob.glob("gbasis/*.cpp"),
            depends=glob.glob("gbasis/*.h") + glob.glob("gbasis/*.pxd"),
            include_dirs=[np.get_include()],
            libraries=["int2"],
            extra_compile_args=get_cxxflags(),
            language="c++",
        ),
        Extension(
            "gbasis.test.cext",
            sources=["gbasis/test/cext.pyx"] + glob.glob("gbasis/*.cpp"),
            depends=glob.glob("gbasis/*.h") + glob.glob("gbasis/*.pxd"),
            include_dirs=[np.get_include(), "gbasis/"],
            libraries=["int2"],
            extra_compile_args=get_cxxflags(),
            language="c++",
        ),
        Extension(
            "gbasis.cext_fns",
            sources=["gbasis/cext_fns.pyx"] + ["gbasis/fns.cpp", "gbasis/calc.cpp",
                                               "gbasis/common.cpp", "gbasis/iter_pow.cpp"],
            depends=["gbasis/fns.h", "gbasis/calc.h", "gbasis/common.h", "gbasis/iter_pow.h"]
                    + ["gbasis/fns.pxd, gbasis/cext.pxd"],
            include_dirs=[np.get_include()],
            libraries=["int2"],
            extra_compile_args=get_cxxflags(),
            language="c++",
        ),
        Extension(
            "gbasis.cext_ints",
            sources=["gbasis/cext_ints.pyx"] + ["gbasis/ints.cpp", "gbasis/calc.cpp"],
            depends=["gbasis/ints.h", "gbasis/calc.h"] + ["gbasis/ints.pxd, gbasis/cext.pxd"],
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
