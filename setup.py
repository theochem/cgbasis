#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import glob
import json
import subprocess
import Cython.Build
import numpy as np
from setuptools import setup
from distutils.extension import Extension


def get_version():
    """Load the version from version.py, without importing it.
    This function assumes that the last line in the file contains a variable defining the
    version string with single quotes.
    """
    with open('meanfield/version.py', 'r') as f:
        return f.read().split('=')[-1].replace('\'', '').strip()


def get_readme():
    with open('README.rst') as f:
        return f.read()


# Library configuration functions
# -------------------------------
lib_config_keys = ["include_dirs", "library_dirs", "libraries", "extra_objects",
                   "extra_compile_args", "extra_link_args"]


def print_lib_config(heading, lib_config):
    """Print (partial) lib_config"""
    print("   %s" % heading)
    if len(lib_config) == 0:
        print("      -")
    else:
        for key, value in sorted(lib_config.items()):
            if len(value) > 0:
                print("      %s: %s" % (key, value))


def get_lib_config_setup(prefix, fn_setup_cfg):
    """Get library configuration from a setup.cfg"""
    lib_config = {}
    if os.path.isfile(fn_setup_cfg):
        config = ConfigParser.ConfigParser()
        config.read(fn_setup_cfg)
        if config.has_section(prefix):
            for key in lib_config_keys:
                if config.has_option(prefix, key):
                    value = config.get(prefix, key).strip()
                    if value is not None and len(value) > 0:
                        lib_config[key] = value.split(":")
        print_lib_config("From %s" % fn_setup_cfg, lib_config)
    else:
        print("   File %s not found. Skipping." % fn_setup_cfg)
    return lib_config


def get_lib_config_env(prefix):
    """Read library config from the environment variables"""
    lib_config = {}
    for key in lib_config_keys:
        varname = ("%s_%s" % (prefix, key)).upper()
        value = os.getenv(varname)
        if value is not None:
            lib_config[key] = value.split(":")
    print_lib_config("From environment variables", lib_config)
    return lib_config


class PkgConfigError(Exception):
    pass


def run_pkg_config(libname, option):
    """Safely try to call pkg-config"""
    try:
        return subprocess.check_output(["pkg-config", libname, "--" + option],
                                       stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        raise PkgConfigError("pkg-config did not exit properly")
    except OSError:
        raise PkgConfigError("pkg-config not installed")


def get_lib_config_pkg(libname):
    """Get library config from the pkg-config program"""
    lib_config = {
        "include_dirs": [word[2:] for word in run_pkg_config(libname, "cflags-only-I").split()],
        "library_dirs": [word[2:] for word in run_pkg_config(libname, "libs-only-L").split()],
        "libraries": [word[2:] for word in run_pkg_config(libname, "libs-only-l").split()],
        "extra_compile_args": run_pkg_config(libname, "cflags-only-other").split(),
        "extra_link_args": run_pkg_config(libname, "libs-only-other").split(),
    }
    print_lib_config("From pkg-config", lib_config)
    return lib_config


def all_empty(lib_config):
    """Test if all lib_config fields are empty"""
    if len(lib_config) == 0:
        return True
    return all(len(value) == 0 for value in lib_config.values())


def all_exist(lib_config):
    """Test if all paths in the lib_config exist"""
    for key, value in lib_config.items():
        for path in value:
            if not os.path.exists(path):
                return False
    return True


def lib_config_magic(prefix, libname, static_config={}, known_include_dirs=[]):
    """Detect the configuration of a given library

    Parameters
    ----------

    prefix : str
        The prefix for this library. This is a name that HORTON uses to refer to the
        library.

    libname : str
        The library name as it is known to the compiler and to pkg-config. For example, if
        the shared object is libfoo.so, then the library name is foo.

    static_config : dict
        If given, this static library configuration is attempted. Ignored when empty, or
        when it contains non-existing files.

    known_include_dirs : list of str
        When all other methods of finding the library settings fail, the first existing
        directory in this list is added to the include path. This is useful when header
        files are commonly installed in a place that is not considered by default by most
        compilers.
    """
    print("%s Configuration" % prefix.upper())

    # Start out empty
    lib_config = dict((key, []) for key in lib_config_keys)

    # Update with info from setup.cfg
    lib_config.update(get_lib_config_setup(prefix, "setup.cfg"))

    # Override with environment variables
    lib_config.update(get_lib_config_env(prefix))

    # If no environment variables were set, attempt to use the static config.
    if all_empty(lib_config):
        if all_empty(static_config):
            print("   No static config available for this library")
        elif not all_exist(static_config):
            print_lib_config("Static lib not found in ${QAWORKDIR}", static_config)
        else:
            # If the static build is present, use it.
            print_lib_config("Static lib config in ${QAWORKDIR}", static_config)
            lib_config.update(static_config)

    # If also the static config did not work, try pkg-config
    if all_empty(lib_config):
        try:
            # Try to get dynamic link info from pkg-config
            lib_config.update(get_lib_config_pkg(libname))
        except PkgConfigError:
            print("   pkg-config failed.")

    # Uber-dumb fallback. It works most of the times.
    if all_empty(lib_config):
        lib_config["libraries"] = [libname]
        for include_dir in known_include_dirs:
            if os.path.isdir(include_dir):
                lib_config["include_dirs"] = [include_dir]
                break
        print_lib_config("Last resort fallback plan", lib_config)

    print_lib_config("Final", lib_config)
    return lib_config

# Locate ${QAWORKDIR}
# -------------------
qaworkdir = "tools/libs"

# Configuration of LibInt2
# ------------------------
# Load dependency information
with open("dependencies.json") as f:
    dependencies = json.load(f)
dependencies = dict((d['name'], d) for d in dependencies)
# Static build info in the QAWORKDIR:
libint2_dir = "%s/cached/libint-%s" % (qaworkdir, str(dependencies["libint"]["version_ci"]))
libint2_static_config = {
    "extra_objects": ["%s/lib/libint2.a" % libint2_dir],
    "include_dirs": ["%s/include/libint2" % libint2_dir],
}
# Common include dirs that are not considered by the compiler by default:
known_libint2_include_dirs = ["/usr/include/libint2", "/opt/local/include/libint2"]
libint2_config = lib_config_magic(
    "libint2", "int2", libint2_static_config, known_libint2_include_dirs)


setup(
    name="gbasis",
    version=get_version(),
    description="",
    author="Toon Verstraelen",
    author_email="Toon.Verstraelen@UGent.be",
    url="https://github.com/theochem/gbasis",
    package_dir={"gbasis": "gbasis"},
    packages=["gbasis", "gbasis.test"],
    cmdclass={"build_ext": Cython.Build.build_ext},
    ext_modules=[Extension(
        "gbasis.cext",
        sources=["gbasis/cext.pyx"] + glob.glob("gbasis/*.cpp"),
        depends=glob.glob("gbasis/*.h") + glob.glob("gbasis/*.h"),
        include_dirs=[np.get_include(), "."] +
                     libint2_config["include_dirs"],
        library_dirs=libint2_config["library_dirs"],
        libraries=libint2_config["libraries"],
        extra_objects=libint2_config["extra_objects"],
        extra_compile_args=libint2_config["extra_compile_args"] +
                           ["-std=c++11"],
        extra_link_args=libint2_config["extra_link_args"],
        language="c++"),
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
    requires=["numpy", "scipy", "setuptools", "distutils", "Cython"],
    )
