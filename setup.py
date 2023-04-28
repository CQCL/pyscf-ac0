#!/usr/bin/env python
# Copyright 2014-2020 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

NAME = "pyscf-ac0"
AUTHOR = "Michal Krompiec"
AUTHOR_EMAIL = "michal.krompiec@quantinuum.com"
DESCRIPTION = "Interface to and port of fragments of GAMMCOR code developed by Pernal group at Lodz University of Technology"
SO_EXTENSIONS = {}
DEPENDENCIES = ["pyscf", "numpy"]

#######################################################################
# Unless not working, nothing below needs to be changed.
import subprocess

metadata = globals()
import os
import sys
from setuptools import setup, find_namespace_packages, Extension

topdir = os.path.abspath(os.path.join(__file__, ".."))
modules = find_namespace_packages(include=["pyscf.*"])


def guess_version():
    for module in modules:
        module_path = os.path.join(topdir, *module.split("."))
        for version_file in ["__init__.py", "_version.py"]:
            version_file = os.path.join(module_path, version_file)
            if os.path.exists(version_file):
                with open(version_file, "r") as f:
                    for line in f.readlines():
                        if line.startswith("__version__"):
                            delim = '"' if '"' in line else "'"
                            return line.split(delim)[1]
    raise ValueError("Version string not found")


if not metadata.get("VERSION", None):
    VERSION = guess_version()

pyscf_lib_dir = os.path.join(topdir, "pyscf", "lib")

settings = {
    "name": metadata.get("NAME", None),
    "version": VERSION,
    "description": metadata.get("DESCRIPTION", None),
    "author": metadata.get("AUTHOR", None),
    "author_email": metadata.get("AUTHOR_EMAIL", None),
    "install_requires": metadata.get("DEPENDENCIES", []),
}

if sys.platform.startswith("darwin"):  # OSX
    if not "LDFLAGS" in os.environ:
        os.environ["LDFLAGS"] = ""
    os.environ["LDFLAGS"] = (
        os.environ["LDFLAGS"]
        + " -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"
    )
previous_path = os.getcwd()
path, _ = os.path.split(os.path.realpath(__file__))
os.chdir("pyscf/cas_ac0")
subprocess.run(
    # ["f2py", "--quiet", "-c", "-m", "pyscf/cas_ac0/ac0_lib", "pyscf/cas_ac0/accas_lib.f90"],
    ["f2py", "-c", "-m", "ac0_lib", "accas_lib.f90"],
    check=True,
    # stderr=subprocess.DEVNULL,
)
os.chdir(previous_path)

setup(
    include_package_data=True, 
    packages=modules, 
    entry_points={
        'console_scripts': ['rdm_ac0=pyscf.cas_ac0._cli:rdm_ac0']
    }, 
    **settings
)
