from numpy.distutils.core import setup, Extension
import os
import sys

if sys.platform.startswith("darwin"):  # OSX
    if not "LDFLAGS" in os.environ:
        os.environ["LDFLAGS"] = ""
    os.environ["LDFLAGS"] = (
        os.environ["LDFLAGS"]
        + " -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"
    )

setup_args = dict(
    ext_modules=[
        Extension(
            name='pyscf.cas_ac0.ac0_lib',
            sources=['pyscf/cas_ac0/accas_lib.f90'],
        )
    ],
)
setup(**setup_args)
