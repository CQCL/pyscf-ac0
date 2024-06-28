# Copyright 2023 Quantinuum
#
# You may not use this file except in compliance with the Licence.
# You may obtain a copy of the Licence in the LICENCE file accompanying
# these documents
#

"""Module for defining command-line interface to cas_ac0 functionality."""

from pyscf.cas_ac0.accas import get_ac0_corr_energy_from_file
import sys


def rdm_ac0():
    """Command-line interface to compute ac0 core energy from input h5 file."""
    if len(sys.argv) != 2:
        raise ValueError("Usage: rdm_ac0 file.h5")
    else:
        print(get_ac0_corr_energy_from_file(sys.argv[1]))
