#!/usr/bin/env python3

from pyscf.cas_ac0 import accas
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: rdm_ac0 file.h5")
    else:
        print(accas.get_ac0_corr_energy_from_file(sys.argv[1]))
