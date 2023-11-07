import pyscf
from pyscf.cas_ac0 import accas
import numpy as np
import pytest

def test_accas_O2():
    mol = pyscf.M(atom="O 0 0 0; O 0 0 1.2", basis="ccpvdz", spin=2)
    myhf = mol.RHF().run()
    ncas, nelecas = (6, (5, 3))
    mycas = myhf.CASSCF(ncas, nelecas)
    mycas.natorb = True
    mycas.run()
    e = accas.get_cas_ac0_energy(myhf, mycas)
    # print(f"CAS-AC0: {e}")
    assert e == pytest.approx(-149.96063718895718)


def test_accas_H2():
    mol = pyscf.M(atom="H 0 0 0; H 0.7 0 0", basis="cc-pvdz")
    myhf = mol.RHF().run()
    ncas, nelecas = (2, 2)
    mycas = myhf.CASSCF(ncas, nelecas)
    mycas.natorb = True
    mycas.run()
    e = accas.get_cas_ac0_energy(myhf, mycas)
    # print(f"CAS-AC0: {e}")
    assert e == pytest.approx(-1.1572426436084569, abs=1e-8)
