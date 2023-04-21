# Copyright 2023 Quantinuum
#
# You may not use this file except in compliance with the Licence.
# You may obtain a copy of the Licence in the LICENCE file accompanying
# these documents
#

"""AC0 code, based on GAMMCOR https://github.com/pernalk/GAMMCOR."""

import scipy
import numpy as np
from pyscf import ao2mo, lib

# disable 'Too many arguments'
# pylint: disable=R0913
# disable 'Too many local variables'
# pylint: disable=R0914
# disable redefined-outer-name (Redefining name 'ac0_not_found' from outer scope)
# pylint: disable=W0621
try:
    from pyscf.cas_ac0 import ac0_lib
    from ac0_lib import accas_lib as ac0_lib
except (ModuleNotFoundError, ImportError) as ac0_not_found:
    import subprocess
    import os

    previous_path = os.getcwd()
    path, filename = os.path.split(os.path.realpath(__file__))
    os.chdir(path)
    import platform

    if platform.system() == "Darwin":
        if not "LDFLAGS" in os.environ:
            os.environ["LDFLAGS"] = ""
        os.environ["LDFLAGS"] = (
            os.environ["LDFLAGS"]
            + " -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"
        )
    subprocess.run(
        ["f2py", "--quiet", "-c", "-m", "ac0_lib", "accas_lib.f90"],
        check=True,
        stderr=subprocess.DEVNULL,
    )
    try:
        # pylint: disable=C0412
        from pyscf.cas_ac0 import ac0_lib
        from ac0_lib import accas_lib as ac0_lib
    except (ModuleNotFoundError, ImportError) as ac0_not_found:
        raise RuntimeError(
            "Cannot import ac0_lib, most likely due to an f2py compilation failure. "
            + "If on Mac, make sure "
            + "-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib is in the LDFLAGS"
        ) from ac0_not_found

    os.chdir(previous_path)


def _get_no_trafo_mask(act_rdm1, orbitals, core_mask, active_mask):
    cas_occ, ucas = scipy.linalg.eigh(-act_rdm1)
    idx = np.argmax(abs(ucas.real), axis=0)
    ucas[:, ucas[idx, np.arange(len(cas_occ))].real < 0] *= -1
    casorb_idx = np.argsort(cas_occ.round(9), kind="mergesort")
    if not np.all(casorb_idx[:-1] <= casorb_idx[1:]):
        cas_occ = cas_occ[casorb_idx]
        ucas = ucas[:, casorb_idx]
    cas_occ = -cas_occ
    mo_occ = np.zeros(orbitals.shape[1])
    mo_occ[core_mask] = 2
    mo_occ[active_mask] = cas_occ

    natorbs = orbitals.copy()
    natorbs[:, active_mask] = orbitals[:, active_mask] @ ucas
    if getattr(orbitals, "orbsym", None) is not None:
        orbsym = np.copy(orbitals.orbsym)
        orbsym[active_mask] = orbsym[active_mask][casorb_idx]
        natorbs = lib.tag_array(natorbs, orbsym=orbsym)

    return natorbs, mo_occ, cas_occ, ucas


def _get_xone(hcore, mos):
    h1eff = np.einsum("ki,kl,lj->ij", mos, hcore, mos)
    return h1eff[np.tril_indices(h1eff.shape[0])]


def get_cas_ac0_energy(mf, mc):
    nvirt = mc.mo_coeff.shape[0] - mc.ncore - mycas.ncas
    active_mask = [False] * mc.ncore + [True] * mc.ncas + [False] * nvirt
    dm1, dm2 = mc.fcisolver.make_rdm12(mc.ci, mc.ncas, mc.nelecas)
    core_mask = [((not active_mask[i]) and o > 0) for i, o in enumerate(mf.mo_occ)]
    if mc.natorb:
        natural_orbitals = mc.mo_coeff
    else:
        natural_orbitals, _, _, _ = _get_no_trafo_mask(
            dm1, mc.mo_coeff, core_mask, active_mask
        )
    raw_integrals = ao2mo.outcore.full_iofree(
        mf.mol, natural_orbitals, aosym=1
    )  # NO basis

    return get_ac0_energy(
        h_core=mf.get_hcore(),
        e_nuc=mf.energy_nuc(),
        orbitals=mc.mo_coeff,
        cas_orbs=mc.ncas,
        cas_elec=sum(mc.nelecas),
        core_mask=core_mask,
        active_mask=active_mask,
        dm1=dm1,
        dm2=dm2,
        integrals=raw_integrals,
    )


def get_ac0_energy(
    h_core,
    e_nuc,
    orbitals,
    cas_orbs,
    cas_elec,
    core_mask,
    active_mask,
    dm1,
    dm2,
    integrals,
):
    nbasis = orbitals.shape[0]
    natorbs, nat_occ, _, ucas = _get_no_trafo_mask(
        dm1, orbitals, core_mask, active_mask
    )

    rdm2_nat = ac0_lib.trrdm2(dm2, ucas.T, ucas.shape[0])  # transform to NO basis

    integrals = integrals.reshape([nbasis] * 4)
    twono = ac0_lib.get_two_el(
        integrals, ac0_lib.get_two_el_size(nbasis), nbasis
    )  # NO basis
    occ = nat_occ / 2

    rdm2act = ac0_lib.get_rdm2_act(
        ac0_lib.getnrdm2act(cas_orbs), rdm2_nat, cas_orbs
    )  # NO basis

    xone = _get_xone(h_core, natorbs)  # NO basis

    e_tot = ac0_lib.accas(
        e_nuc,
        twono,
        np.eye(nbasis),  # still a unit matrix because everything is in NO basis
        occ,
        xone,
        rdm2act,
        int(round(np.sum(occ) * 2)),
        cas_orbs,
        cas_elec,
        nbasis,
        xone.shape[0],
        twono.shape[0],
    )
    return e_tot


if __name__ == "__main__":
    import pyscf

    mol = pyscf.M(atom="O 0 0 0; O 0 0 1.2", basis="ccpvdz", spin=2)
    myhf = mol.RHF().run()
    ncas, nelecas = (6, (5, 3))
    mycas = myhf.CASSCF(ncas, nelecas).run()
    print(f"CAS-AC0: {get_cas_ac0_energy(myhf, mycas)}")
