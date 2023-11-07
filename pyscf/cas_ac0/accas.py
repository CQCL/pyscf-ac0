# Copyright 2023 Quantinuum
#
# You may not use this file except in compliance with the Licence.
# You may obtain a copy of the Licence in the LICENCE file accompanying
# these documents
#

"""AC0 code, based on GAMMCOR https://github.com/pernalk/GAMMCOR."""

import numpy as np
import h5py
from pyscf import ao2mo
from pyscf.cas_ac0.ac0_lib import accas_lib as ac0

# disable 'Too many arguments'
# pylint: disable=R0913
# disable 'Too many local variables'
# pylint: disable=R0914
# disable C0103: Argument name ... doesn't conform to snake_case naming style
# pylint: disable=C0103


def _get_xone(hcore, mos):
    h1eff = np.einsum("ki,kl,lj->ij", mos, hcore, mos)
    return h1eff[np.tril_indices(h1eff.shape[0])]


def get_cas_ac0_energy(mf, mc):
    # nvirt = mc.mo_coeff.shape[0] - mc.ncore - mc.ncas
    # active_mask = [False] * mc.ncore + [True] * mc.ncas + [False] * nvirt
    _, dm2 = mc.fcisolver.make_rdm12(mc.ci, mc.ncas, mc.nelecas)
    # core_mask = [((not active_mask[i]) and o > 0) for i, o in enumerate(mf.mo_occ)]
    natural_orbitals = mc.mo_coeff
    nbasis = natural_orbitals.shape[1]
    if not mc.natorb:
        raise ValueError(
            "CAS-AC0: need natural orbitals, please set CASSCF.natorb=True"
        )
    integrals = ao2mo.outcore.full_iofree(mf.mol, natural_orbitals, aosym=1)
    integrals = integrals.reshape([nbasis] * 4)
    twono = ac0.get_two_el(integrals, ac0.get_two_el_size(nbasis), nbasis)

    cas_orbs = mc.ncas
    cas_elec = sum(mc.nelecas)

    nat_occ = mc.mo_occ
    occ = nat_occ / 2

    rdm2act = ac0.get_rdm2_act(ac0.getnrdm2act(cas_orbs), dm2, cas_orbs)  # NO basis

    xone = _get_xone(mf.get_hcore(), natural_orbitals)
    e_tot, _ = ac0.accas(
        mf.energy_nuc(),
        twono,
        np.eye(nbasis),  # UCAS is a unit matrix because everything is in NO basis
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


def get_cas_ac0_correction(mf, mc):
    # nvirt = mc.mo_coeff.shape[0] - mc.ncore - mc.ncas
    # active_mask = [False] * mc.ncore + [True] * mc.ncas + [False] * nvirt
    _, dm2 = mc.fcisolver.make_rdm12(mc.ci, mc.ncas, mc.nelecas)
    # core_mask = [((not active_mask[i]) and o > 0) for i, o in enumerate(mf.mo_occ)]
    natural_orbitals = mc.mo_coeff
    nbasis = natural_orbitals.shape[1]
    if not mc.natorb:
        raise ValueError(
            "CAS-AC0: need natural orbitals, please set CASSCF.natorb=True"
        )
    integrals = ao2mo.outcore.full_iofree(mf.mol, natural_orbitals, aosym=1)
    integrals = integrals.reshape([nbasis] * 4)
    twono = ac0.get_two_el(integrals, ac0.get_two_el_size(nbasis), nbasis)

    cas_orbs = mc.ncas
    cas_elec = sum(mc.nelecas)

    nat_occ = mc.mo_occ
    occ = nat_occ / 2

    rdm2act = ac0.get_rdm2_act(ac0.getnrdm2act(cas_orbs), dm2, cas_orbs)  # NO basis

    xone = _get_xone(mf.get_hcore(), natural_orbitals)
    _, e_corr = ac0.accas(
        mf.energy_nuc(),
        twono,
        np.eye(nbasis),  # UCAS is a unit matrix because everything is in NO basis
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
    return e_corr


def get_ac0_corr_energy_from_file(filename: str):
    data_file = h5py.File(
        filename,
        "r",
    )
    h_core = np.asarray(data_file["h_core"])
    e_nuc = np.asarray(data_file["e_nuc"])[0]
    natural_orbitals = np.asarray(data_file["natural_orbitals"])
    nat_occ = np.asarray(data_file["nat_occ"])
    ucas = np.asarray(data_file["ucas"])
    cas_dimensions = np.asarray(data_file["cas_dimensions"])
    cas_orbs = cas_dimensions[0]
    cas_elec = cas_dimensions[1]
    dm2 = np.asarray(data_file["dm2"])
    integrals = np.asarray(data_file["integrals"])

    nbasis = natural_orbitals.shape[1]
    rdm2_nat = ac0.trrdm2(dm2, ucas.T, ucas.shape[0])  # transform to NO basis

    integrals = integrals.reshape([nbasis] * 4)
    twono = ac0.get_two_el(integrals, ac0.get_two_el_size(nbasis), nbasis)
    occ = nat_occ / 2.0
    rdm2act = ac0.get_rdm2_act(ac0.getnrdm2act(cas_orbs), rdm2_nat, cas_orbs)

    xone = _get_xone(h_core, natural_orbitals)  # transform to NO basis

    _, e_corr = ac0.accas(
        e_nuc,
        twono,
        np.eye(nbasis),  # a unit matrix because everything is in NO basis
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
    return e_corr
