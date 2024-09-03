from pyscf import lib, ao2mo
import scipy.linalg
import numpy as np
import itertools
import functools

einsum = functools.partial(np.einsum, optimize=True)


def eigxy(apb, amb):
    """Solve the RPA eigenproblem.
    """
    w, v = np.linalg.eigh(apb)
    apb_sqrt = (v * np.sqrt(w)[None]) @ v.T

    w, v = np.linalg.eigh(apb_sqrt @ amb @ apb_sqrt)
    mask = w > 1e-10
    w = np.sqrt(w[mask])
    v = v[:, mask]

    y = apb_sqrt @ v
    norm = einsum("ij,in,jn->n", amb, y, y)
    norm = 2 * norm / w
    y /= np.sqrt(norm)[None]

    x = einsum("ij,jx,x->ix", amb, y, 1 / w)

    return w, x, y


def _get_veff(mc):
    """Get a the blocks of the effective potential.

    Returns a list of the blocks of the effective potential, where the blocks correspond to summing
    over the correlated, active, and virtual density matrices, respectively.

    Args:
        mc: The CASSCF object.

    Returns:
        The blocks of the effective potential.
    """
    # Get the spaces
    ncor = mc.ncore
    nact = mc.ncas
    nvir = mc.mo_occ.size - ncor - nact
    cor = slice(0, ncor)
    act = slice(ncor, ncor + nact)
    vir = slice(ncor + nact, None)

    # Get the density matrices
    dm = np.array([
        einsum("i,pi,qi->pq", mc.mo_occ[slc], mc.mo_coeff[:, slc], mc.mo_coeff[:, slc])
        for slc in (cor, act, vir)
    ])

    # Evaluate the effective potential in the AO basis
    veff = mc._scf.get_veff(mc._scf._eri, dm)

    # Rotate the effective potential to the NO basis
    veff = einsum("...pq,pi,qj->...ij", veff, mc.mo_coeff, mc.mo_coeff)

    return list(veff)


def _build_act_act_0(mc, h1e, rdm2):
    """Build the active-active block of the α=0 A+B and A-B matrices.

    Args:
        mc: The CASSCF object.
        h1e: The 1-electron α=0 Hamiltonian in a basis of natural orbitals.
        rdm2: The 2-particle CAS reduced density matrix in a basis of natural orbitals, for the
            active space.
    """
    # Get the spaces
    ncor = mc.ncore
    nact = mc.ncas
    act = slice(ncor, ncor + nact)
    mo_coeff_act = mc.mo_coeff[:, act]

    # Get the occupancy measures
    mo_occ = mc.mo_occ[act] * 0.5
    mo_occ_sqrt = np.sqrt(mo_occ)
    mo_occ_sqrt[mo_occ < 0.5] *= -1
    mo_occ_difs = np.subtract.outer(mo_occ, mo_occ)

    # Get the Hamiltonian blocks
    h2e_aaaa = ao2mo.kernel(mc._scf._eri, mo_coeff_act, compact=False).reshape((nact,) * 4)

    # One-body terms
    a = einsum("ps,qs,pr->pqrs", mo_occ_difs, h1e[act, act], np.eye(nact))

    # Two-body terms
    a += einsum("sqtu,purt->pqrs", h2e_aaaa, rdm2)
    a += einsum("sutq,putr->pqrs", h2e_aaaa, rdm2)
    a -= einsum("ptsu,tuqr->pqrs", h2e_aaaa, rdm2)

    # Two-body terms with identity
    tmp = einsum("twpu,wutr->pr", h2e_aaaa, rdm2)
    tmp += einsum("tupw,wurt->pr", h2e_aaaa, rdm2)
    a -= einsum("pr,qs->pqrs", tmp, np.eye(nact)) * 0.5

    # Symmetrise
    a += a.transpose(1, 0, 3, 2)

    # Get the index helpers
    pairs = np.abs(np.subtract.outer(mo_occ, mo_occ)) > 1e-8
    p, q = np.tril_indices(nact)
    mask = pairs[p, q]
    p, q = p[mask], q[mask]
    tril = np.tril_indices(p.size)
    triu = np.triu_indices(p.size)

    # Pack the A+B matrix
    apb = np.zeros((p.size, p.size))
    apb[triu] += a[p, q][:, p, q][triu]
    apb[triu] += a[q, p][:, p, q][triu]
    apb[tril] = apb.T[tril]
    m = np.add.outer(mo_occ_sqrt, mo_occ_sqrt)[p, q]
    m[m == 0] = 1
    apb /= np.multiply.outer(m, m)

    # Pack the A-B matrix
    amb = np.zeros((p.size, p.size))
    amb[triu] += a[p, q][:, p, q][triu]
    amb[triu] -= a[q, p][:, p, q][triu]
    amb[tril] = amb.T[tril]
    m = np.subtract.outer(mo_occ_sqrt, mo_occ_sqrt)[p, q]
    m[m == 0] = 1
    amb /= np.multiply.outer(m, m)

    return apb, amb


def _build_act_cor_0(mc, h1e, rdm2, veffs):
    """Build the active-core block of the α=0 A+B and A-B matrices.

    Args:
        mc: The CASSCF object.
        h1e: The 1-electron α=0 Hamiltonian in a basis of natural orbitals.
        rdm2: The 2-particle CAS reduced density matrix in a basis of natural orbitals, for the
            active space.
        veffs: The contributions to the effective potential.
    """
    # Get the spaces
    ncor = mc.ncore
    nact = mc.ncas
    cor = slice(0, ncor)
    act = slice(ncor, ncor + nact)
    occ = slice(0, ncor + nact)
    mo_coeff_act = mc.mo_coeff[:, act]
    mo_coeff_cor = mc.mo_coeff[:, cor]

    # Get the density matrices and occupancy measures
    mo_occ = mc.mo_occ * 0.5
    mo_occ_sqrt = np.sqrt(mo_occ)
    mo_occ_sqrt[mo_occ < 0.5] *= -1
    mo_occ_difs = np.subtract.outer(mo_occ, mo_occ)

    # Get the Hamiltonian blocks
    h2e_cccc = ao2mo.kernel(mc._scf._eri, mo_coeff_cor, compact=False).reshape((ncor,) * 4)
    h2e_aaaa = ao2mo.kernel(mc._scf._eri, mo_coeff_act, compact=False).reshape((nact,) * 4)

    # Note: We only need A_{pqrq}, so re-order as A_{qpr} for efficiency

    # One-body terms
    a = einsum("pq,qq,pr->qpr", mo_occ_difs[act, cor], h1e[cor, cor], np.eye(nact))
    a += einsum("qr,pr,qq->qpr", mo_occ_difs[cor, act], h1e[act, act], np.eye(ncor))

    # One-body terms with V_{eff}
    tmp = veffs[0][cor, cor]
    a += einsum("pr,qq->qpr", np.diag(mo_occ[act]), tmp)

    # Two-body terms
    a += einsum("ttpr,q,t->qpr", h2e_aaaa, mo_occ[cor], mo_occ[act]) * 2
    a -= einsum("trpt,q,t->qpr", h2e_aaaa, mo_occ[cor], mo_occ[act])

    # Two-body terms with identity
    tmp = einsum("twpu,wutr->pr", h2e_aaaa, rdm2)
    tmp += einsum("tupw,wurt->pr", h2e_aaaa, rdm2)
    a -= einsum("pr,qq->qpr", tmp, np.eye(ncor)) * 0.5
    tmp = einsum("ttpr,t,r->pr", h2e_cccc, mo_occ[cor], mo_occ[cor]) * 4
    tmp -= einsum("trpt,r,t->pr", h2e_cccc, mo_occ[cor], mo_occ[cor]) * 2
    a -= einsum("pr,qq->qpr", np.eye(nact), tmp) * 0.5

    apb = []
    amb = []
    for q in range(ncor):
        # Get the index helpers
        mask = np.abs(mo_occ[act] - mo_occ[cor][q]) > 1e-8
        if not np.any(mask):
            continue
        p = np.arange(nact)[mask]
        tril = np.tril_indices(p.size)
        triu = np.triu_indices(p.size)

        # Pack the A+B matrix
        apb_q = np.zeros((p.size, p.size))
        apb_q[triu] = a[q].reshape(p.size, p.size)[triu]
        apb_q[tril] = a[q].reshape(p.size, p.size).T[tril]
        m = mo_occ_sqrt[act] + mo_occ_sqrt[cor][q]
        m[m == 0] = 1
        apb_q /= np.multiply.outer(m, m)
        apb.append(apb_q)

        # Pack the A-B matrix
        amb_q = np.zeros((p.size, p.size))
        amb_q[triu] = a[q].reshape(p.size, p.size)[triu]
        amb_q[tril] = a[q].reshape(p.size, p.size).T[tril]
        m = mo_occ_sqrt[act] - mo_occ_sqrt[cor][q]
        m[m == 0] = 1
        amb_q /= np.multiply.outer(m, m)
        amb.append(amb_q)

    return apb, amb


def _build_vir_act_0(mc, h1e, rdm2):
    """Build the virtual-active block of the α=0 A+B and A-B matrices.

    Args:
        mc: The CASSCF object.
        h1e: The 1-electron α=0 Hamiltonian in a basis of natural orbitals.
        rdm2: The 2-particle CAS reduced density matrix in a basis of natural orbitals, for the
            active space.
    """
    # Get the spaces
    ncor = mc.ncore
    nact = mc.ncas
    nvir = mc.mo_occ.size - ncor - nact
    act = slice(ncor, ncor + nact)
    vir = slice(ncor + nact, None)
    mo_coeff_act = mc.mo_coeff[:, act]

    # Get the occupancy measures
    mo_occ = mc.mo_occ * 0.5
    mo_occ_sqrt = np.sqrt(mo_occ)
    mo_occ_sqrt[mo_occ < 0.5] *= -1
    mo_occ_difs = np.subtract.outer(mo_occ, mo_occ)

    # Get the Hamiltonian blocks
    h2e_aaaa = ao2mo.kernel(mc._scf._eri, mo_coeff_act, compact=False).reshape((nact,) * 4)

    # Note: We only need A_{pqps}, so re-order as A_{pqs} for efficiency

    # One-body terms
    a = einsum("ps,qs,pp->pqs", mo_occ_difs[vir, act], h1e[act, act], np.eye(nvir))
    a += einsum("qp,pp,sq->pqs", mo_occ_difs[act, vir], h1e[vir, vir], np.eye(nact))

    # Two-body terms with identity
    tmp = einsum("twpu,wutr->pr", h2e_aaaa, rdm2)
    tmp += einsum("tupw,wurt->pr", h2e_aaaa, rdm2)
    a -= einsum("pp,qs->pqs", np.eye(nvir), tmp) * 0.5

    apb = []
    amb = []
    for p in range(nvir):
        # Get the index helpers
        mask = np.abs(mo_occ[vir][p] - mo_occ[act]) > 1e-8
        if not np.any(mask):
            continue
        q = np.arange(nact)[mask]
        tril = np.tril_indices(q.size)
        triu = np.triu_indices(q.size)

        # Pack the A+B matrix
        apb_p = np.zeros((q.size, q.size))
        apb_p[triu] = a[p].reshape(q.size, q.size)[triu]
        apb_p[tril] = a[p].reshape(q.size, q.size).T[tril]
        m = mo_occ_sqrt[vir][p] + mo_occ_sqrt[act]
        m[m == 0] = 1
        apb_p /= np.multiply.outer(m, m)
        apb.append(apb_p)

        # Pack the A-B matrix
        amb_p = np.zeros((q.size, q.size))
        amb_p[triu] = a[p].reshape(q.size, q.size)[triu]
        amb_p[tril] = a[p].reshape(q.size, q.size).T[tril]
        m = mo_occ_sqrt[vir][p] - mo_occ_sqrt[act]
        m[m == 0] = 1
        amb_p /= np.multiply.outer(m, m)
        amb.append(amb_p)

    return apb, amb


def _build_vir_cor_0(mc, h1e, rdm2):
    """Build the virtual-core block of the α=0 A+B and A-B matrices.

    Args:
        mc: The CASSCF object.
        h1e: The 1-electron α=0 Hamiltonian in a basis of natural orbitals.
        rdm2: The 2-particle CAS reduced density matrix in a basis of natural orbitals, for the
            active space.
    """
    # Get the spaces
    ncor = mc.ncore
    nact = mc.ncas
    nvir = mc.mo_occ.size - ncor - nact
    cor = slice(0, ncor)
    vir = slice(ncor + nact, None)
    mo_coeff_cor = mc.mo_coeff[:, cor]

    # Get the occupancy measures
    mo_occ = mc.mo_occ * 0.5
    mo_occ_difs = np.subtract.outer(mo_occ, mo_occ)

    # Get the Hamiltonian blocks
    h2e_cccc = ao2mo.kernel(mc._scf._eri, mo_coeff_cor, compact=False).reshape((ncor,) * 4)

    # Note: We only need A_{pqpq}, so re-order as A_{pq} for efficiency

    # One-body terms
    a = einsum("pq,qq->pq", mo_occ_difs[vir, cor], h1e[cor, cor])
    a += einsum("qp,pp->pq", mo_occ_difs[cor, vir], h1e[vir, vir])

    # Two-body terms with identity
    tmp = einsum("ttpr,t,r->pr", h2e_cccc, mo_occ[cor], mo_occ[cor]) * 4
    tmp -= einsum("trpt,r,t->pr", h2e_cccc, mo_occ[cor], mo_occ[cor]) * 2
    a -= einsum("pp,qq->pq", np.eye(nvir), tmp) * 0.5

    # Pack the A+B and A-B matrices
    apb = amb = [np.array([[x]]) for x in a.ravel()]

    return apb, amb


def _calculate_energy(mc, h1e, rdm2, veffs, w_0, x_0, y_0, block_size=16):
    """Build a block of the α=1 A+B and A-B matrices.

    Args:
        mc: The CASSCF object.
        h1e: The 1-electron α=1 Hamiltonian in a basis of natural orbitals.
        rdm2: The 2-particle CAS reduced density matrix in a basis of natural orbitals, for the
            active space.
        veffs: The contributions to the effective potential.
        w_0: The eigenvalues of the α=0 A+B and A-B matrices.
        x_0: The X eigenvectors of the α=0 A+B and A-B matrices.
        y_0: The Y eigenvectors of the α=0 A+B and A-B matrices.
        block_size: The size of the blocks to iterate over. Lower sizes will use less memory but
            may be slower.
    """
    # Get the spaces
    ncor = mc.ncore
    nact = mc.ncas
    nvir = mc.mo_occ.size - ncor - nact
    norb = mc.mo_occ.size
    cor = slice(0, ncor)
    act = slice(ncor, ncor + nact)
    vir = slice(ncor + nact, norb)
    not_cor = slice(ncor, norb)
    not_vir = slice(0, ncor + nact)

    # Get the occupancy measures
    mo_occ = mc.mo_occ * 0.5
    mo_occ_sqrt = np.sqrt(mo_occ)
    mo_occ_sqrt[mo_occ < 0.5] *= -1
    mo_occ_difs = np.subtract.outer(mo_occ, mo_occ)

    # Get the index helpers
    pairs = np.abs(mo_occ_difs) > 1e-8
    pq = sorted((p, q) for q, p in itertools.combinations(range(act.start, act.stop), 2))
    pq += [x[::-1] for x in itertools.product(range(cor.start, cor.stop), range(act.start, act.stop))]
    pq += list(itertools.product(range(vir.start, vir.stop), range(act.start, act.stop)))
    pq += list(itertools.product(range(vir.start, vir.stop), range(cor.start, cor.stop)))
    pq = [(p, q) for p, q in pq if pairs[p, q]]
    p = np.array([x[0] for x in pq])
    q = np.array([x[1] for x in pq])
    tril = np.tril_indices(p.size)
    triu = np.triu_indices(p.size)

    def _get_contribution(p):
        # Get the Hamiltonian blocks
        h2e = ao2mo.kernel(
            mc._scf._eri,
            (mc.mo_coeff[:, p], mc.mo_coeff, mc.mo_coeff, mc.mo_coeff),
            compact=False,
        ).reshape((p.stop - p.start, norb, norb, norb))

        a = einsum("ps,qs,pr->pqrs", mo_occ_difs[p], h1e, np.eye(norb)[p])

        tmp = einsum("pqrs,q,s->pqrs", h2e[:, cor, :, not_cor], mo_occ[cor], mo_occ[not_cor]) * 2
        tmp -= einsum("prqs,q,s->pqrs", h2e[:, :, cor, not_cor], mo_occ[cor], mo_occ[not_cor])
        a[:, cor, :, not_cor] += tmp

        tmp = einsum("pqrs,q,s->pqrs", h2e[:, not_cor, :, cor], mo_occ[not_cor], mo_occ[cor]) * 2
        tmp -= einsum("prqs,q,s->pqrs", h2e[:, :, not_cor, cor], mo_occ[not_cor], mo_occ[cor])
        a[:, not_cor, :, cor] += tmp

        tmp = einsum("pqrs,q,r->pqrs", h2e[:, cor, not_cor, :], mo_occ[cor], mo_occ[not_cor]) * 2
        tmp -= einsum("prsq,q,r->pqrs", h2e[:, not_cor, :, cor], mo_occ[cor], mo_occ[not_cor])
        a[:, cor, not_cor, :] -= tmp

        tmp = einsum("pqrs,q,r->pqrs", h2e[:, not_cor, cor, :], mo_occ[not_cor], mo_occ[cor]) * 2
        tmp -= einsum("prsq,q,r->pqrs", h2e[:, cor, :, not_cor], mo_occ[not_cor], mo_occ[cor])
        a[:, not_cor, cor, :] -= tmp

        a[:, cor, :, cor] += einsum("pqrs,q,s->pqrs", h2e[:, cor, :, cor], mo_occ[cor], mo_occ[cor]) * 2
        a[:, cor, :, cor] -= einsum("prqs,q,s->pqrs", h2e[:, :, cor, cor], mo_occ[cor], mo_occ[cor])

        a[:, cor, cor, :] -= einsum("pqrs,q,r->pqrs", h2e[:, cor, cor, :], mo_occ[cor], mo_occ[cor]) * 2
        a[:, cor, cor, :] += einsum("prsq,q,r->pqrs", h2e[:, cor, :, cor], mo_occ[cor], mo_occ[cor])

        tmp = veffs[0][p, :] + veffs[1][p, :]
        a[:, cor, :, cor] += einsum("pr,qs->pqrs", tmp, np.diag(mo_occ[cor]))

        tmp = veffs[0][p, :]
        a[:, not_cor, :, not_cor] += einsum("pr,qs->pqrs", tmp, np.diag(mo_occ[not_cor]))

        a[:, act, :, act] += einsum("prtu,qust->pqrs", h2e[:, :, act, act], rdm2)
        a[:, act, :, act] += einsum("ptur,quts->pqrs", h2e[:, act, act, :], rdm2)
        a[:, act, act, :] -= einsum("ptus,tuqr->pqrs", h2e[:, act, act, :], rdm2)

        tmp = np.zeros((p.stop - p.start, norb))
        tmp[:, act] += einsum("putw,wutr->pr", h2e[:, act, act, act], rdm2)
        tmp[:, act] += einsum("pwtu,wurt->pr", h2e[:, act, act, act], rdm2)
        for slc, sign in ((not_vir, 1), (act, -1)):
            # Avoids double counting of the active part
            tmp[:, slc] += einsum("prtt,t,r->pr", h2e[:, slc, slc, slc], mo_occ[slc], mo_occ[slc]) * 2 * sign
            tmp[:, slc] -= einsum("pttr,r,t->pr", h2e[:, slc, slc, slc], mo_occ[slc], mo_occ[slc]) * sign
            tmp[:, slc] += einsum("prtt,r,t->pr", h2e[:, slc, slc, slc], mo_occ[slc], mo_occ[slc]) * 2 * sign
            tmp[:, slc] -= einsum("pttr,t,r->pr", h2e[:, slc, slc, slc], mo_occ[slc], mo_occ[slc]) * sign
        a -= einsum("pr,qs->pqrs", tmp, np.eye(norb)) * 0.5

        return a

    apb = np.zeros((p.size, p.size))
    amb = np.zeros((p.size, p.size))
    for p0 in range(0, norb, block_size):
        p1 = min(p0 + block_size, norb)
        a = _get_contribution(slice(p0, p1))

        # For (p, q) contributions
        imask = np.logical_and(p >= p0, p < p1)
        pi = p[imask] - p0
        qi = q[imask]

        # For (q, p) contributions
        jmask = np.logical_and(q >= p0, q < p1)
        pj = p[jmask]
        qj = q[jmask] - p0

        apb[imask] += a[pi, qi][:, p, q]
        apb[jmask] += a[qj, pj][:, p, q]
        apb[jmask] += a[qj, pj][:, q, p]
        apb[imask] += a[pi, qi][:, q, p]

        amb[imask] += a[pi, qi][:, p, q]
        amb[jmask] -= a[qj, pj][:, p, q]
        amb[jmask] += a[qj, pj][:, q, p]
        amb[imask] -= a[pi, qi][:, q, p]

    apb[tril] = apb.T[tril]
    m = np.add.outer(mo_occ_sqrt, mo_occ_sqrt)[p, q]
    m[m == 0] = 1
    apb /= np.multiply.outer(m, m)

    amb[tril] = amb.T[tril]
    m = np.subtract.outer(mo_occ_sqrt, mo_occ_sqrt)[p, q]
    m[m == 0] = 1
    amb /= np.multiply.outer(m, m)

    # Rotate the A+B and A-B matrices into the X and Y basis
    apb = x_0.T @ apb @ x_0
    amb = y_0.T @ amb @ y_0

    # Find the perturbative term
    xm = 2 * (apb - amb) / np.add.outer(w_0, w_0)
    xm = y_0 @ xm @ y_0.T
    m = np.add.outer(mo_occ_sqrt, mo_occ_sqrt)[p, q]
    xm *= np.multiply.outer(m, m)

    # Calculate the energy
    e_corr = 0.0
    block = scipy.linalg.block_diag(*[np.ones((n, n)) for n in [ncor, nact, nvir]])
    for p0 in range(0, norb, block_size):
        p1 = min(p0 + block_size, norb)

        h2e = ao2mo.kernel(
            mc._scf._eri,
            (mc.mo_coeff[:, p0:p1], mc.mo_coeff, mc.mo_coeff, mc.mo_coeff),
            compact=False,
        ).reshape((p1 - p0, norb, norb, norb))
        for slc in (cor, act, vir):
            h2e[max(slc.start, p0) : min(slc.stop, p1), slc, slc, slc] = 0.0

        mask = np.logical_and(p >= p0, p < p1)
        pi = p[mask] - p0
        qi = q[mask]

        e_corr += np.sum(xm[mask] * h2e[pi, qi][:, p, q])

    return e_corr


def run(mf, mc, block_size=16):
    nact = mc.ncas
    nelec = mc.nelecas

    rdm2 = mc.fcisolver.make_rdm2(mc.ci, mc.ncas, mc.nelecas)

    norb = mf.mo_occ.size
    ncor = mc.ncore
    nact = mc.ncas
    nvir = norb - ncor - nact

    cor = slice(0, ncor)
    act = slice(ncor, ncor + nact)
    vir = slice(ncor + nact, norb)

    # Get the 1-electron Hamiltonian
    hcore = einsum("pq,pi,qj->ij", mf.get_hcore(), mc.mo_coeff, mc.mo_coeff)
    veffs = _get_veff(mc)
    fock = hcore + veffs[0] + veffs[1] + veffs[2]
    rdm2 = rdm2.transpose(1, 3, 0, 2) * 0.5
    block = scipy.linalg.block_diag(*[np.ones((n, n)) for n in [ncor, nact, nvir]])

    # Get the α=0 1e Hamiltonian
    h1e = fock * block
    for i, slc in enumerate((cor, act, vir)):
        h1e[slc, slc] -= veffs[i][slc, slc]

    # Active-active
    apb_act_act_0, amb_act_act_0 = _build_act_act_0(mc, h1e, rdm2)
    w_act_act_0, x_act_act_0, y_act_act_0 = eigxy(apb_act_act_0, amb_act_act_0)

    # Active-core
    apb_act_cor_0, amb_act_cor_0 = _build_act_cor_0(mc, h1e, rdm2, veffs)
    if apb_act_cor_0:
        w_act_cor_0, x_act_cor_0, y_act_cor_0 = zip(*[eigxy(apb, amb) for apb, amb in zip(apb_act_cor_0, amb_act_cor_0)])
    else:
        w_act_cor_0, x_act_cor_0, y_act_cor_0 = [], [], []

    # Virtual-active
    apb_vir_act_0, amb_vir_act_0 = _build_vir_act_0(mc, h1e, rdm2)
    if apb_vir_act_0:
        w_vir_act_0, x_vir_act_0, y_vir_act_0 = zip(*[eigxy(apb, amb) for apb, amb in zip(apb_vir_act_0, amb_vir_act_0)])
    else:
        w_vir_act_0, x_vir_act_0, y_vir_act_0 = [], [], []

    # Virtual-core
    apb_vir_cor_0, amb_vir_cor_0 = _build_vir_cor_0(mc, h1e, rdm2)
    w_vir_cor_0 = [np.array(a.ravel()) for a in apb_vir_cor_0]
    x_vir_cor_0 = [np.ones_like(a) / np.sqrt(2.0) for a in apb_vir_cor_0]
    y_vir_cor_0 = [np.ones_like(a) / np.sqrt(2.0) for a in apb_vir_cor_0]

    # Collect
    w_0 = np.concatenate([w_act_act_0, *w_act_cor_0, *w_vir_act_0, *w_vir_cor_0])
    x_0 = scipy.linalg.block_diag(x_act_act_0, *x_act_cor_0, *x_vir_act_0, *x_vir_cor_0)
    y_0 = scipy.linalg.block_diag(y_act_act_0, *y_act_cor_0, *y_vir_act_0, *y_vir_cor_0)

    # Get the α=1 1e Hamiltonian
    h1e = hcore * (1 - block)
    for eslc in (cor, act, vir):
        for i, islc in enumerate((cor, act, vir)):
            if eslc != islc:
                h1e[eslc, eslc] -= veffs[i][eslc, eslc]

    # Calculate the energy
    e_corr = _calculate_energy(mc, h1e, rdm2, veffs, w_0, x_0, y_0, block_size=block_size)

    return e_corr


if __name__ == "__main__":
    import time
    from memory_profiler import memory_usage
    from pyscf import gto, scf, mcscf
    from pyscf.cas_ac0.accas import get_cas_ac0_energy

    print(f"               | {'time (ms)':^17s} | {'memory (ms)':^17s}")
    print(f"ncor nact nvir | {'old':>8s} {'new':>8s} | {'old':>8s} {'new':>8s}")
    for mol, cas in [
        (gto.M(atom="O 0 0 0; H 0 0 1; H 0 1 0", basis="6-31g", verbose=0), (6, 4)),
        (gto.M(atom="Li 0 0 0; H 0 0 2", basis="sto3g", verbose=0), (4, 2)),
        (gto.M(atom="Li 0 0 0; H 0 0 2", basis="sto3g", verbose=0), (4, 4)),
        (gto.M(atom="Li 0 0 0; H 0 0 2", basis="cc-pvdz", verbose=0), (4, 2)),
        (gto.M(atom="Li 0 0 0; Li 0 0 2", basis="6-31g", verbose=0), (6, 6)),
        (gto.M(atom="O 0 0 0; H 0 0 1; H 0 1 0", basis="sto3g", verbose=0), (2, 0)),
        (gto.M(atom="O 0 0 0; H 0 0 1; H 0 1 0", basis="cc-pvdz", verbose=0), (2, 0)),
        (gto.M(atom="O 0 0 0; O 0 0 1", basis="aug-cc-pvdz", verbose=0), (2, 0)),
        #(gto.M(atom="O 0 0 0; O 0 0 1; O 0 0 2", basis="aug-cc-pvdz", verbose=0), (6, 6)),
    ]:
        mol.max_memory = 1e10

        mf = scf.RHF(mol)
        mf.conv_tol = 1e-14
        mf.kernel()

        mc = mcscf.CASSCF(mf, *cas)
        mc.fcisolver.conv_tol = 1e-10
        mc.natorb = True
        mc.kernel()

        f1 = lambda: run(mf, mc) + mc.e_tot
        t0 = time.time()
        e1 = f1()
        t1 = time.time() - t0
        m1 = max(memory_usage(f1, interval=1e-4))

        f2 = lambda: get_cas_ac0_energy(mf, mc)
        t0 = time.time()
        e2 = f2()
        t2 = time.time() - t0
        m2 = max(memory_usage(f2, interval=1e-4))

        assert abs(e1 - e2) < 1e-10
        print(f"{mc.ncore:4d} {mc.ncas:4d} {mol.nao-mc.ncore-mc.ncas:4d} | {t2 * 1000:8.2f} {t1 * 1000:8.2f} | {m2:8.1f} {m1:8.1f}", flush=True)
