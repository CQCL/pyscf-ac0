from pyscf import lib, ao2mo
import scipy.linalg
import numpy as np
import itertools


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
    norm = lib.einsum("ij,in,jn->n", amb, y, y)
    norm = 2 * norm / w
    y /= np.sqrt(norm)[None]

    x = lib.einsum("ij,jx,x->ix", amb, y, 1 / w)

    return w, x, y


def get_veff(mc, external=slice(None), internal=slice(None)):
    """Get a block of the effective potential.

    Args:
        mc: The CASSCF object.
        external: The block of the effective potential to sum over density matrix elements. If a
            tuple, specifies the blocks for each of the two external indices.
        internal: The block of the effective potential to return.

    Returns:
        The block of the effective potential.
    """
    mo_coeff_int = mc.mo_coeff[:, internal]
    mo_occ_int = mc.mo_occ[internal]
    if isinstance(external, tuple):
        mo_coeff_ext = (mc.mo_coeff[:, external[0]], mc.mo_coeff[:, external[1]])
    else:
        mo_coeff_ext = (mc.mo_coeff[:, external],) * 2

    # Evaluate the effective potential in the AO basis
    dm = lib.einsum("i,pi,qi->pq", mo_occ_int, mo_coeff_int, mo_coeff_int)
    veff = mc._scf.get_veff(mc._scf._eri, dm)

    # Rotate the effective potential to the NO basis
    veff = lib.einsum("pq,pi,qj->ij", veff, mo_coeff_ext[0], mo_coeff_ext[1])

    return veff


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

    a = lib.einsum("ps,qs,pr->pqrs", mo_occ_difs, h1e[act, act], np.eye(nact))
    a += lib.einsum("qr,pr,sq->pqrs", mo_occ_difs, h1e[act, act], np.eye(nact))

    a += lib.einsum("sqtu,purt->pqrs", h2e_aaaa, rdm2)
    a += lib.einsum("sutq,putr->pqrs", h2e_aaaa, rdm2)

    a += lib.einsum("utpr,stqu->pqrs", h2e_aaaa, rdm2)
    a += lib.einsum("urpt,stuq->pqrs", h2e_aaaa, rdm2)

    a -= lib.einsum("ptsu,tuqr->pqrs", h2e_aaaa, rdm2)

    a -= lib.einsum("tqur,sput->pqrs", h2e_aaaa, rdm2)

    w_aa = lib.einsum("twpu,wutr->pr", h2e_aaaa, rdm2)
    w_aa += lib.einsum("tupw,wurt->pr", h2e_aaaa, rdm2)
    a -= lib.einsum("pr,qs->pqrs", w_aa, np.eye(nact)) * 0.5
    a -= lib.einsum("pr,qs->pqrs", np.eye(nact), w_aa) * 0.5

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


def _build_act_cor_0(mc, h1e, rdm2):
    """Build the active-core block of the α=0 A+B and A-B matrices.

    Args:
        mc: The CASSCF object.
        h1e: The 1-electron α=0 Hamiltonian in a basis of natural orbitals.
        rdm2: The 2-particle CAS reduced density matrix in a basis of natural orbitals, for the
            active space.
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
    mo_occ = mc.mo_occ[occ] * 0.5
    mo_occ_sqrt = np.sqrt(mo_occ)
    mo_occ_sqrt[mo_occ < 0.5] *= -1
    mo_occ_difs = np.subtract.outer(mo_occ, mo_occ)

    # Get the Hamiltonian blocks
    h2e_cccc = ao2mo.kernel(mc._scf._eri, mo_coeff_cor, compact=False).reshape((ncor,) * 4)
    h2e_aaaa = ao2mo.kernel(mc._scf._eri, mo_coeff_act, compact=False).reshape((nact,) * 4)

    # Note: We only need A_{pqrq}, so re-order as A_{qpr} for efficiency

    a = lib.einsum("pq,qq,pr->qpr", mo_occ_difs[act, cor], h1e[cor, cor], np.eye(nact))
    a += lib.einsum("qr,pr,qq->qpr", mo_occ_difs[cor, act], h1e[act, act], np.eye(ncor))

    tmp = get_veff(mc, cor, cor)
    a += lib.einsum("pr,qq->qpr", np.diag(mo_occ[act]), tmp)

    a += lib.einsum("utpr,qq,tu->qpr", h2e_aaaa, np.diag(mo_occ)[cor, cor], np.diag(mo_occ)[act, act]) * 2
    a -= lib.einsum("utpr,qu,tq->qpr", h2e_aaaa, np.diag(mo_occ)[cor, act], np.diag(mo_occ)[act, cor])
    a += lib.einsum("urpt,qu,tq->qpr", h2e_aaaa, np.diag(mo_occ)[cor, act], np.diag(mo_occ)[act, cor]) * 2
    a -= lib.einsum("urpt,qq,tu->qpr", h2e_aaaa, np.diag(mo_occ)[cor, cor], np.diag(mo_occ)[act, act])

    w_aa = lib.einsum("twpu,wutr->pr", h2e_aaaa, rdm2)
    w_aa += lib.einsum("tupw,wurt->pr", h2e_aaaa, rdm2)
    a -= lib.einsum("pr,qq->qpr", w_aa, np.eye(ncor)) * 0.5

    w_cc = lib.einsum("twpu,wt,ur->pr", h2e_cccc, np.diag(mo_occ[cor]), np.diag(mo_occ[cor])) * 2
    w_cc -= lib.einsum("twpu,wr,ut->pr", h2e_cccc, np.diag(mo_occ[cor]), np.diag(mo_occ[cor]))
    w_cc += lib.einsum("tupw,wr,ut->pr", h2e_cccc, np.diag(mo_occ[cor]), np.diag(mo_occ[cor])) * 2
    w_cc -= lib.einsum("tupw,wt,ur->pr", h2e_cccc, np.diag(mo_occ[cor]), np.diag(mo_occ[cor]))
    a -= lib.einsum("pr,qq->qpr", np.eye(nact), w_cc) * 0.5

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
    mo_coeff_vir = mc.mo_coeff[:, vir]

    # Get the occupancy measures
    mo_occ_act = mc.mo_occ[act] * 0.5
    mo_occ_vir = mc.mo_occ[vir] * 0.5
    mo_occ_act_sqrt = np.sqrt(mo_occ_act)
    mo_occ_vir_sqrt = np.sqrt(mo_occ_vir)
    mo_occ_act_sqrt[mo_occ_act < 0.5] *= -1
    mo_occ_vir_sqrt[mo_occ_vir < 0.5] *= -1
    mo_occ_vir_act_difs = np.subtract.outer(mo_occ_vir, mo_occ_act)
    mo_occ_act_vir_difs = np.subtract.outer(mo_occ_act, mo_occ_vir)

    # Get the Hamiltonian blocks
    h2e_aaaa = ao2mo.kernel(mc._scf._eri, mo_coeff_act, compact=False).reshape((nact,) * 4)

    # Note: We only need A_{pqps}, so re-order as A_{pqs} for efficiency

    a = lib.einsum("ps,qs,pp->pqs", mo_occ_vir_act_difs, h1e[act, act], np.eye(nvir))
    a += lib.einsum("qp,pp,sq->pqs", mo_occ_act_vir_difs, h1e[vir, vir], np.eye(nact))

    w_aa = lib.einsum("twpu,wutr->pr", h2e_aaaa, rdm2)
    w_aa += lib.einsum("tupw,wurt->pr", h2e_aaaa, rdm2)
    a -= lib.einsum("pp,qs->pqs", np.eye(nvir), w_aa) * 0.5

    apb = []
    amb = []
    for p in range(nvir):
        # Get the index helpers
        mask = np.abs(mo_occ_vir[p] - mo_occ_act) > 1e-8
        if not np.any(mask):
            continue
        q = np.arange(nact)[mask]
        tril = np.tril_indices(q.size)
        triu = np.triu_indices(q.size)

        # Pack the A+B matrix
        apb_p = np.zeros((q.size, q.size))
        apb_p[triu] = a[p].reshape(q.size, q.size)[triu]
        apb_p[tril] = a[p].reshape(q.size, q.size).T[tril]
        m = mo_occ_vir_sqrt[p] + mo_occ_act_sqrt
        m[m == 0] = 1
        apb_p /= np.multiply.outer(m, m)
        apb.append(apb_p)

        # Pack the A-B matrix
        amb_p = np.zeros((q.size, q.size))
        amb_p[triu] = a[p].reshape(q.size, q.size)[triu]
        amb_p[tril] = a[p].reshape(q.size, q.size).T[tril]
        m = mo_occ_vir_sqrt[p] - mo_occ_act_sqrt
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
    mo_coeff_vir = mc.mo_coeff[:, vir]

    # Get the occupancy measures
    mo_occ_cor = mc.mo_occ[cor] * 0.5
    mo_occ_vir = mc.mo_occ[vir] * 0.5
    mo_occ_cor_sqrt = np.sqrt(mo_occ_cor)
    mo_occ_vir_sqrt = np.sqrt(mo_occ_vir)
    mo_occ_cor_sqrt[mo_occ_cor < 0.5] *= -1
    mo_occ_vir_sqrt[mo_occ_vir < 0.5] *= -1
    mo_occ_vir_cor_difs = np.subtract.outer(mo_occ_vir, mo_occ_cor)
    mo_occ_cor_vir_difs = np.subtract.outer(mo_occ_cor, mo_occ_vir)

    # Get the Hamiltonian blocks
    h2e_cccc = ao2mo.kernel(mc._scf._eri, mo_coeff_cor, compact=False).reshape((ncor,) * 4)

    # Note: We only need A_{pqpq}, so re-order as A_{pq} for efficiency

    a = lib.einsum("pq,qq->pq", mo_occ_vir_cor_difs, h1e[cor, cor])
    a += lib.einsum("qp,pp->pq", mo_occ_cor_vir_difs, h1e[vir, vir])

    w_cc = lib.einsum("twpu,wt,ur->pr", h2e_cccc, np.diag(mo_occ_cor), np.diag(mo_occ_cor)) * 2
    w_cc -= lib.einsum("twpu,wr,ut->pr", h2e_cccc, np.diag(mo_occ_cor), np.diag(mo_occ_cor))
    w_cc += lib.einsum("tupw,wr,ut->pr", h2e_cccc, np.diag(mo_occ_cor), np.diag(mo_occ_cor)) * 2
    w_cc -= lib.einsum("tupw,wt,ur->pr", h2e_cccc, np.diag(mo_occ_cor), np.diag(mo_occ_cor))
    a -= lib.einsum("pp,qq->pq", np.eye(nvir), w_cc) * 0.5

    # Pack the A+B and A-B matrices
    apb = amb = [np.array([[x]]) for x in a.ravel()]

    return apb, amb


def _calculate_energy(mc, h1e, rdm2, w_0, x_0, y_0):
    """Build a block of the α=1 A+B and A-B matrices.

    Args:
        mc: The CASSCF object.
        h1e: The 1-electron α=1 Hamiltonian in a basis of natural orbitals.
        rdm2: The 2-particle CAS reduced density matrix in a basis of natural orbitals, for the
            active space.
        w_0: The eigenvalues of the α=0 A+B and A-B matrices.
        x_0: The X eigenvectors of the α=0 A+B and A-B matrices.
        y_0: The Y eigenvectors of the α=0 A+B and A-B matrices.
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
    pq = []
    for p in range(nact):
        for q in range(p):
            if pairs[p+act.start, q+act.start]:
                pq.append((p+act.start, q+act.start))
    for q in range(ncor):
        for p in range(nact):
            if pairs[p+act.start, q+cor.start]:
                pq.append((p+act.start, q+cor.start))
    for p in range(nvir):
        for q in range(nact):
            if pairs[p+vir.start, q+act.start]:
                pq.append((p+vir.start, q+act.start))
    for p in range(nvir):
        for q in range(ncor):
            if pairs[p+vir.start, q+cor.start]:
                pq.append((p+vir.start, q+cor.start))
    p, q = zip(*pq)
    p = np.array(p)
    q = np.array(q)
    r, s = p, q
    tril = np.tril_indices(p.size)
    triu = np.triu_indices(p.size)

    # Get the Hamiltonian blocks
    h2e = ao2mo.kernel(mc._scf._eri, mc.mo_coeff, compact=False).reshape((norb,) * 4)

    # Contractions
    a = lib.einsum("ps,qs,pr->pqrs", mo_occ_difs, h1e, np.eye(norb))

    tmp = lib.einsum("pqrs,p,r->pqrs", h2e[cor, :, not_cor, :], mo_occ[cor], mo_occ[not_cor]) * 2
    tmp -= lib.einsum("prqs,p,r->pqrs", h2e[cor, not_cor, :, :], mo_occ[cor], mo_occ[not_cor])
    a[cor, :, not_cor, :] += tmp
    a[:, not_cor, :, cor] += tmp.transpose(3, 2, 1, 0)

    tmp = lib.einsum("qprs,q,r->pqrs", h2e[cor, :, not_cor, :], mo_occ[cor], mo_occ[not_cor]) * 2
    tmp -= lib.einsum("qsrp,q,r->pqrs", h2e[cor, :, not_cor, :], mo_occ[cor], mo_occ[not_cor])
    a[:, cor, not_cor, :] -= tmp
    a[not_cor, :, :, cor] -= tmp.transpose(2, 3, 0, 1)

    a[cor, :, cor, :] += lib.einsum("pqrs,p,r->pqrs", h2e[cor, :, cor, :], mo_occ[cor], mo_occ[cor]) * 2
    a[cor, :, cor, :] -= lib.einsum("prqs,p,r->pqrs", h2e[cor, cor, :, :], mo_occ[cor], mo_occ[cor])

    a[:, cor, cor, :] -= lib.einsum("qprs,q,r->pqrs", h2e[cor, :, cor, :], mo_occ[cor], mo_occ[cor]) * 2
    a[:, cor, cor, :] += lib.einsum("rpqs,q,r->pqrs", h2e[cor, :, cor, :], mo_occ[cor], mo_occ[cor])

    tmp = lib.einsum("t,pqtt->pq", mo_occ[not_vir], h2e[:, :, not_vir, not_vir]) * 2  #FIXME
    tmp -= lib.einsum("t,pttq->pq", mo_occ[not_vir], h2e[:, not_vir, not_vir, :])
    a[cor, :, cor, :] += lib.einsum("pr,qs->pqrs", np.diag(mo_occ[cor]), tmp)

    tmp = lib.einsum("t,pqtt->pq", mo_occ[cor], h2e[:, :, cor, cor]) * 2  #FIXME
    tmp -= lib.einsum("t,pttq->pq", mo_occ[cor], h2e[:, cor, cor, :])
    a[not_cor, :, not_cor, :] += lib.einsum("pr,qs->pqrs", np.diag(mo_occ[not_cor]), tmp)

    a[act, :, act, :] += lib.einsum("tusq,purt->pqrs", h2e[act, act, :, :], rdm2)
    a[act, :, act, :] += lib.einsum("ustq,putr->pqrs", h2e[act, :, act, :], rdm2)

    a[:, act, act, :] -= lib.einsum("tpus,tuqr->pqrs", h2e[act, :, act, :], rdm2)

    tmp = np.zeros((norb, norb))
    tmp[:, act] += lib.einsum("twup,wutr->pr", h2e[act, act, act, :], rdm2)
    tmp[:, act] += lib.einsum("tuwp,wurt->pr", h2e[act, act, act, :], rdm2)
    for slc, sign in ((not_vir, 1), (act, -1)):
        # Avoids double counting of the active part
        tmp[:, slc] += lib.einsum("twup,wt,ur->pr", h2e[slc, slc, slc, :], np.diag(mo_occ[slc]), np.diag(mo_occ[slc])) * 2 * sign
        tmp[:, slc] -= lib.einsum("twup,wr,ut->pr", h2e[slc, slc, slc, :], np.diag(mo_occ[slc]), np.diag(mo_occ[slc])) * sign
        tmp[:, slc] += lib.einsum("tuwp,wr,ut->pr", h2e[slc, slc, slc, :], np.diag(mo_occ[slc]), np.diag(mo_occ[slc])) * 2 * sign
        tmp[:, slc] -= lib.einsum("tuwp,wt,ur->pr", h2e[slc, slc, slc, :], np.diag(mo_occ[slc]), np.diag(mo_occ[slc])) * sign
    a -= lib.einsum("pr,qs->pqrs", tmp, np.eye(norb)) * 0.5

    # Symmetrise
    a = a + a.transpose(1, 0, 3, 2)

    # Find A+B
    apb = np.zeros((p.size, p.size))
    apb[triu] += a[p, q][:, r, s][triu]
    apb[triu] += a[q, p][:, r, s][triu]
    apb[tril] = apb.T[tril]
    m = np.add.outer(mo_occ_sqrt, mo_occ_sqrt)[p, q]
    m[m == 0] = 1
    apb /= np.multiply.outer(m, m)

    # Find A-B
    amb = np.zeros((p.size, p.size))
    amb[triu] += a[p, q][:, r, s][triu]
    amb[triu] -= a[q, p][:, r, s][triu]
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
    block = scipy.linalg.block_diag(*[np.ones((n, n)) for n in [ncor, nact, nvir]])
    eri_s4 = h2e[p, q][:, r, s] - lib.einsum("pqrs,pq,qr,rs->pqrs", h2e, block, block, block)[p, q][:, r, s]
    e_corr = np.sum(xm * eri_s4)

    return e_corr


def run(mf, mc):
    nact = mc.ncas
    nelec = mc.nelecas

    rdm2 = mc.fcisolver.make_rdm2(mc.ci, mc.ncas, mc.nelecas)

    norb = mf.mo_occ.size
    nocc = mc.ncore
    nact = mc.ncas
    nvir = norb - nocc - nact

    occ = slice(0, nocc)
    act = slice(nocc, nocc + nact)
    vir = slice(nocc + nact, norb)

    # Get the 1-electron Hamiltonian
    hcore = lib.einsum("pq,pi,qj->ij", mf.get_hcore(), mc.mo_coeff, mc.mo_coeff)
    fock = hcore + get_veff(mc)
    rdm2 = rdm2.transpose(1, 3, 0, 2) * 0.5
    block = scipy.linalg.block_diag(*[np.ones((n, n)) for n in [nocc, nact, nvir]])

    # Get the α=0 1e Hamiltonian
    h1e = fock * block
    for slc in (occ, act, vir):
        h1e[slc, slc] -= get_veff(mc, slc, slc)

    # Active-active
    apb_act_act_0, amb_act_act_0 = _build_act_act_0(mc, h1e, rdm2)
    w_act_act_0, x_act_act_0, y_act_act_0 = eigxy(apb_act_act_0, amb_act_act_0)

    # Active-occupied
    apb_act_occ_0, amb_act_occ_0 = _build_act_cor_0(mc, h1e, rdm2)
    if apb_act_occ_0:
        w_act_occ_0, x_act_occ_0, y_act_occ_0 = zip(*[eigxy(apb, amb) for apb, amb in zip(apb_act_occ_0, amb_act_occ_0)])
    else:
        w_act_occ_0, x_act_occ_0, y_act_occ_0 = [], [], []

    # Virtual-active
    apb_vir_act_0, amb_vir_act_0 = _build_vir_act_0(mc, h1e, rdm2)
    if apb_vir_act_0:
        w_vir_act_0, x_vir_act_0, y_vir_act_0 = zip(*[eigxy(apb, amb) for apb, amb in zip(apb_vir_act_0, amb_vir_act_0)])
    else:
        w_vir_act_0, x_vir_act_0, y_vir_act_0 = [], [], []

    # Virtual-occupied
    apb_vir_occ_0, amb_vir_occ_0 = _build_vir_cor_0(mc, h1e, rdm2)
    w_vir_occ_0 = [np.array(a.ravel()) for a in apb_vir_occ_0]
    x_vir_occ_0 = [np.ones_like(a) / np.sqrt(2.0) for a in apb_vir_occ_0]
    y_vir_occ_0 = [np.ones_like(a) / np.sqrt(2.0) for a in apb_vir_occ_0]

    # Collect
    w_0 = np.concatenate([w_act_act_0, *w_act_occ_0, *w_vir_act_0, *w_vir_occ_0])
    x_0 = scipy.linalg.block_diag(x_act_act_0, *x_act_occ_0, *x_vir_act_0, *x_vir_occ_0)
    y_0 = scipy.linalg.block_diag(y_act_act_0, *y_act_occ_0, *y_vir_act_0, *y_vir_occ_0)

    # Get the α=1 1e Hamiltonian
    h1e = hcore * (1 - block)
    for eslc in (occ, act, vir):
        for islc in (occ, act, vir):
            if eslc != islc:
                h1e[eslc, eslc] -= get_veff(mc, eslc, islc)

    # Calculate the energy
    e_corr = _calculate_energy(mc, h1e, rdm2, w_0, x_0, y_0)

    return e_corr


if __name__ == "__main__":
    import time
    from memory_profiler import memory_usage
    from pyscf import gto, scf, mcscf
    from pyscf.cas_ac0.accas import get_cas_ac0_energy

    for mol, cas in [
        (gto.M(atom="O 0 0 0; H 0 0 1; H 0 1 0", basis="6-31g", verbose=0), (6, 4)),
        (gto.M(atom="Li 0 0 0; H 0 0 2", basis="sto3g", verbose=0), (4, 2)),
        (gto.M(atom="Li 0 0 0; H 0 0 2", basis="cc-pvdz", verbose=0), (4, 2)),
        (gto.M(atom="Li 0 0 0; Li 0 0 2", basis="6-31g", verbose=0), (6, 6)),
        (gto.M(atom="O 0 0 0; H 0 0 1; H 0 1 0", basis="sto3g", verbose=0), (2, 0)),
        (gto.M(atom="O 0 0 0; H 0 0 1; H 0 1 0", basis="cc-pvdz", verbose=0), (2, 0)),
        (gto.M(atom="O 0 0 0; O 0 0 1", basis="aug-cc-pvdz", verbose=0), (2, 0)),
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
        print("old: %6.2f ms  %6.1f mb   new: %6.2f ms  %6.1f mb" % (t2 * 1000, m2, t1 * 1000, m1))
