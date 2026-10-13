"""Radical-pair spin model: singlet yield for a singlet-born pair with Haberkorn recombination.

Two independent implementations of the yield:
  * ``yield_eig``       eigenbasis Lorentzian formula (fast, vectorised over directions)
  * ``yield_liouville`` direct Liouville-space Haberkorn solution (slow, also supports relaxation)

Units: hyperfine/exchange/dipolar couplings are given in mT; angular frequency = GAMMA_E * (mT value).
Molecular frame: z = flavin symmetry axis.  Hyperfine principal axes coincide with x, y, z (ASSUMED).
"""
from __future__ import annotations

import numpy as np

GAMMA_E = 1.76085963023e8  # rad s^-1 mT^-1 (free electron, CODATA 2018)

# --- verified tensor set (Hiscock et al., PNAS 113:4634 (2016), main text, section "Origin of the Spike in PhiS")
H16 = {
    "N5": (-0.087, -0.100, 1.757),
    "N10": (-0.014, -0.024, 0.605),
}
# --- recalled in the study brief; UNVERIFIED and different from the primary source (sensitivity set only)
RECALLED = {
    "N5": (-0.0989, -0.0987, 1.7569),
    "N10": (-0.0035, -0.0021, 0.4265),
}


def spin_ops(s: float):
    d = int(round(2 * s + 1))
    m = np.arange(s, -s - 1, -1.0)
    sz = np.diag(m).astype(complex)
    sp = np.zeros((d, d), complex)
    for i in range(1, d):
        mm = m[i]
        sp[i - 1, i] = np.sqrt(s * (s + 1) - mm * (mm + 1))
    sx = (sp + sp.conj().T) / 2
    sy = (sp - sp.conj().T) / (2j)
    return sx, sy, sz


class RPModel:
    """nuc_a / nuc_b: list of (spin, (Axx,Ayy,Azz)) or (spin, 3x3 array) in mT."""

    def __init__(self, nuc_a, nuc_b=(), J=0.0, D=0.0, name=""):
        self.name = name
        spins = [0.5, 0.5] + [s for s, _ in nuc_a] + [s for s, _ in nuc_b]
        dims = [int(round(2 * s + 1)) for s in spins]
        self.dim = int(np.prod(dims))
        self.M = self.dim // 4

        def op(i, o):
            out = np.array([[1.0 + 0j]])
            for j, d in enumerate(dims):
                out = np.kron(out, o if j == i else np.eye(d))
            return out

        ops = [spin_ops(s) for s in spins]
        self.SA = [op(0, ops[0][a]) for a in range(3)]
        self.SB = [op(1, ops[1][a]) for a in range(3)]
        H0 = np.zeros((self.dim, self.dim), complex)
        idx = 2
        for side, nucs in ((self.SA, nuc_a), (self.SB, nuc_b)):
            for s, A in nucs:
                A = np.diag(A).astype(float) if np.ndim(A) == 1 else np.asarray(A, float)
                I = [op(idx, ops[idx][a]) for a in range(3)]
                for a in range(3):
                    for b in range(3):
                        if A[a, b] != 0.0:
                            H0 += GAMMA_E * A[a, b] * side[a] @ I[b]
                idx += 1
        sasb = sum(self.SA[a] @ self.SB[a] for a in range(3))
        self.PS = 0.25 * np.eye(self.dim) - sasb
        if J != 0.0:
            H0 += -GAMMA_E * J * (2 * sasb + 0.5 * np.eye(self.dim))
        if D != 0.0:  # axial point-dipole form about molecular z (ASSUMED axis)
            H0 += -GAMMA_E * D * (3 * self.SA[2] @ self.SB[2] - sasb)
        self.H0 = H0
        self.Z = [self.SA[a] + self.SB[a] for a in range(3)]
        self._liou = None

    # ------------------------------------------------------------------ eigenbasis implementation
    def hamiltonian(self, n, B0_mT):
        n = np.atleast_2d(n)
        zee = GAMMA_E * B0_mT * sum(n[:, a, None, None] * self.Z[a][None] for a in range(3))
        return self.H0[None] + zee

    def yield_eig(self, n, B0_mT, k, chunk=None):
        n = np.atleast_2d(np.asarray(n, float))
        if chunk is None:  # bound memory for large Hilbert spaces
            chunk = int(max(1, min(256, 4e6 // self.dim**2)))
        out = np.empty(len(n))
        for s in range(0, len(n), chunk):
            H = self.hamiltonian(n[s : s + chunk], B0_mT)
            E, V = np.linalg.eigh(H)
            P = np.conj(np.swapaxes(V, 1, 2)) @ self.PS[None] @ V
            w = np.abs(P) ** 2
            dE = E[:, :, None] - E[:, None, :]
            out[s : s + chunk] = (w * (k * k / (k * k + dE**2))).sum((1, 2)) / self.M
        return out

    # ------------------------------------------------------------------ Liouville implementation
    def _prep_liouville(self):
        d = self.dim
        I = np.eye(d)

        def comm(A):  # row-major vec: vec(A rho - rho A) = (A kron I - I kron A^T) vec(rho)
            return -1j * (np.kron(A, I) - np.kron(I, A.T))

        self._liou = {
            "L0": comm(self.H0),
            "LZ": [comm(self.Z[a]) for a in range(3)],
            "rho0": (self.PS / self.M).reshape(-1),
            "p": self.PS.T.reshape(-1),
            "relax": sum(np.kron(S, S.T) for S in (self.SA + self.SB)) - 1.5 * np.eye(d * d),
        }

    def yield_liouville(self, n, B0_mT, k, relax=0.0):
        """Phi = k Tr[P_S (k - L)^-1 rho0], L = -i[H,.] + relax*D (depolarising, both electrons)."""
        if self._liou is None:
            self._prep_liouville()
        c = self._liou
        n = np.atleast_2d(np.asarray(n, float))
        out = np.empty(len(n))
        for i, nn in enumerate(n):
            L = c["L0"] + GAMMA_E * B0_mT * sum(nn[a] * c["LZ"][a] for a in range(3))
            if relax:
                L = L + relax * c["relax"]
            x = np.linalg.solve(k * np.eye(self.dim**2) - L, c["rho0"])
            out[i] = (k * (c["p"] @ x)).real
        return out


def yield_from_hamiltonian(H, PS, M, k):
    """Single dense Hamiltonian version (used by the H = 0 check)."""
    E, V = np.linalg.eigh(H)
    P = V.conj().T @ PS @ V
    dE = E[:, None] - E[None, :]
    return float((np.abs(P) ** 2 * (k * k / (k * k + dE**2))).sum() / M)


def primary_model(tensors=H16):
    return RPModel([(1, tensors["N5"]), (1, tensors["N10"])], name="primary")


def isotropic_version(tensors=H16):
    iso = {kk: (sum(v) / 3,) * 3 for kk, v in tensors.items()}
    return RPModel([(1, iso["N5"]), (1, iso["N10"])], name="isotropic")
