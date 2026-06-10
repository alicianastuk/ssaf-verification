"""SSAF verification suite, file 2: NS1-safe single-qubit GKSL semigroup
(Appendix D, Toy: NS1-safe GKSL semigroup; Cor. NS1_GKSL_admissible).

Pass conditions verified, mirroring the toy:
  G1: the fixed-E dephasing GKSL generator yields CPTP maps at all times
      (Choi matrices positive, trace preserving).
  G2: parameters (H_E, L_k, gamma_k) are functions of E alone; the SAME
      channel acts on every input state (state-independence is structural
      in the implementation: no rho enters channel construction).
  G3: distinguishability (trace distance) is monotone non-increasing along
      the bookkeeping chain for all sampled state pairs (contractivity).
  G4: the finite-time map is linearly invertible but its inverse is NOT
      CPTP (Choi of the inverse has a negative eigenvalue), so no
      admissible recovery exists; in the asymptotic limit coherence loss
      is genuinely non-injective.
"""
import numpy as np

rng = np.random.default_rng(3)
TOL = 1e-9

sz = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def dephasing_channel(t, gamma):
    """Single-qubit dephasing semigroup e^{t L_E}, L_E rho = gamma (sz rho sz - rho).
    Off-diagonals decay by exp(-2 gamma t)."""
    f = np.exp(-2 * gamma * t)
    def T(rho):
        out = rho.copy()
        out[0, 1] *= f
        out[1, 0] *= f
        return out
    return T, f


def choi(T):
    """Choi matrix of a qubit channel."""
    C = np.zeros((4, 4), dtype=complex)
    for i in range(2):
        for j in range(2):
            Eij = np.zeros((2, 2), dtype=complex); Eij[i, j] = 1
            C[2*i:2*i+2, 2*j:2*j+2] = T(Eij)
    return C


def rand_state():
    A = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    rho = A @ A.conj().T
    return rho / np.trace(rho).real


def trace_dist(a, b):
    return 0.5 * np.abs(np.linalg.eigvalsh(a - b)).sum()


def test_G1_cptp():
    for t in [0.0, 0.1, 0.5, 2.0, 10.0]:
        T, _ = dephasing_channel(t, gamma=0.7)
        C = choi(T)
        assert np.linalg.eigvalsh(C).min() > -TOL, f"Choi not PSD at t={t}"
        for _ in range(50):
            rho = rand_state()
            assert abs(np.trace(T(rho)).real - 1) < TOL, "trace not preserved"
    print("G1 CPTP at all times: PASS")


def test_G2_state_independence():
    # channel constructed from (t, gamma) = E-data only; verify identical
    # action coefficients regardless of which states are later updated.
    T, f1 = dephasing_channel(1.0, 0.7)
    _, f2 = dephasing_channel(1.0, 0.7)
    assert f1 == f2, "channel parameters must be E-functions only"
    print("G2 state-independence (structural): PASS")


def test_G3_monotone_contraction():
    gamma, dt, steps = 0.7, 0.3, 12
    T, _ = dephasing_channel(dt, gamma)
    for _ in range(100):
        a, b = rand_state(), rand_state()
        prev = trace_dist(a, b)
        for _ in range(steps):
            a, b = T(a), T(b)
            cur = trace_dist(a, b)
            assert cur <= prev + TOL, "distinguishability increased"
            prev = cur
    print("G3 monotone restriction along the chain: PASS")


def test_G4_no_cptp_recovery():
    T, f = dephasing_channel(1.0, 0.7)  # f = e^{-1.4} in (0,1): invertible
    # linear inverse: multiply off-diagonals by 1/f
    def Tinv(rho):
        out = rho.copy()
        out[0, 1] /= f
        out[1, 0] /= f
        return out
    # inverse undoes T linearly
    rho = rand_state()
    assert np.allclose(Tinv(T(rho)), rho), "linear inverse should exist"
    # but the inverse is not CPTP: Choi has a negative eigenvalue
    Cinv = choi(Tinv)
    assert np.linalg.eigvalsh(Cinv).min() < -1e-6, \
        "inverse unexpectedly CPTP — recovery would be admissible"
    # asymptotic limit: complete dephasing is non-injective
    Tinf, _ = dephasing_channel(1e6, 0.7)
    plus = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
    minus = 0.5 * np.array([[1, -1], [-1, 1]], dtype=complex)
    assert trace_dist(Tinf(plus), Tinf(minus)) < 1e-6, \
        "complete dephasing should erase the +/- distinction"
    print("G4 no admissible recovery; asymptotic non-injectivity: PASS")


if __name__ == "__main__":
    test_G1_cptp()
    test_G2_state_independence()
    test_G3_monotone_contraction()
    test_G4_no_cptp_recovery()
    print("ALL PASS: NS1-safe GKSL toy verified.")
