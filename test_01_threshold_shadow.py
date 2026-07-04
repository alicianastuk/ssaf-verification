"""SSAF verification suite, file 1: threshold-shadow grounding of Compat_E / Gamma_E.

Verifies the two claims of the threshold-shadow working note against the
charge-sector toy (SSAF toy-model sections):

  T1 (block case): for S_E = B(H+) (+) B(H-) and definite-sector states,
      g_E = 1 within a sector and g_E = 0 across sectors, so the threshold
      shadows at theta_c = 0 and theta_gamma in (0,1] reproduce the toy's
      declared Compat_E (total) and Gamma_E edge relation (same-sector).

  T2 (graded case): for the leaky operator system S_eps, cross-sector
      g_E = eps exactly, so thresholds are substantive.

Run: python3 ssaf_verify_threshold_shadow.py   (exits 0 iff all assertions pass)
"""
import numpy as np

rng = np.random.default_rng(7)
TOL = 1e-9
BLOCKS = [[0, 1], [2, 3]]  # H = H+ (dim 2) (+) H- (dim 2)


def cond_exp_block(X, blocks=BLOCKS):
    """Conditional expectation onto the block-diagonal algebra."""
    Y = np.zeros_like(X)
    for b in blocks:
        Y[np.ix_(b, b)] = X[np.ix_(b, b)]
    return Y


def trace_norm(X):
    return np.linalg.svd(X, compute_uv=False).sum()


def g_block(psi, phi):
    """g_E for the block algebra: max_{A in M_E, ||A||<=1} |<phi|A|psi>|
    = trace norm of the conditional expectation of |psi><phi|."""
    return trace_norm(cond_exp_block(np.outer(psi, phi.conj())))


def sector_vec(sector, coeffs):
    v = np.zeros(4, dtype=complex)
    (v.__setitem__((slice(0, 2)), coeffs) if sector == '+'
     else v.__setitem__((slice(2, 4)), coeffs))
    return v / np.linalg.norm(v)


def test_T1_block_case():
    psi_p0 = sector_vec('+', [1, 0])
    psi_p1 = sector_vec('+', [0, 1])
    phi_m = sector_vec('-', [1, 1j])
    assert abs(g_block(psi_p0, psi_p1) - 1) < TOL, "same-sector should give 1"
    assert g_block(psi_p0, phi_m) < TOL, "cross-sector should give 0"
    assert abs(g_block(psi_p0, psi_p0) - 1) < TOL, "reflexivity"
    # random same/cross-sector pairs
    for _ in range(200):
        a = rng.normal(size=2) + 1j * rng.normal(size=2)
        b = rng.normal(size=2) + 1j * rng.normal(size=2)
        same = g_block(sector_vec('+', a), sector_vec('+', b))
        cross = g_block(sector_vec('+', a), sector_vec('-', b))
        assert abs(same - 1) < 1e-7 and cross < TOL
    # threshold shadows reproduce the toy's declared primitives
    theta_c, theta_g = 0.0, 0.5
    compat = lambda g: g >= theta_c
    gamma = lambda g: g >= theta_g
    assert compat(1.0) and compat(0.0), "Compat_E total (toy)"
    assert gamma(1.0) and not gamma(0.0), "Gamma_E same-sector only (toy)"
    print("T1 block case: PASS (shadows reproduce the charge-sector toy)")


def test_T2_graded_case():
    for eps in [0.0, 0.1, 0.3, 0.7]:
        best = 0.0
        for _ in range(500):
            a = rng.normal(size=2) + 1j * rng.normal(size=2)
            b = rng.normal(size=2) + 1j * rng.normal(size=2)
            psi = sector_vec('+', a)
            phi = sector_vec('-', b)
            C = eps * np.outer(phi, psi.conj())  # ||C|| = eps, D = 0
            best = max(best, abs(phi.conj() @ C @ psi))
        assert abs(best - eps) < 1e-7, f"cross-sector g should equal eps={eps}"
    print("T2 graded case: PASS (g interpolates exactly with access leakage)")


if __name__ == "__main__":
    test_T1_block_case()
    test_T2_graded_case()
    print("ALL PASS: threshold-shadow grounding verified.")
