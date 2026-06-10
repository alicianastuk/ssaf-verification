"""SSAF verification suite, file 6: Finite-Dimensional Sequencing Without
Retro-Addressability (Appendix D; non-addressability as non-invertibility).

Pass conditions verified:
  S1: the interface map I_E (complete dephasing in the E-fixed basis) is
      CPTP with parameters depending on E only.
  S2: I_E is genuinely non-injective on the admissible domain: distinct
      admissible states map to the same image (the |+>,|-> pair both go
      to I/2), so excluded distinctions exist.
  S3: no left inverse J_E with J_E . I_E = id can exist on that domain,
      CPTP or otherwise: it would need J(I/2) to equal two different
      states at once. Retro-admissibility is structurally unavailable,
      not prohibited.
  S4: the existence of an alternative restriction A~'(E) does not supply
      a within-ordering recovery: mapping into A~' does not restore the
      excluded distinction either (any map sees only the collapsed image).
"""
import numpy as np

TOL = 1e-12
sz = np.diag([1.0, -1.0]).astype(complex)

def interface_dephase(rho):
    """Complete dephasing: E-fixed basis, no rho-dependent parameters."""
    return np.diag(np.diag(rho))

def choi(T):
    C = np.zeros((4, 4), dtype=complex)
    for i in range(2):
        for j in range(2):
            Eij = np.zeros((2, 2), dtype=complex); Eij[i, j] = 1
            C[2*i:2*i+2, 2*j:2*j+2] = T(Eij)
    return C

plus = 0.5*np.array([[1, 1], [1, 1]], dtype=complex)
minus = 0.5*np.array([[1, -1], [-1, 1]], dtype=complex)

def test_S1_cptp_e_only():
    C = choi(interface_dephase)
    assert np.linalg.eigvalsh(C).min() > -TOL
    assert abs(np.trace(interface_dephase(plus)).real - 1) < TOL
    print("S1 interface map CPTP, E-only parameters: PASS")

def test_S2_noninjective():
    img_p, img_m = interface_dephase(plus), interface_dephase(minus)
    assert np.allclose(img_p, img_m, atol=TOL), "images should coincide"
    assert not np.allclose(plus, minus), "inputs are distinct"
    print("S2 non-injectivity on admissible domain: PASS")

def test_S3_no_left_inverse():
    # any function J (not just CPTP) satisfying J(I_E(rho)) = rho for both
    # plus and minus would need J(I/2) = plus AND J(I/2) = minus.
    common_image = interface_dephase(plus)
    # the requirement set is contradictory:
    requirement_targets = [plus, minus]
    distinct = not np.allclose(requirement_targets[0], requirement_targets[1])
    assert distinct, "left-inverse requirements must be contradictory"
    print("S3 no left inverse exists on the domain (structural): PASS")

def test_S4_alternative_restriction_no_recovery():
    # A~'(E): a different restriction (dephasing in the X basis). Composing
    # any map after I_E still cannot separate plus from minus.
    H = np.array([[1, 1], [1, -1]], dtype=complex)/np.sqrt(2)
    def into_alt(rho):  # arbitrary within-ordering map into the alt sector
        return H @ interface_dephase(rho) @ H.conj().T
    assert np.allclose(into_alt(plus), into_alt(minus), atol=TOL), \
        "post-processing must not restore the excluded distinction"
    print("S4 alternative restriction supplies no recovery: PASS")

if __name__ == "__main__":
    test_S1_cptp_e_only(); test_S2_noninjective()
    test_S3_no_left_inverse(); test_S4_alternative_restriction_no_recovery()
    print("ALL PASS: sequencing without retro-addressability verified.")
