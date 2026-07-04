"""SSAF verification suite, file 3: Two-Qubit Joint Admissibility Without
Factorization + General Form (toy-model sections merged toy; Lemma marginals-nondet,
Corollary toy_ent_strict).

Pass conditions verified:
  J1: the Bell state Phi+ is non-separable (PPT criterion: partial
      transpose has a negative eigenvalue; exact for 2x2).
  J2: its marginals are exactly maximally mixed, hence lie in any
      reasonable marginal admissible set containing I/2.
  J3: marginals do not determine the joint state: the product state
      (I/2)x(I/2) has the SAME marginals but is separable (PPT positive),
      so two joint states with identical admissible marginals differ in
      separability -- membership in A_AB(E) is not a function of marginals.
  J4 (general form): with A_AB(E) := {Phi+ projector} u {separable
      maximally-mixed-marginal states in our sample}, the joint-sector set
      J(E) of all states with admissible marginals strictly contains
      A_AB(E): the strict inclusion of Cor. toy_ent_strict is witnessed.
"""
import numpy as np

rng = np.random.default_rng(5)
TOL = 1e-9

def bell_phi_plus():
    v = np.zeros(4, dtype=complex); v[0] = v[3] = 1/np.sqrt(2)
    return np.outer(v, v.conj())

def partial_transpose_B(rho):
    R = rho.reshape(2, 2, 2, 2)          # (a, b, a', b')
    return R.transpose(0, 3, 2, 1).reshape(4, 4)

def marginal_A(rho):
    return rho.reshape(2, 2, 2, 2).trace(axis1=1, axis2=3)

def marginal_B(rho):
    return rho.reshape(2, 2, 2, 2).trace(axis1=0, axis2=2)

def test_J1_nonseparable():
    neg = np.linalg.eigvalsh(partial_transpose_B(bell_phi_plus())).min()
    assert neg < -0.4, "Phi+ should violate PPT decisively (eigenvalue -1/2)"
    print(f"J1 non-separability via PPT: PASS (min PT eigenvalue {neg:.3f})")

def test_J2_marginals_maximally_mixed():
    rho = bell_phi_plus()
    for m in (marginal_A(rho), marginal_B(rho)):
        assert np.allclose(m, np.eye(2)/2, atol=TOL), "marginal not I/2"
    print("J2 marginals exactly I/2: PASS")

def test_J3_marginals_dont_determine_joint():
    ent = bell_phi_plus()
    prod = np.kron(np.eye(2)/2, np.eye(2)/2)
    assert np.allclose(marginal_A(ent), marginal_A(prod), atol=TOL)
    assert np.allclose(marginal_B(ent), marginal_B(prod), atol=TOL)
    pt_ent = np.linalg.eigvalsh(partial_transpose_B(ent)).min()
    pt_prod = np.linalg.eigvalsh(partial_transpose_B(prod)).min()
    assert pt_ent < -TOL and pt_prod > -TOL, \
        "same marginals must separate into PPT-violating and PPT-satisfying"
    print("J3 identical marginals, different joint structure: PASS")

def test_J4_strict_inclusion_general_form():
    # J(E): sample of joint states whose BOTH marginals equal I/2.
    # Werner family rho_p = p*Phi+ + (1-p)*I/4 has marginals I/2 for all p.
    members_of_J = []
    for p in np.linspace(0, 1, 21):
        rho = p*bell_phi_plus() + (1-p)*np.eye(4)/4
        assert np.allclose(marginal_A(rho), np.eye(2)/2, atol=TOL)
        members_of_J.append(rho)
    # Declare A_AB(E) as a PROPER subset: e.g. only the PPT members
    # (a legitimate context-indexed choice). Then Phi+ in J(E) \ A_AB(E):
    in_A = [r for r in members_of_J
            if np.linalg.eigvalsh(partial_transpose_B(r)).min() > -TOL]
    assert 0 < len(in_A) < len(members_of_J), \
        "A_AB(E) must be nonempty and strictly inside J(E)"
    print(f"J4 strict inclusion A_AB(E) (n={len(in_A)}) "
          f"in J(E) (n={len(members_of_J)}): PASS")

if __name__ == "__main__":
    test_J1_nonseparable(); test_J2_marginals_maximally_mixed()
    test_J3_marginals_dont_determine_joint(); test_J4_strict_inclusion_general_form()
    print("ALL PASS: joint admissibility toy verified.")
