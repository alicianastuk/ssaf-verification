"""SSAF verification suite, file 4: Correlation Without Influence Under
Local CPTP Updating (toy-model sections; Lemma Correlation Without Influence,
Prop. no_steerable_marginals).

Pass conditions verified:
  C1: for the Bell state, applying ANY sampled CPTP map on B alone leaves
      the marginal on A exactly invariant (no influence, no steering).
  C2: correlations persist: after local dephasing on B, the joint state
      remains correlated (Z-Z correlator nonzero) while marginals stay I/2.
  C3: the same marginal invariance holds for every sampled joint state,
      not just Bell states (linearity + partial-trace structure).
"""
import numpy as np

rng = np.random.default_rng(9)
TOL = 1e-9

def bell_phi_plus():
    v = np.zeros(4, dtype=complex); v[0] = v[3] = 1/np.sqrt(2)
    return np.outer(v, v.conj())

def rand_kraus_qubit(n_ops=3):
    """Random CPTP map on one qubit via Stinespring: isometry columns."""
    A = rng.normal(size=(2*n_ops, 2)) + 1j*rng.normal(size=(2*n_ops, 2))
    Q, _ = np.linalg.qr(A)            # 2n x 2 isometry
    return [Q[2*k:2*k+2, :] for k in range(n_ops)]

def apply_on_B(rho, kraus):
    out = np.zeros_like(rho)
    for K in kraus:
        KB = np.kron(np.eye(2), K)
        out += KB @ rho @ KB.conj().T
    return out

def marginal_A(rho):
    return rho.reshape(2, 2, 2, 2).trace(axis1=1, axis2=3)

def rand_joint_state():
    A = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
    rho = A @ A.conj().T
    return rho / np.trace(rho).real

def test_C1_marginal_invariance_bell():
    rho = bell_phi_plus()
    mA0 = marginal_A(rho)
    for _ in range(200):
        K = rand_kraus_qubit()
        assert np.allclose(marginal_A(apply_on_B(rho, K)), mA0, atol=1e-10), \
            "local update on B steered the A marginal"
    print("C1 A-marginal invariant under 200 random local B-updates: PASS")

def test_C2_correlation_persists():
    sz = np.diag([1.0, -1.0]).astype(complex)
    ZZ = np.kron(sz, sz)
    deph = [np.sqrt(0.5)*np.eye(2), np.sqrt(0.5)*sz]   # full dephasing on B
    rho = apply_on_B(bell_phi_plus(), deph)
    corr = np.trace(ZZ @ rho).real
    assert abs(corr - 1.0) < TOL, "Z-Z correlation should survive B-dephasing"
    assert np.allclose(marginal_A(rho), np.eye(2)/2, atol=TOL)
    print(f"C2 correlation persists (ZZ = {corr:.3f}) with untouched marginal: PASS")

def test_C3_generic_states():
    for _ in range(100):
        rho = rand_joint_state()
        mA0 = marginal_A(rho)
        K = rand_kraus_qubit()
        assert np.allclose(marginal_A(apply_on_B(rho, K)), mA0, atol=1e-10)
    print("C3 marginal invariance on 100 generic joint states: PASS")

if __name__ == "__main__":
    test_C1_marginal_invariance_bell(); test_C2_correlation_persists()
    test_C3_generic_states()
    print("ALL PASS: correlation-without-influence toy verified.")
