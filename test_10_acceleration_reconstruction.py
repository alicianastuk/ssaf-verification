"""SSAF verification suite, file 10: Apparent Late-Time Acceleration
Without a Dark-Energy Primitive (toy-model sections; Eqs. a_tt_general,
a_tt_linear_tau; reconstructed diagnostics H(t), q(t)).

Pass conditions verified, mirroring the toy:
  A1: the ordering-indexed baseline is exactly linear, a(tau) = a0 + v*tau,
      so a''(tau) = 0 identically.
  A2: the clock factor alpha(tau) is an E-indexed declared function with
      alpha'(tau) < 0 (slowing operational clock); it never sees a state.
  A3: the reconstructed second derivative d2a/dt2 computed numerically on
      the t-grid is positive everywhere, and matches the paper's closed
      form  d2a/dt2 = -(alpha'/alpha^3) v  pointwise.
  A4: the reconstructed deceleration parameter q(t) = -(a * a_tt)/(a_t)^2
      is negative (acceleration) along the chain.
  A5: with a constant clock (alpha' = 0) the reconstructed acceleration
      vanishes: the effect is entirely a property of the reconstruction
      map, exactly as the toy claims.
"""
import numpy as np

a0, v = 1.0, 0.05
alpha0, k = 1.0, 0.15          # alpha(tau) = alpha0 * exp(-k tau), E-fixed
TAU = np.linspace(0.0, 20.0, 200001)

def alpha(tau):    return alpha0*np.exp(-k*tau)
def alpha_p(tau):  return -k*alpha0*np.exp(-k*tau)
def a_of_tau(tau): return a0 + v*tau

def reconstruct():
    t = np.concatenate([[0.0], np.cumsum(np.diff(TAU)*0.5*(alpha(TAU[1:]) + alpha(TAU[:-1])))])
    a = a_of_tau(TAU)
    a_t = np.gradient(a, t)
    a_tt = np.gradient(a_t, t)
    return t, a, a_t, a_tt

def test_A1_baseline_linear():
    a = a_of_tau(TAU)
    second = np.diff(a, 2)/np.diff(TAU)[:-1]**2
    assert np.max(np.abs(second)) < 1e-5, "a''(tau) must vanish to numerical precision"
    print("A1 baseline a''(tau) = 0: PASS")

def test_A2_clock_e_indexed_and_slowing():
    assert np.all(alpha_p(TAU) < 0), "clock factor must be strictly slowing"
    # structural NS1: alpha is a closed-form function of the ordering index
    # alone; nothing in this module constructs alpha from a state.
    print("A2 alpha E-indexed, alpha' < 0 everywhere: PASS")

def test_A3_reconstructed_acceleration_matches_closed_form():
    t, a, a_t, a_tt = reconstruct()
    inner = slice(1000, -1000)   # avoid gradient edge artifacts
    assert np.all(a_tt[inner] > 0), "reconstructed d2a/dt2 must be positive"
    closed = -(alpha_p(TAU)/alpha(TAU)**3)*v
    rel = np.abs(a_tt[inner] - closed[inner])/closed[inner]
    assert np.max(rel) < 1e-3, f"closed-form mismatch, max rel err {np.max(rel):.2e}"
    print(f"A3 d2a/dt2 > 0 and matches -(alpha'/alpha^3)v "
          f"(max rel err {np.max(rel):.1e}): PASS")

def test_A4_deceleration_parameter_negative():
    t, a, a_t, a_tt = reconstruct()
    inner = slice(1000, -1000)
    q = -(a[inner]*a_tt[inner])/(a_t[inner]**2)
    assert np.all(q < 0), "q(t) must be negative (acceleration)"
    print(f"A4 q(t) < 0 along the chain (mean q = {q.mean():.3f}): PASS")

def test_A5_constant_clock_null_control():
    tau = TAU
    t = alpha0*tau               # alpha' = 0 reconstruction
    a = a_of_tau(tau)
    a_tt = np.gradient(np.gradient(a, t), t)
    assert np.max(np.abs(a_tt[1000:-1000])) < 1e-5, \
        "constant clock must reconstruct zero acceleration"
    print("A5 constant-clock control reconstructs zero acceleration: PASS")

if __name__ == "__main__":
    test_A1_baseline_linear(); test_A2_clock_e_indexed_and_slowing()
    test_A3_reconstructed_acceleration_matches_closed_form()
    test_A4_deceleration_parameter_negative(); test_A5_constant_clock_null_control()
    print("ALL PASS: acceleration-as-reconstruction verified.")
