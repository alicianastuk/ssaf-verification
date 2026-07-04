"""SSAF verification suite, file 9: Phase-Circle Envelope Ordering
(toy-model sections; worked instance Delta_0 = pi/2, Delta_pi = pi/6, and the
non-closure certification under a shrunk-window regime).

Pass conditions verified:
  P1: the principal-value offset table matches the paper's twelve listed
      values exactly.
  P2: the membership predicate at (pi/2, pi/6) yields exactly the paper's
      mutual edges 1<->2, 3<->4, 1<->3, 2<->4 and no edge for (1,4),(2,3).
  P3: window membership is symmetric here (PV negation), so the listed
      edges are mutual, consistent with the toy's evaluation.
  P4 (regime change / non-closure): under a context with shrunk windows
      (pi/4, pi/12), the pair (1,2) is certified non-closed: neither the
      offset nor its antipode lies in the window. Failure is membership
      only; nothing is retuned to repair it.
  P5: the antipode map is an involution on PV space, and for an
      A-invariant window the two disjuncts of the closure predicate
      coincide (Lemma antipode invariance, conditional form).
"""
import numpy as np

PHI = {1: 0.0, 2: np.pi/3, 3: np.pi, 4: 4*np.pi/3}

def pv(x):
    """principal value in (-pi, pi]"""
    y = (x + np.pi) % (2*np.pi) - np.pi
    return np.pi if np.isclose(y, -np.pi) else y

def offset(i, j):
    return pv(PHI[j] - PHI[i])

def in_window(x, d0, dpi):
    return abs(x) <= d0 + 1e-12 or abs(x) >= np.pi - dpi - 1e-12

def antipode(x):
    return pv(x + np.pi)

PAPER_TABLE = {
    (1,2): np.pi/3,  (1,3): np.pi,      (1,4): -2*np.pi/3,
    (2,1): -np.pi/3, (2,3): 2*np.pi/3,  (2,4): np.pi,
    (3,1): np.pi,    (3,2): -2*np.pi/3, (3,4): np.pi/3,
    (4,1): 2*np.pi/3,(4,2): np.pi,      (4,3): -np.pi/3,
}

def test_P1_offset_table():
    for (i, j), val in PAPER_TABLE.items():
        assert np.isclose(offset(i, j), val), f"offset {i}{j} mismatch"
    print("P1 all twelve PV offsets match the paper: PASS")

def test_P2_edge_set():
    d0, dpi = np.pi/2, np.pi/6
    edges = {(i, j) for i in PHI for j in PHI if i != j
             and in_window(offset(i, j), d0, dpi)}
    expected = {(1,2),(2,1),(3,4),(4,3),(1,3),(3,1),(2,4),(4,2)}
    assert edges == expected, f"edge set {sorted(edges)} != paper's"
    for pair in [(1,4),(2,3)]:
        assert pair not in edges
    print("P2 directed participation edges match the worked instance: PASS")

def test_P3_mutuality():
    d0, dpi = np.pi/2, np.pi/6
    for i in PHI:
        for j in PHI:
            if i != j:
                assert in_window(offset(i, j), d0, dpi) == \
                       in_window(offset(j, i), d0, dpi)
    print("P3 membership symmetric under direction reversal here: PASS")

def test_P4_nonclosure_regime():
    d0, dpi = np.pi/4, np.pi/12     # shrunk-window context E'
    x = offset(1, 2)                 # pi/3
    assert not in_window(x, d0, dpi), "offset should fall outside"
    assert not in_window(antipode(x), d0, dpi), "antipode should too"
    print("P4 certified window non-closure under shrunk windows: PASS")

def test_P5_antipode_involution_and_invariance():
    for x in np.linspace(-np.pi + 1e-6, np.pi, 50):
        assert np.isclose(antipode(antipode(x)), pv(x))
    # A-invariant window: d0 = dpi = pi/4 makes W = near-0 u near-pi with
    # A mapping each lobe onto the other; the closure disjuncts coincide.
    d = np.pi/4
    for x in np.linspace(-np.pi + 1e-6, np.pi, 720):
        assert in_window(x, d, d) == in_window(antipode(x), d, d)
    print("P5 antipode involution; invariant-window disjuncts coincide: PASS")

if __name__ == "__main__":
    test_P1_offset_table(); test_P2_edge_set(); test_P3_mutuality()
    test_P4_nonclosure_regime(); test_P5_antipode_involution_and_invariance()
    print("ALL PASS: phase-circle envelope ordering verified.")
