"""SSAF verification suite, file 5: Non-Total Ordering Without Dynamics +
Admissibility Entailment at the quotient level (Appendix D ordering pair;
Lemmas toy_no_completion and toy_no_completion_entailment).

The toy's declared base constraints (Eq. toy_prereq_base):
    A <| C,  B <| C,  C <| D,  B <| F
with prec_E the strict part of the reflexive-transitive closure.

Pass conditions verified:
  O1: prec_E is irreflexive, transitive, and acyclic.
  O2: prec_E is not total: incomparable pairs exist (e.g. A,B and D,F).
  O3 (no completion, configuration level): the generated prerequisite set
      contains neither direction for the incomparable pairs -- the ordering
      claim is undefined, and any assignment would add structure.
  O4 (no completion, entailment/quotient level): running the identical
      closure on the addressability-class labels yields the same
      incomparability -- passing to the quotient does not complete the
      order (stability under coarse-graining).
"""
from itertools import product

CLASSES = ["A", "B", "C", "D", "F"]
BASE = {("A", "C"), ("B", "C"), ("C", "D"), ("B", "F")}

def transitive_closure(base, elems):
    rel = set(base)
    changed = True
    while changed:
        changed = False
        for (x, y), (u, v) in product(list(rel), list(rel)):
            if y == u and (x, v) not in rel:
                rel.add((x, v)); changed = True
    return rel

PREC = transitive_closure(BASE, CLASSES)   # strict relation (no reflexive pairs added)

def test_O1_strict_partial_order():
    assert all((x, x) not in PREC for x in CLASSES), "irreflexivity failed"
    for (x, y), (u, v) in product(PREC, PREC):
        if y == u:
            assert (x, v) in PREC, "transitivity failed"
    assert all(not ((x, y) in PREC and (y, x) in PREC)
               for x in CLASSES for y in CLASSES), "cycle found"
    print(f"O1 strict partial order (|prec| = {len(PREC)}): PASS")

def test_O2_not_total():
    incomparable = [(x, y) for x in CLASSES for y in CLASSES if x < y
                    and (x, y) not in PREC and (y, x) not in PREC]
    assert ("A", "B") in incomparable and ("D", "F") in incomparable
    print(f"O2 incomparable pairs exist {incomparable}: PASS")

def test_O3_no_completion_config_level():
    # the ONLY licensed facts are the generated set; both directions absent
    for pair in [("A", "B"), ("D", "F")]:
        x, y = pair
        assert (x, y) not in PREC and (y, x) not in PREC, \
            "an order was entailed where the toy declares undefinedness"
    # adding either direction is consistent but is EXTRA structure:
    # verify it is not entailed (i.e., closure without it never produces it)
    print("O3 no internal completion at configuration level: PASS")

def test_O4_entailment_level_no_completion():
    # entailment toy: lift to CLAIMS. "K <| L" reads: licensing L is
    # structurally committed to licensing K. So at the claim level,
    # Adm(L) entails Adm(K) whenever K prec L: the entailment relation is
    # the CONVERSE of prec, closed reflexively-transitively.
    entails = {(y, x) for (x, y) in PREC} | {(x, x) for x in CLASSES}
    # closure sanity: transitive by construction (converse of transitive).
    for (x, y), (u, v) in product(entails, entails):
        if y == u:
            assert (x, v) in entails, "entailment closure not transitive"
    # licensing D commits to C, A, B (prerequisite chain):
    for k in ["C", "A", "B"]:
        assert ("D", k) in entails
    # no completion at the claim level: for incomparable pairs, neither
    # class-admissibility claim entails the other.
    for x, y in [("A", "B"), ("D", "F")]:
        assert (x, y) not in entails and (y, x) not in entails, \
            "claim-level entailment appeared where the toy declares undefinedness"
    print("O4 no entailment completion at the claim level: PASS")

if __name__ == "__main__":
    test_O1_strict_partial_order(); test_O2_not_total()
    test_O3_no_completion_config_level(); test_O4_entailment_level_no_completion()
    print("ALL PASS: ordering pair verified.")
