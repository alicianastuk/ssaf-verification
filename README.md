# SSAF Verification Suite

Executable verification of the finite toy constructions in the SSAF
paper's toy-model sections. Every file `test_NN_*.py` encodes one toy's pass conditions as
assertions; `python3 run_all.py` runs everything and exits nonzero on any
failure.

Requires: Python 3, numpy.

Verified against: SSAF manuscript v1.0 (December 2025, revised July
2026). Formal results are cited by label name in each test header, so
references are stable under section renumbering.

## Coverage

| File | Toy | Status |
|---|---|---|
| test_01_threshold_shadow.py | Declared Compatibility as a Threshold Shadow of Graded Access | PASS |
| test_02_ns1_safe_gksl.py | NS1-Safe Single-Qubit GKSL Semigroup | PASS |
| test_03_joint_admissibility.py | Two-Qubit Joint Admissibility + General Form (PPT, strict inclusion) | PASS |
| test_04_correlation_without_influence.py | Correlation Without Influence Under Local CPTP Updating | PASS |
| test_05_ordering_no_completion.py | Non-Total Ordering + Admissibility Entailment (no-completion at both levels) | PASS |
| test_06_sequencing_no_retro.py | Finite-Dim Sequencing Without Retro-Addressability | PASS |
| test_07_coherence_basin.py | Coherence Basin Depth Functional Contraction | PASS |
| test_08_charge_sectors.py | Connectivity Strictly Coarser via Charge Sectors | PASS |
| test_09_phase_circle.py | Phase-Circle Envelope Ordering (worked instance + non-closure regime) | PASS |
| test_10_acceleration_reconstruction.py | Late-Time Acceleration Without Dark Energy (closed-form match + null control) | PASS |

## Conventions

- Each test header cites the toy and the formal results it checks.
- Assertions mirror the paper's stated pass conditions; no test asserts
  anything the paper does not claim.
- Randomized checks use fixed seeds for reproducibility.
