# Working note: Compat_E and Γ_E as threshold shadows of a graded access functional

**Status: integrated into the SSAF manuscript v1.0. The construction below
appears as the toy "Declared Compatibility as a Threshold Shadow of Graded
Access", open question (v) was updated to partially resolved, and the
coherence inclusion corollary appears in the primitives section.**
**Numerically verified: `test_01_threshold_shadow.py` (both regimes pass).**

## What this attacks

Open question (v) asks whether the declared structural objects can be recovered
from a graded overlap or operational-distinguishability structure under E. The
paper's own candidate phrasing: Compat_E may be recoverable as the threshold
shadow of a graded compatibility functional, with Γ_E, N_E, µ_E recoverable
from the same graded structure.

This note constructs such a functional from E-accessible operations only, and
proves that its threshold shadows reproduce exactly the (Compat_E, Γ_E) pair
declared in the charge-sector toy. The toy thereby becomes the unit test for
the grounding.

## Definition (graded access functional)

Fix E. Let S_E ⊆ B(H) be the E-accessible operator system: *-closed, unital
(I ∈ S_E), with convex unit ball. (In the charge-sector toy, S_E is the
block algebra B(H₊) ⊕ B(H₋); the general definition needs only the operator-
system structure.)

For ρ, σ in the admissible set, define

    g_E(ρ, σ) := sup { |⟨φ| A |ψ⟩| : A ∈ S_E, ‖A‖ ≤ 1,
                        ψ ∈ supp ρ, φ ∈ supp σ, ‖ψ‖ = ‖φ‖ = 1 }.

g_E depends only on (supp ρ, supp σ, S_E). All free data — S_E and the
thresholds below — are fixed by E.

## Lemma 1 (elementary properties)

(i) g_E ∈ [0, 1] (Cauchy–Schwarz with ‖A‖ ≤ 1).
(ii) g_E(ρ, σ) = g_E(σ, ρ) (*-closure: |⟨φ|A|ψ⟩| = |⟨ψ|A†|φ⟩|).
(iii) g_E(ρ, ρ) = 1 (unitality: take A = I, φ = ψ).

No transitivity, metric, or topological structure is claimed — deliberately
matching the paper's statement that Compat_E assumes none of these.

## Theorem (two-threshold recovery, block case)

Let S_E = B(H₊) ⊕ B(H₋) and let A(E) be the definite-charge admissible set of
the charge-sector toy (supp ρ contained in a single sector for every
admissible ρ). Then

    g_E(ρ, σ) = 1   if sec(ρ) = sec(σ),
    g_E(ρ, σ) = 0   if sec(ρ) ≠ sec(σ).

*Proof.* Same sector k: both supports lie in H_k; the partial isometry
|φ⟩⟨ψ| belongs to B(H_k) ⊆ S_E, has norm 1, and achieves |⟨φ|φ⟩| · ‖ψ‖ = 1.
Different sectors: every A ∈ S_E is block-diagonal, so ⟨φ|A|ψ⟩ picks out the
off-block compression of A, which vanishes. ∎

**Consequence.** Choosing thresholds θ_c = 0 and any θ_γ ∈ (0, 1]:

    { (ρ,σ) : g_E ≥ θ_c }  =  total relation        =  the toy's Compat_E,
    { (ρ,σ) : g_E ≥ θ_γ }  =  same-sector relation  =  the toy's R_E (Γ_E edge set).

The toy's two declared primitives are the two threshold shadows of one
operational functional. Declared structure is recovered, not posited.

## Corollary (coherence inclusion becomes a theorem)

For any θ_γ ≥ θ_c, superlevel sets nest: {g_E ≥ θ_γ} ⊆ {g_E ≥ θ_c}, i.e.
Γ_E ⊆ Compat_E automatically. In any two-threshold grounding, the
compatibility coherence axiom added in the primitives refactor is not an
independent assumption but a consequence of threshold ordering.

## Proposition (the functional is genuinely graded)

Let S_ε := { D + C : D block-diagonal, C cross-block, ‖D + C‖ ≤ 1, ‖C‖ ≤ ε }
(leaky superselection: cross-sector access of strength ε). Then for
definite-sector pairs in different sectors, g_E = ε exactly (achieved by
C = ε|φ⟩⟨ψ|), while same-sector pairs retain g_E = 1.

So thresholds are substantive, not decorative: under degraded sector
protection, θ ≤ ε licenses cross-sector compatibility and θ > ε enforces
separation. Numerically confirmed for ε ∈ {0, 0.1, 0.3, 0.7}.

## NS1 compliance

(S_E, θ_c, θ_γ) are fixed by E. The arguments ρ, σ enter g_E only through
their supports as the relata of the predicate being evaluated; the induced
relations {g_E ≥ θ} are state-independent context-indexed predicates, fixed
once by (S_E, θ). Evaluating a fixed predicate at a pair is not retuning,
in the same sense already used for the declared Compat_E.

## Honest scope (to be stated wherever this is inserted)

1. **What is achieved:** Compat_E and Γ_E, in the charge-sector instantiation,
   are derivable as threshold shadows of a single functional built from
   E-accessible operations. This resolves the *operational-access half* of
   open question (v) in the affirmative, by construction, for this class of
   models. N_E admits a natural candidate from the same object
   (superlevel neighbourhoods N_E(ρ) = {σ : g_E(ρ,σ) ≥ τ}), offered as a
   remark, not developed.
2. **What is not achieved:** the construction lives on B(H). The
   sub-inner-product half of question (v) — grounding in primitives strictly
   weaker than inner-product geometry — remains open, exactly as the
   revised item (v) already states. µ_E is not recovered. ~E is deliberately
   untouched (g_E sees supports, not states), preserving the separation
   between linkage primitives and the indistinguishability quotient.
3. **Sufficiency, not uniqueness:** this exhibits *a* grounding; it does not
   claim the grounding is unique or canonical.

## Proposed placement (for discussion)

- New Appendix D toy: "Axiom Unit Test: Declared Compatibility as a Threshold
  Shadow of Graded Access" — definition, Lemma 1, Theorem, Corollary,
  Proposition, NS1 remark, scope. Cross-wired to the charge-sector toy
  (its structure is the test target) and from the main-body grounding
  discussion (the passage that currently says "This paper does not establish
  such a derivation" would gain: "a finite-dimensional instance of exactly
  this derivation is worked in Toy~\ref{...}").
- Open question (v): updated from fully open to partially resolved, with the
  sub-inner-product residue stated as what remains.
- Axiom dependency table: no change required (this grounds primitives, not
  axioms), though the row-IV unit-test cell could optionally gain the new toy.
