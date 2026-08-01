# Risk spike 14: Cubic CPD KKT/QP/SOCP recovery

Issue: [#14](https://github.com/qingsongyukuai/GeoRBF/issues/14)

Evidence seams: T06, T07

Requirements: IR-001, IR-007, IR-011, IR-014, KER-006, KER-007, NUM-004, NUM-006, VAL-005, VAL-013

Probe: [`spikes/cubic-cpd-recovery`](../../spikes/cubic-cpd-recovery/README.md)

## Verdict

**The manufactured mathematical experiment succeeds.** One physical Canonical Problem IR built from the same ten generalized value/first-derivative functionals is recovered from an Equality augmented KKT, a reduced QP, and a reduced SOCP. The complete `Pi_1 = span{1,x,y,z}` is retained, the Cubic side condition and reduced strict positivity hold, all three routes agree on the canonical observables, and no ridge, jitter, polynomial deletion, automatic anchor, or kernel fallback is used.

This closes the mathematical risk tested by #14; it does not implement a product Canonical Problem IR, solver adapter, API, diagnostics system, or production dependency decision. The temporary pinned backend realization succeeds for this experiment only. The independent production dependency gate remains governed by #13 and its unresolved `faer` factor-workspace capacity evidence.

## Manufactured case

The case contains ten deterministic, single-support scalar functionals

\[
F_i(f)=a_i f(x_i)+b_i^\mathsf T\nabla f(x_i),
\]

comprising five value functionals, three coordinate derivatives, one value/derivative contraction, and one directional derivative. For Cubic `k(x,y)=||x-y||^3`, the probe assembles the symmetric generalized-functional pairing `K` from the analytic value, first, and mixed jets.

The physical canonical problem owns these functionals, eleven equality relations, one affine upper bound, one SOC, and one semantic latent. The equalities fix every functional value and the latent, so the upper bound and SOC are independently verified as redundant. Equality may elide both without changing the feasible set; QP retains the bound; SOCP retains the cone. The forms therefore derive from one canonical semantics rather than inventing truth-relative constraints below the lowering seam.

The polynomial pairing is constructed in the specified normalized coordinates and has shape `10 x 4`. Its physical/normalized coefficient map is explicit and reversible. A fixed nonzero physical coefficient vector satisfying `P^T c=0` manufactures functional values, the semantic latent, hard residual truth, FieldEnergy, and objective. Rank and reduced-positivity evidence use the specified separated acceptance/rejection thresholds and return a numerical decision gray zone between them.

## Mathematical evidence

| Property | Observed | Acceptance |
| --- | ---: | ---: |
| `rank(P)` | `4` | `4` |
| singular values of `P` | `2.5120, 1.8438, 1.3163, 1.2910` | full rank outside gray zone |
| `||P^T Z||_max` | `1.6653e-16` | side condition satisfied |
| normalized reduced symmetry defect | `3.7050e-16` | `<= 3.4106e-13` (`256 u d`) |
| smallest eigenvalue of symmetrized `Z^T K Z` | `2.3302` | strictly positive |
| affine reproduction error | `3.3307e-16` | `<= 1e-11` |
| Equality KKT inertia | `15 / 15 / 0` | `15 / 15 / 0` |
| Equality normalized backward error | `3.5034e-17` | `<= 1e-11` |
| cross-route canonical observable error | `3.5527e-15` | `<= 1e-8` |
| recovered FieldEnergy | `0.8188262034758015` | manufactured truth |
| recovered total objective | `0.4094131017379007` | `1/2 FieldEnergy` |

The Equality form is a `30 x 30` symmetric saddle system over 15 primal variables and 15 hard constraints. Its primal variables retain ten field coefficients, four normalized polynomial coefficients, and one semantic latent; recovery returns the polynomial in physical coordinates.

The QP and SOCP use a `10 -> 6` null-space map stored as `faer` Householder QR factors. `Z` is never retained as a dense matrix: the probe expands or projects transient vectors through the reflectors and only materializes backend-required reduced rows and the `6 x 6` Hessian. Cholesky and eigenvalue evidence are checked before Clarabel is called.

## Numerical-form and recovery evidence

All forms receive exactly eight rounds of GeoRBF-owned max-norm Ruiz scaling. Each factor is quantized to the nearest power of two, clipped per round and cumulatively, the KKT uses diagonal congruence, and the SOCP block is scaled as a unit.

| Evidence | Reduced QP | Reduced SOCP | Limit |
| --- | ---: | ---: | ---: |
| scaled primal | `3.9280e-12` | `3.3324e-12` | `1e-8` |
| scaled dual cone violation | `0` | `0` | `1e-8` |
| scaled stationarity | `1.5382e-11` | `1.5876e-11` | `1e-8` |
| scaled complementarity | `6.8326e-12` | `7.8546e-12` | `1e-8` |
| scaled relative gap | `2.0534e-11` | `2.1995e-11` | `1e-8` |
| reduction/scaling/polynomial round-trip | `1.3878e-17` | `1.3878e-17` | `1e-11` |
| physical slack equation / hard violation | `1.5345e-11` | `1.3853e-11` | `1e-8` |
| manufactured truth error | `1.7764e-15` | `9.9920e-16` | `1e-8` |

Recovery reconstructs field coefficients, the physical polynomial, the shared-level latent, all canonical residuals and all four canonical slacks, the CPD side condition, hard violation, FieldEnergy, and objective. Independently unscaled backend slacks are checked against `A_physical x + s = b` and the corresponding canonical slack. Those physical observables—not backend row layouts, pivots, dual choices, or iteration trajectories—are compared across routes.

## Failure evidence

Three deliberately damaged cases fail at the intended boundary:

- flattening all support/derivative `z` information gives `rank(P)=3` and is rejected before any solve;
- subtracting known negative curvature along one valid null-space vector makes the reduced pairing nonpositive and is rejected before Clarabel;
- deliberately corrupting the Householder recovery map for an otherwise acceptable QP candidate violates the side condition, canonical hard relations, and map round-trip; the shared recovery validator rejects it as a recovery-verification failure.

Every failure records `hidden_regularization_applied=false`. The probe does not attempt to rescue them with a smaller polynomial, ridge/jitter, an anchor, a kernel change, or relaxed tolerance.

## Mathematical conclusions versus backend boundary

The following conclusions are properties of the manufactured algebra and independently recomputed observables: complete polynomial rank; `range(Z)=ker(P^T)` to floating-point evidence; reduced positive definiteness; affine reproduction; expected KKT inertia; equality feasibility; physical side condition; common recovered FieldEnergy/objective; and visibility of the three damaged cases.

The following conclusions are only facts about the temporary implementation: `faer = 0.24.4` successfully realizes the SVD, Householder QR, Cholesky/eigendecomposition, and LBLT solve used here; `clarabel = 0.11.1` with `qdldl`, one requested thread, and its fixed internal equilibration/regularization profile returns acceptable candidates for these reduced QP/SOCP instances. No statement is made about other versions, platform-wide behavior, capacity, retry behavior, general conditioning, or production eligibility.

## Replay

From `spikes/cubic-cpd-recovery` under Rust 1.85.0:

```text
cargo check --locked --all-targets
cargo test --locked --all-targets
cargo run --locked
```

The local replay passed all ten T06/T07 integration tests and emitted the evidence above. `.github/workflows/risk-spike-14.yml` repeats the locked replay on Ubuntu; a configured workflow is not itself remote-run proof.
