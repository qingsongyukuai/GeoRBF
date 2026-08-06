# Issue 40: Verified energy-orthonormal quotient basis

Issue: [#40](https://github.com/qingsongyukuai/GeoRBF/issues/40)

Primary evidence seam: the existing public `ProblemBuilder` to snapshot to
`fit` to diagnostics/model-query workflow. A crate-private conformance seam is
used only for LLT pivot intervals and cancellation cases that cannot be
manufactured through a public matrix API.

## One unregularized quotient LLT

After the implicit Householder congruence from issue 39, the representation
forms the complete solver-facing quotient FieldEnergy Gram matrix. It submits
that matrix once to faer's LLT with both dynamic-regularization parameters set
to zero. A successful LLT does not trigger quotient RRQR, SVD, eigenvalue
truncation, ridge, or jitter. All quotient dimensions remain present in the
factor and in the solver form.

The retained lower factor is accepted only when

\[
\eta_G =
\frac{\lVert G-\hat L\hat L^T\rVert_\infty}
{\lVert G\rVert_\infty+\lVert |\hat L||\hat L|^T\rVert_\infty}
\le 10^{-11}.
\]

Each pivot is also recomputed from the retained lower-factor row with
outward-rounded square, sum, and subtraction bounds. Every interval lower
bound must be strictly positive. When LLT rejects, the failed pivot interval is
recomputed from its row. An interval that crosses zero is classified as
requiring precision rescue; an interval wholly at or below zero is a numerical
non-positive-factorization failure. Neither is reported as rank deficient. The
pending-rescue marker is not promoted to the terminal Numerical Decision Gray
Zone before the bounded double-double rescue supplied by issue 42.

## Reversible energy coordinates

For `G = L L^T`, solver coordinates `u` recover Householder quotient
coordinates through `z = L^{-T}u`, then canonical representer coefficients
through the existing implicit `Q2` expansion. Canonical functional responses
use the contragredient map `L^{-1}`. Mapping recovered canonical coefficients
back through `Q2^T` and `L^T` reconstructs `u`.

The solver-independent quotient Hessian is therefore the exact identity
matrix. Before publication, the representation verifies the complete
`L^{-1} G L^{-T}` identity defect, the transformed complete-Pi1 side condition,
a full-support solver-coordinate recovery probe, and round trips for every
canonical representer response. Each defect has a `1e-11` policy limit.

## Auditable evidence

Successful public Cubic analysis exposes the actual quotient and retained-mode
dimensions, truncated-mode count, LLT and post-success full-spectrum counts,
`eta_G`, every outward-rounded pivot interval, the energy/side/recovery/response
defects, and explicit booleans for kernel ridge, Gram jitter, and mode
truncation. The normal successful path records one LLT, zero quotient
full-spectrum analyses, zero truncated modes, and no problem regularization.

The public close-support regression fits and queries two distinct nearby
supports, checks that both quotient modes remain, and records a small but
strictly positive pivot ratio. The internal conformance corpus separately
proves that a reliably positive diagonal mode is retained even at a `1e-30`
relative scale, while a cancellation-ambiguous pivot is reserved for precision
rescue.
