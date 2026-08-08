# Issue 41: Precision-rescue oracle and double-double arithmetic

Issue: [#41](https://github.com/qingsongyukuai/GeoRBF/issues/41)

Evidence seam: the parent specification's crate-private precision conformance
seam. No matrix or precision-control API is public.

## Independent oracle identity

`validation/oracle/precision-rescue-v1` contains declarative cases, a
dependency-free CPython generator, the generated fixture, and a source
manifest. The generator uses only Python `decimal` at 160 working decimal
digits with ties-to-even rounding; it imports neither GeoRBF production code
nor faer, Clarabel, LAPACK, Surfe, or a native multiprecision runtime. The
manifest records SHA-256 identities for the declaration, generator, and
fixture. Rust tests separately pin the bytes of all three and the manifest, so
a changed input, oracle implementation, result, or provenance fails
deterministically.

The generated fixture retains 120 significant decimal digits and a canonical
two-word binary64 projection for every result. Its cases cover addition,
subtraction, multiplication, division, square root, a multi-term anisotropic
Cubic generalized-functional pairing, and symmetric Schur accumulation.

## Pure-Rust arithmetic

`precision_rescue` implements a normalized unevaluated sum of two binary64
words. Error-free `two_sum` and fused `two_product` transforms support the
basic operations, while division and square root use deterministic residual
correction. There is no new Cargo dependency and therefore no native
multiprecision runtime dependency. The conformance tests compare the combined
two-word result with independent oracle values at an error scale proportional
to `f64::EPSILON²`, exercising the approximately 106-bit representation rather
than comparing only its leading word.

The same module recomputes Cubic value/gradient generalized-functional
pairings from canonical f64 inputs entirely in double-double arithmetic. The
Cubic jet scales displacement before its quadratic form, preserving finite
mixed derivatives even when the unscaled squared distance would underflow. It
also provides a deterministic symmetric Schur entry accumulation
`pairing - sum(left[k] * right[k])`, suitable for the bounded rescue module in
issue #42.

## Distinguishable modal conclusions

The Schur corpus contains independent truths for a strictly positive
`2^-100` residual, exact algebraic zero, negative `-2^-100` curvature, and a
general two-product accumulation. Tests compare each value with the oracle and
independently assert the three signs, so a small positive mode, true zero, and
negative curvature cannot collapse into the same conclusion.
