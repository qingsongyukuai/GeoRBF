# Precision-rescue oracle v1

This corpus is the external truth source for GeoRBF issue #41. The checked-in
declarations are evaluated by CPython's independent `decimal` stack at 160
decimal digits with ties-to-even rounding. The generator imports no GeoRBF
code, numerical backend, or native multiprecision library.

Regenerate with `python3 validation/oracle/precision-rescue-v1/generate.py`.
The source manifest fixes the declaration, generator, and fixture byte hashes;
Rust conformance tests additionally pin the tracked manifest and fixture
identities.

The corpus covers double-double addition, subtraction, multiplication,
division and square root; a multi-term anisotropic Cubic canonical pairing;
and symmetric Schur entries representing a strictly small positive mode, a
true algebraic zero, negative curvature, and a general off-diagonal-style
accumulation.
