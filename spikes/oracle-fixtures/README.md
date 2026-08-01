# Independent high-precision oracle fixture spike

This disposable Python experiment is the T13 evidence seam for [GeoRBF issue #15](https://github.com/qingsongyukuai/GeoRBF/issues/15). It does not import, link, execute, or derive values from GeoRBF production code or a production numerical backend. The generator uses only CPython's independent `decimal` stack at 120 decimal digits with round-to-nearest, ties-to-even.

## Representative cases

The deliberately small declarative corpus proves the pipeline rather than freezing a production corpus:

- `cubic-general-jet` evaluates `k(x,y)=r_M^3`, both first jets, and the mixed jet at a non-origin point under a non-identity determinant-one metric;
- `cubic-origin` exercises the analytic origin branch and emits canonical positive zero for every Cubic jet component;
- `cubic-generalized-functional` contracts value and first-derivative coefficients on both kernel arguments and records manufactured affine-field functional truth.

Every result scalar stores 110 significant decimal digits and the exact hexadecimal encoding of its correctly rounded `f64`. The CaseId is the SHA-256 identity of the canonical declarative case. Each fixture separately hashes its result and semantic content; the manifest hashes every declaration and fixture and records provenance.

## Locked environment

| Boundary | Fixed value |
| --- | --- |
| Interpreter | CPython 3.12.3 |
| Arithmetic | standard-library `decimal`, 120 digits, `ROUND_HALF_EVEN` |
| Third-party closure | empty, recorded by `requirements.lock` |
| OCI image | `python:3.12.3-slim-bookworm@sha256:afc139a0a640942491ec481ad8dda10f2c5b753f5c969393b12480155fe15a63` |

The OCI reference is a digest, not a mutable tag. `manifest.json` also records hashes of the generator, independent verifier, and lockfile. It intentionally excludes timestamps, host paths, locale, temporary directories, and VCS state.

## Replay

Under exact CPython 3.12.3:

```text
python -m pip install --require-hashes -r requirements.lock
python -m unittest discover -s tests -v
python generate.py --check
python generate.py --check
python verify.py
```

The independent consumer validates manifest coverage, source/lock/case/fixture/content/output hashes, stable CaseIds, canonical input identity, environment pins, and decimal-to-hexadecimal `f64` round-trip. The tests copy and tamper an artifact and require a stable verification failure.

The pinned container can replay the same seam with:

```text
docker build -t georbf-oracle-spike .
docker run --rm georbf-oracle-spike
```

The repository workflow also runs the checked-out artifacts twice in the exact digest-pinned base image, so the tag used in the `Dockerfile` cannot redirect the execution.

## Boundary of the conclusion

This spike proves a deterministic, independently consumable fixture pipeline for the three representative Cubic cases. It does not freeze a production kernel, CPD, optimization, diagnostic, tolerance, or ValidationProfile corpus; admit an oracle dependency to product code; or establish that future corpus formulas are correct merely because the pipeline is reproducible.
