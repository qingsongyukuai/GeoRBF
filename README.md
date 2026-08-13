# GeoRBF

GeoRBF is a pure Rust implementation of the frozen Surfe non-visual modelling
core. It provides safe geological constraint types, nine radial kernels, five
model types, pure-Rust dense LU and quadratic solvers, and immutable scalar and
gradient evaluation. Normal Cargo builds do not use C++, Eigen, Qt, VTK,
CMake, bindgen, native BLAS, or the external Surfe reference.

The compatibility target is Surfe commit
`290dbe0ab344f4258a4935f05cad0f153f0f69a4`. GeoRBF intentionally excludes
GUI, visualization, `geo_builder`, file-format display helpers, and undefined
C++ behavior.

## Requirements

- Rust 1.82 or newer
- No system libraries or native toolchain for normal build, test, or use

The crate is currently marked `publish = false`; consume it from this
repository or a vendored source checkout.

## Quick start

The public lifecycle is `Builder` -> `fit()` -> immutable `FittedModel`:

```rust
use georbf::{Builder, ModelType, Point, RbfKernel};

let mut builder = Builder::new(ModelType::SingleSurface);
builder
    .set_rbf_kernel(RbfKernel::Cubic)
    .set_polynomial_order(1);

for [x, y, z] in [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0, 1.0, 0.0],
] {
    builder.add_interface_xyz(x, y, z, 0.0)?;
}
builder.add_planar_normal(0.5, 0.5, 0.0, 0.0, 0.0, 1.0)?;

let model = builder.fit()?;
let query = Point::new(0.25, 0.75, 0.5)?;
assert_eq!(model.evaluate_scalar(&query)?, 0.5);
assert_eq!(model.evaluate_gradient(&query)?, [0.0, 0.0, 1.0]);

# Ok::<(), Box<dyn std::error::Error>>(())
```

`FittedModel` owns a snapshot of the configuration and constraints. Editing
the builder after fitting does not mutate an earlier model. Fitted models are
`Send + Sync`; scalar and gradient batch evaluation preserves input order and
does not share mutable kernel state.

## Model types

`ModelType` exposes the five frozen model branches:

- `SingleSurface`: linear equality, ordinary inequality/QP, or restricted
  range depending on the supplied constraints and parameters.
- `LajaunieApproach`: exact-level reference points and same-level increments.
- `StratigraphicHorizons`: ordered horizons and lithostratigraphic relations.
- `ContinuousProperty`: the actually reachable frozen interface-value path.
- `VectorField`: planar-Hessian potential and gradient fitting.

Kernel names and model codes are deliberately strict Surfe-compatible values;
use the typed `RbfKernel` and `ModelType` enums where possible. Public failures
offer `surfe_category()` for stable programmatic classification. A `None`
category means GeoRBF safely rejected a source case that had no stable Surfe
exception, such as undefined or non-finite C++ behavior.

Restricted Single Surface, Lajaunie, and Stratigraphic public evaluation
follows the frozen reachable path and evaluates its Modified-Kernel field. The
separately migrated explicit reconstruction operation remains available as a
lower-level API, but frozen `ComputeInterpolant` did not call it automatically.
Greedy settings are retained, while the frozen public API's zero-round behavior
is reported honestly rather than presenting unreachable hooks as a feature.

## Build and verify

```text
cargo fmt --all --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
cargo test --doc --all-features
cargo build --release --all-features
python3 tools/audit_pure_rust.py
cargo package --locked
```

The checked-in global parity fixture is consumed without discovering or
starting a C++ oracle. Recreating the oracle evidence or rerunning the
cross-language performance comparison is an explicit maintainer operation and
requires the separately ignored frozen reference.

## Compatibility and release evidence

- [Compatibility and deliberate safety differences](docs/port/compatibility.md)
- [Per-module frozen source traceability](docs/port/source-traceability.md)
- [Global behavior parity report](docs/port/parity-report.md)
- [Same-machine performance report](docs/port/performance-report.md)
- [Final release audit](docs/port/release-audit.md)

The repository CI definition runs formatting, strict linting, all tests,
documentation, release build, dependency/native audits, and package creation
on Linux, macOS, and Windows.

## License

GeoRBF is MIT licensed. Rust translations and adaptations retain the frozen
Surfe Government of Canada MIT notice in [NOTICE](NOTICE); the project license
is in [LICENSE](LICENSE). Canada wordmarks, related graphics, Surfe third-party
submodules, and excluded visual code are not distributed.
