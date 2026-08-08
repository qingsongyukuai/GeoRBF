# Issue 45: source-localized conflict witnesses

Issue: [#45](https://github.com/qingsongyukuai/GeoRBF/issues/45)

Failed fits diagnosed as `DirectInputConflict` or `InfeasibleProblem` now expose
one common `FitReport::conflict_witness()` audit surface. The witness identifies
the original caller `SourceId`s, records the canonical hard-relation
multipliers, and publishes the independently recomputed canonical residual and
strict separation margin. Relations and distinct source identities have stable
ordering; the API does not promise a globally minimum conflict set.

Direct equal-left-side conflicts and contradictory shared-level graph cycles
produce their witness before backend execution. The canonical hard recovery
graph also recognizes an exactly reconstructed left-side dependency whose
target is not exactly reconstructed, rebuilds the proof from original source
relations, and verifies that the source combination has exactly zero canonical
left side and a nonzero target before returning `DirectInputConflict`.

For general convex infeasibility, a backend termination is still only a
candidate. GeoRBF first recovers and normalizes the ray, recomputes stationarity,
dual-cone membership, strict separation, scaling round-trip, and complete
provenance against its immutable QP form, then repeats the certificate check
using only caller-owned canonical rows. Published residuals and separation are
the raw combination metrics for the published multipliers; scale-aware limits
are published separately. Derived nonnegativity rows therefore cannot leak into
the source witness. Only a validated ray produces an `InfeasibleProblem`
witness. Zero-multiplier constraints are omitted from the proof combination.
Every original source attached to an active compressed hard row is nevertheless
expanded through the hard recovery provenance, so aliases remain locatable
without being misrepresented as additional algebraic terms.

Invalid rays, limit terminations, representation failures, poor conditioning,
rank-decision gray zones, and other numerical failures do not receive a hard
conflict witness. Equality and QP conflict failures remain `FitFailure` values;
their reports contain no canonical acceptance, recovered hard-relation fields,
or candidate model.
