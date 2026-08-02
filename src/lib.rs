//! GeoRBF product crate.
//!
//! Version 0.1.0 establishes the internal Cubic generalized-functional Equality
//! core, versioned numerical policy, bounded KKT attempts, and physical Recover
//! and Verify boundary. Domain inputs and fitting APIs are intentionally not
//! public until they are complete.

#![forbid(unsafe_code)]

// These modules are the v0.1.0 product spine. They remain crate-internal until
// a later milestone can expose the complete domain-facing fit workflow.
#[allow(dead_code)]
mod capacity;
#[allow(dead_code)]
mod cubic;
#[allow(dead_code)]
mod cubic_equality;
#[allow(dead_code)]
mod faer_backend;
#[allow(dead_code)]
mod functional;
#[allow(dead_code)]
mod kkt;
#[allow(dead_code)]
mod math;
#[allow(dead_code)]
mod numerical;

#[cfg(test)]
mod oracle_fixture;
