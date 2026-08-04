#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]

pub mod advanced;
#[allow(dead_code)]
mod capacity;
#[allow(dead_code)]
mod clarabel_backend;
#[allow(dead_code)]
mod cubic;
#[allow(dead_code)]
mod cubic_equality;
#[allow(dead_code)]
mod cubic_execution;
mod cubic_solver_form;
pub mod diagnostics;
#[allow(dead_code)]
mod faer_backend;
pub mod fit;
#[allow(dead_code)]
mod functional;
pub mod geometry;
pub mod kernel;
#[allow(dead_code)]
mod kkt;
#[allow(dead_code)]
mod math;
pub mod model;
#[allow(dead_code)]
mod numerical;
pub mod observation;
pub mod problem;
pub mod relation;

pub use geometry::{Point3, Vector3};
pub use model::SolvedModel;
pub use problem::{GroupId, ProblemBuilder, ProblemSnapshot, SourceId};

#[cfg(test)]
mod oracle_fixture;
