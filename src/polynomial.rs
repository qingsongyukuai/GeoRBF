//! Polynomial values and first derivatives in frozen Surfe term order.
//!
//! Sources:
//! `surfe_lib/basis.{h,cpp}@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! (`Polynomial_Basis`, `Poly_Zero`, `Poly_First`, and `Poly_Second`).

use crate::Point;

#[path = "polynomial/lagrangian.rs"]
mod lagrangian;

pub use lagrangian::LagrangianPolynomialBasis;

/// Polynomial orders implemented by frozen Surfe.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PolynomialOrder {
    Zero,
    First,
    Second,
}

/// A complete or constant-truncated polynomial basis.
///
/// Complete term order is:
///
/// - order zero: `[1]`;
/// - order one: `[x, y, z, 1]`;
/// - order two: `[x², y², z², xy, xz, yz, x, y, z, 1]`.
///
/// Frozen Surfe's truncated form removes only the final constant term. This
/// yields lengths zero, three, and nine for orders zero, one, and two. The
/// Lajaunie and Stratigraphic Surfaces polynomial blocks use that form.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PolynomialBasis {
    order: PolynomialOrder,
    truncated: bool,
}

impl PolynomialBasis {
    /// Construct a complete polynomial basis.
    pub const fn complete(order: PolynomialOrder) -> Self {
        Self {
            order,
            truncated: false,
        }
    }

    /// Construct Surfe's constant-truncated polynomial basis.
    pub const fn truncated(order: PolynomialOrder) -> Self {
        Self {
            order,
            truncated: true,
        }
    }

    pub const fn order(self) -> PolynomialOrder {
        self.order
    }

    pub const fn is_truncated(self) -> bool {
        self.truncated
    }

    /// Number of terms returned by values and each derivative.
    pub const fn term_count(self) -> usize {
        match (self.order, self.truncated) {
            (PolynomialOrder::Zero, false) => 1,
            (PolynomialOrder::Zero, true) => 0,
            (PolynomialOrder::First, false) => 4,
            (PolynomialOrder::First, true) => 3,
            (PolynomialOrder::Second, false) => 10,
            (PolynomialOrder::Second, true) => 9,
        }
    }

    /// Evaluate polynomial terms at a point.
    ///
    /// The fourth Surfe coordinate `c` is intentionally unused.
    pub fn values(self, point: &Point) -> Vec<f64> {
        let mut values = match self.order {
            PolynomialOrder::Zero => vec![1.0],
            PolynomialOrder::First => vec![point.x(), point.y(), point.z(), 1.0],
            PolynomialOrder::Second => vec![
                point.x() * point.x(),
                point.y() * point.y(),
                point.z() * point.z(),
                point.x() * point.y(),
                point.x() * point.z(),
                point.y() * point.z(),
                point.x(),
                point.y(),
                point.z(),
                1.0,
            ],
        };
        values.truncate(self.term_count());
        values
    }

    /// Evaluate first derivatives with respect to `x` in matching term order.
    pub fn dx(self, point: &Point) -> Vec<f64> {
        let mut values = match self.order {
            PolynomialOrder::Zero => vec![0.0],
            PolynomialOrder::First => vec![1.0, 0.0, 0.0, 0.0],
            PolynomialOrder::Second => vec![
                2.0 * point.x(),
                0.0,
                0.0,
                point.y(),
                point.z(),
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
            ],
        };
        values.truncate(self.term_count());
        values
    }

    /// Evaluate first derivatives with respect to `y` in matching term order.
    pub fn dy(self, point: &Point) -> Vec<f64> {
        let mut values = match self.order {
            PolynomialOrder::Zero => vec![0.0],
            PolynomialOrder::First => vec![0.0, 1.0, 0.0, 0.0],
            PolynomialOrder::Second => vec![
                0.0,
                2.0 * point.y(),
                0.0,
                point.x(),
                0.0,
                point.z(),
                0.0,
                1.0,
                0.0,
                0.0,
            ],
        };
        values.truncate(self.term_count());
        values
    }

    /// Evaluate first derivatives with respect to `z` in matching term order.
    pub fn dz(self, point: &Point) -> Vec<f64> {
        let mut values = match self.order {
            PolynomialOrder::Zero => vec![0.0],
            PolynomialOrder::First => vec![0.0, 0.0, 1.0, 0.0],
            PolynomialOrder::Second => vec![
                0.0,
                0.0,
                2.0 * point.z(),
                0.0,
                point.x(),
                point.y(),
                0.0,
                0.0,
                1.0,
                0.0,
            ],
        };
        values.truncate(self.term_count());
        values
    }
}
