//! Surfe-compatible geological constraint value types and orientation math.
//!
//! Sources:
//! - `surfe_lib/modelling_input.h@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! - `surfe_lib/modelling_input.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`
//! - `surfe_lib/surfe_api.cpp@290dbe0ab344f4258a4935f05cad0f153f0f69a4`

mod grouping;

pub use crate::ordering::CollocationRemoval;
pub use grouping::InterfaceGrouping;

use crate::{
    geometry::{is_zero_vector, require_finite},
    ConstraintError, Point, DEGREES_TO_RADIANS, RADIANS_TO_DEGREES,
};

/// Surfe's two documented polarity values.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[repr(i32)]
pub enum Polarity {
    Upright = 0,
    Overturned = 1,
}

impl Polarity {
    pub const fn surfe_code(self) -> i32 {
        self as i32
    }
}

impl TryFrom<i32> for Polarity {
    type Error = ConstraintError;

    fn try_from(code: i32) -> Result<Self, Self::Error> {
        match code {
            0 => Ok(Self::Upright),
            1 => Ok(Self::Overturned),
            _ => Err(ConstraintError::InvalidPolarity),
        }
    }
}

/// An interface point with an exact level and optional symmetric level bounds.
#[derive(Clone, Debug)]
pub struct Interface {
    point: Point,
    level: f64,
    level_bounds: [f64; 2],
}

impl Interface {
    pub fn new(x: f64, y: f64, z: f64, level: f64) -> Result<Self, ConstraintError> {
        Self::with_c(x, y, z, level, 0.0)
    }

    pub fn with_c(x: f64, y: f64, z: f64, level: f64, c: f64) -> Result<Self, ConstraintError> {
        require_finite(&[level])?;
        Ok(Self {
            point: Point::with_c(x, y, z, c)?,
            level,
            level_bounds: [0.0, 0.0],
        })
    }

    pub const fn point(&self) -> &Point {
        &self.point
    }

    pub const fn level(&self) -> f64 {
        self.level
    }

    pub const fn level_bounds(&self) -> [f64; 2] {
        self.level_bounds
    }

    pub const fn level_lower_bound(&self) -> f64 {
        self.level_bounds[0]
    }

    pub const fn level_upper_bound(&self) -> f64 {
        self.level_bounds[1]
    }

    pub fn set_level_bounds(&mut self, uncertainty: f64) -> Result<(), ConstraintError> {
        require_finite(&[uncertainty])?;
        self.level_bounds = [-uncertainty, uncertainty];
        Ok(())
    }
}

/// An inequality point and its exact source level.
#[derive(Clone, Debug)]
pub struct Inequality {
    point: Point,
    level: f64,
}

impl Inequality {
    pub fn new(x: f64, y: f64, z: f64, level: f64) -> Result<Self, ConstraintError> {
        Self::with_c(x, y, z, level, 0.0)
    }

    pub fn with_c(x: f64, y: f64, z: f64, level: f64, c: f64) -> Result<Self, ConstraintError> {
        require_finite(&[level])?;
        Ok(Self {
            point: Point::with_c(x, y, z, c)?,
            level,
        })
    }

    pub const fn point(&self) -> &Point {
        &self.point
    }

    pub const fn level(&self) -> f64 {
        self.level
    }
}

/// A planar orientation with Surfe strike/dip/polarity conventions.
#[derive(Clone, Debug)]
pub struct Planar {
    point: Point,
    dip: f64,
    strike: f64,
    polarity: Polarity,
    normal: [f64; 3],
    normal_bounds: Option<[[f64; 2]; 3]>,
}

impl Planar {
    pub fn from_normal(
        x: f64,
        y: f64,
        z: f64,
        nx: f64,
        ny: f64,
        nz: f64,
    ) -> Result<Self, ConstraintError> {
        Self::from_normal_with_c(x, y, z, nx, ny, nz, 0.0)
    }

    pub fn from_normal_with_c(
        x: f64,
        y: f64,
        z: f64,
        nx: f64,
        ny: f64,
        nz: f64,
        c: f64,
    ) -> Result<Self, ConstraintError> {
        require_finite(&[nx, ny, nz])?;
        let normal = [nx, ny, nz];
        if is_zero_vector(normal) {
            return Err(ConstraintError::ZeroNormal);
        }
        if !(-1.0..=1.0).contains(&nz) {
            return Err(ConstraintError::NormalZOutOfRange);
        }

        let polarity = if nz < 0.0 {
            Polarity::Overturned
        } else {
            Polarity::Upright
        };
        let dip = nz.acos() * RADIANS_TO_DEGREES;
        let mut dip_direction = ny.atan2(nx) * RADIANS_TO_DEGREES;
        if dip_direction < 0.0 {
            dip_direction += 360.0;
        }
        let strike = 360.0 - dip_direction;

        Ok(Self {
            point: Point::with_c(x, y, z, c)?,
            dip,
            strike,
            polarity,
            normal,
            normal_bounds: None,
        })
    }

    pub fn from_strike_dip_polarity(
        x: f64,
        y: f64,
        z: f64,
        strike: f64,
        dip: f64,
        polarity: Polarity,
    ) -> Result<Self, ConstraintError> {
        Self::from_strike_dip_polarity_with_c(x, y, z, strike, dip, polarity, 0.0)
    }

    pub fn from_strike_dip_polarity_with_c(
        x: f64,
        y: f64,
        z: f64,
        strike: f64,
        dip: f64,
        polarity: Polarity,
        c: f64,
    ) -> Result<Self, ConstraintError> {
        require_finite(&[strike, dip])?;
        let point = Point::with_c(x, y, z, c)?;
        let normal = normal_from_strike_dip_polarity(strike, dip, polarity)?;
        Ok(Self {
            point,
            dip,
            strike,
            polarity,
            normal,
            normal_bounds: None,
        })
    }

    pub fn from_azimuth_dip_polarity(
        x: f64,
        y: f64,
        z: f64,
        azimuth: f64,
        dip: f64,
        polarity: Polarity,
    ) -> Result<Self, ConstraintError> {
        Self::from_azimuth_dip_polarity_with_c(x, y, z, azimuth, dip, polarity, 0.0)
    }

    pub fn from_azimuth_dip_polarity_with_c(
        x: f64,
        y: f64,
        z: f64,
        azimuth: f64,
        dip: f64,
        polarity: Polarity,
        c: f64,
    ) -> Result<Self, ConstraintError> {
        require_finite(&[azimuth])?;
        let strike = if azimuth >= 90.0 {
            azimuth - 90.0
        } else {
            azimuth + 270.0
        };
        Self::from_strike_dip_polarity_with_c(x, y, z, strike, dip, polarity, c)
    }

    pub const fn point(&self) -> &Point {
        &self.point
    }

    pub const fn dip(&self) -> f64 {
        self.dip
    }

    pub const fn strike(&self) -> f64 {
        self.strike
    }

    pub const fn polarity(&self) -> Polarity {
        self.polarity
    }

    pub const fn normal(&self) -> [f64; 3] {
        self.normal
    }

    pub const fn nx(&self) -> f64 {
        self.normal[0]
    }

    pub const fn ny(&self) -> f64 {
        self.normal[1]
    }

    pub const fn nz(&self) -> f64 {
        self.normal[2]
    }

    pub const fn normal_bounds(&self) -> Option<[[f64; 2]; 3]> {
        self.normal_bounds
    }

    pub const fn nx_bounds(&self) -> Option<[f64; 2]> {
        match self.normal_bounds {
            Some(bounds) => Some(bounds[0]),
            None => None,
        }
    }

    pub const fn ny_bounds(&self) -> Option<[f64; 2]> {
        match self.normal_bounds {
            Some(bounds) => Some(bounds[1]),
            None => None,
        }
    }

    pub const fn nz_bounds(&self) -> Option<[f64; 2]> {
        match self.normal_bounds {
            Some(bounds) => Some(bounds[2]),
            None => None,
        }
    }

    pub fn dip_vector(&self) -> [f64; 3] {
        let strike = -self.strike * DEGREES_TO_RADIANS;
        let dip = -self.dip * DEGREES_TO_RADIANS;
        let mut vector = [
            strike.cos() * dip.cos(),
            strike.sin() * dip.cos(),
            dip.sin(),
        ];
        let length = vector
            .into_iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        for component in &mut vector {
            *component /= length;
        }
        vector
    }

    pub fn strike_vector(&self) -> [f64; 3] {
        let strike = -self.strike * DEGREES_TO_RADIANS;
        [-strike.sin(), strike.cos(), 0.0]
    }

    pub fn set_normal_bounds(
        &mut self,
        delta_strike: f64,
        delta_dip: f64,
    ) -> Result<(), ConstraintError> {
        require_finite(&[delta_strike, delta_dip])?;

        let theta = self.strike * DEGREES_TO_RADIANS;
        let phi = self.dip * DEGREES_TO_RADIANS;
        let dtheta = delta_strike * DEGREES_TO_RADIANS;
        let dphi = delta_dip * DEGREES_TO_RADIANS;

        let corners = [
            [
                (dtheta + theta).cos() * (dphi + phi).sin(),
                -(dtheta + theta).sin() * (dphi + phi).sin(),
                (dphi + phi).cos(),
            ],
            [
                -(dtheta - theta).cos() * (dphi - phi).sin(),
                -(dtheta - theta).sin() * (dphi - phi).sin(),
                (dphi - phi).cos(),
            ],
            [
                -(dtheta + theta).cos() * (dphi - phi).sin(),
                (dtheta + theta).sin() * (dphi - phi).sin(),
                (dphi - phi).cos(),
            ],
            [
                (dtheta - theta).cos() * (dphi + phi).sin(),
                (dtheta - theta).sin() * (dphi + phi).sin(),
                (dphi + phi).cos(),
            ],
        ];

        let mut bounds = [[0.0; 2]; 3];
        for axis in 0..3 {
            let mut lower = corners[0][axis];
            let mut upper = corners[0][axis];
            for corner in &corners[1..] {
                if corner[axis] < lower {
                    lower = corner[axis];
                }
                if corner[axis] > upper {
                    upper = corner[axis];
                }
            }
            bounds[axis] = [lower, upper];
        }
        self.normal_bounds = Some(bounds);
        Ok(())
    }
}

fn normal_from_strike_dip_polarity(
    strike_degrees: f64,
    dip_degrees: f64,
    polarity: Polarity,
) -> Result<[f64; 3], ConstraintError> {
    let strike = -strike_degrees * DEGREES_TO_RADIANS;
    let dip = -dip_degrees * DEGREES_TO_RADIANS;

    let vx = strike.cos() * dip.cos();
    let vy = strike.sin() * dip.cos();
    let vz = dip.sin();
    let vpx = -vy;
    let vpy = vx;

    let mut normal = [-vz * vpy, vz * vpx, vx * vpy - vy * vpx];
    let length = normal
        .into_iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    if !length.is_finite() || length == 0.0 {
        return Err(ConstraintError::DegenerateOrientation);
    }
    for component in &mut normal {
        *component /= length;
    }
    if (polarity == Polarity::Overturned && normal[2] > 0.0)
        || (polarity == Polarity::Upright && normal[2] < 0.0)
    {
        for component in &mut normal {
            *component = -*component;
        }
    }
    if normal.into_iter().all(f64::is_finite) {
        Ok(normal)
    } else {
        Err(ConstraintError::DegenerateOrientation)
    }
}

/// A tangent direction and its restricted-range angle envelope.
#[derive(Clone, Debug)]
pub struct Tangent {
    point: Point,
    tangent: [f64; 3],
    angle_bounds: Option<[f64; 2]>,
    inner_product_constraint: f64,
}

impl Tangent {
    pub fn new(x: f64, y: f64, z: f64, tx: f64, ty: f64, tz: f64) -> Result<Self, ConstraintError> {
        Self::with_c(x, y, z, tx, ty, tz, 0.0)
    }

    pub fn with_c(
        x: f64,
        y: f64,
        z: f64,
        tx: f64,
        ty: f64,
        tz: f64,
        c: f64,
    ) -> Result<Self, ConstraintError> {
        require_finite(&[tx, ty, tz])?;
        let tangent = [tx, ty, tz];
        if is_zero_vector(tangent) {
            return Err(ConstraintError::ZeroTangent);
        }
        Ok(Self {
            point: Point::with_c(x, y, z, c)?,
            tangent,
            angle_bounds: None,
            inner_product_constraint: 0.0,
        })
    }

    pub const fn point(&self) -> &Point {
        &self.point
    }

    pub const fn vector(&self) -> [f64; 3] {
        self.tangent
    }

    pub const fn tx(&self) -> f64 {
        self.tangent[0]
    }

    pub const fn ty(&self) -> f64 {
        self.tangent[1]
    }

    pub const fn tz(&self) -> f64 {
        self.tangent[2]
    }

    pub const fn angle_bounds(&self) -> Option<[f64; 2]> {
        self.angle_bounds
    }

    pub const fn angle_lower_bound(&self) -> Option<f64> {
        match self.angle_bounds {
            Some(bounds) => Some(bounds[0]),
            None => None,
        }
    }

    pub const fn angle_upper_bound(&self) -> Option<f64> {
        match self.angle_bounds {
            Some(bounds) => Some(bounds[1]),
            None => None,
        }
    }

    pub const fn inner_product_constraint(&self) -> f64 {
        self.inner_product_constraint
    }

    pub fn set_angle_bounds(&mut self, angle: f64) -> Result<(), ConstraintError> {
        require_finite(&[angle])?;
        let value = ((90.0 - angle) * DEGREES_TO_RADIANS).cos() * 2.0;
        self.angle_bounds = Some(if value < 0.0 {
            [value, 0.0]
        } else {
            [0.0, value]
        });
        Ok(())
    }
}

/// The four constraint categories, kept separate as in frozen Surfe.
#[derive(Clone, Debug, Default)]
pub struct Constraints {
    pub inequalities: Vec<Inequality>,
    pub interfaces: Vec<Interface>,
    pub planars: Vec<Planar>,
    pub tangents: Vec<Tangent>,
}
