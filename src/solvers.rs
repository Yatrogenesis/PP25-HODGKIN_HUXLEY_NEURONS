//! Numerical integration methods for solving differential equations.
//!
//! This module provides various integration schemes for solving the
//! Hodgkin-Huxley equations, including Runge-Kutta 4th order (RK4)
//! and exponential Euler methods optimized for stiff systems.

use crate::error::{HHError, Result};
use nalgebra::DVector;

/// Function type for the system of ODEs: dy/dt = f(t, y)
pub type ODEFunction = dyn Fn(f64, &DVector<f64>) -> DVector<f64>;

/// Fourth-order Runge-Kutta integrator.
///
/// The RK4 method is a popular explicit method that provides good accuracy
/// for smooth systems. For the Hodgkin-Huxley equations, it balances
/// accuracy and computational efficiency.
#[derive(Debug, Clone)]
pub struct RK4Integrator {
    /// Integration time step (ms)
    pub dt: f64,
    /// Maximum allowed value for state variables
    pub max_value: f64,
}

impl RK4Integrator {
    /// Create a new RK4 integrator.
    ///
    /// # Arguments
    ///
    /// * `dt` - Integration time step in ms (typically 0.01-0.1 ms)
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::solvers::RK4Integrator;
    ///
    /// let integrator = RK4Integrator::new(0.01).unwrap();
    /// assert_eq!(integrator.dt, 0.01);
    /// ```
    pub fn new(dt: f64) -> Result<Self> {
        if dt <= 0.0 || !dt.is_finite() {
            return Err(HHError::InvalidParameter {
                parameter: "dt".to_string(),
                value: dt,
                reason: "Time step must be positive and finite".to_string(),
            });
        }

        if dt > 1.0 {
            return Err(HHError::InvalidParameter {
                parameter: "dt".to_string(),
                value: dt,
                reason: "Time step too large for accurate integration".to_string(),
            });
        }

        Ok(Self {
            dt,
            max_value: 1e6,
        })
    }

    /// Perform one integration step using RK4 method.
    ///
    /// # Arguments
    ///
    /// * `t` - Current time
    /// * `y` - Current state vector
    /// * `f` - Function computing dy/dt
    ///
    /// # Returns
    ///
    /// New state vector at time t + dt
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::solvers::RK4Integrator;
    /// use nalgebra::DVector;
    ///
    /// let integrator = RK4Integrator::new(0.01).unwrap();
    /// let y = DVector::from_vec(vec![1.0, 0.0]);
    ///
    /// // Simple harmonic oscillator: dy/dt = [y[1], -y[0]]
    /// let f = |_t: f64, y: &DVector<f64>| -> DVector<f64> {
    ///     DVector::from_vec(vec![y[1], -y[0]])
    /// };
    ///
    /// let y_next = integrator.step(0.0, &y, &f).unwrap();
    /// assert!(y_next.len() == 2);
    /// ```
    pub fn step(&self, t: f64, y: &DVector<f64>, f: &ODEFunction) -> Result<DVector<f64>> {
        // RK4 formula:
        // k1 = f(t, y)
        // k2 = f(t + dt/2, y + dt*k1/2)
        // k3 = f(t + dt/2, y + dt*k2/2)
        // k4 = f(t + dt, y + dt*k3)
        // y_next = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

        let k1 = f(t, y);
        self.check_finite(&k1, "k1")?;

        let y2 = y + &k1 * (self.dt / 2.0);
        let k2 = f(t + self.dt / 2.0, &y2);
        self.check_finite(&k2, "k2")?;

        let y3 = y + &k2 * (self.dt / 2.0);
        let k3 = f(t + self.dt / 2.0, &y3);
        self.check_finite(&k3, "k3")?;

        let y4 = y + &k3 * self.dt;
        let k4 = f(t + self.dt, &y4);
        self.check_finite(&k4, "k4")?;

        let y_next = y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (self.dt / 6.0);
        self.check_finite(&y_next, "y_next")?;
        self.check_bounds(&y_next, t + self.dt)?;

        Ok(y_next)
    }

    /// Check if all values in vector are finite.
    fn check_finite(&self, v: &DVector<f64>, location: &str) -> Result<()> {
        for (i, &val) in v.iter().enumerate() {
            if !val.is_finite() {
                return Err(HHError::NonFiniteValue {
                    location: format!("{}[{}]", location, i),
                    value: val,
                });
            }
        }
        Ok(())
    }

    /// Check if values are within reasonable bounds.
    fn check_bounds(&self, v: &DVector<f64>, t: f64) -> Result<()> {
        for (i, &val) in v.iter().enumerate() {
            if val.abs() > self.max_value {
                return Err(HHError::IntegrationError {
                    time: t,
                    reason: format!("Value at index {} exceeded bounds: {}", i, val),
                });
            }
        }
        Ok(())
    }

    /// Integrate from t0 to t_final with given initial conditions.
    ///
    /// # Arguments
    ///
    /// * `t0` - Initial time
    /// * `t_final` - Final time
    /// * `y0` - Initial state vector
    /// * `f` - Function computing dy/dt
    ///
    /// # Returns
    ///
    /// Vector of (time, state) tuples at each integration step
    pub fn integrate(
        &self,
        t0: f64,
        t_final: f64,
        y0: &DVector<f64>,
        f: &ODEFunction,
    ) -> Result<Vec<(f64, DVector<f64>)>> {
        let n_steps = ((t_final - t0) / self.dt).ceil() as usize;
        let mut results = Vec::with_capacity(n_steps + 1);

        let mut t = t0;
        let mut y = y0.clone();

        results.push((t, y.clone()));

        while t < t_final {
            let dt_actual = self.dt.min(t_final - t);
            let integrator = if (dt_actual - self.dt).abs() < 1e-10 {
                self.clone()
            } else {
                Self { dt: dt_actual, max_value: self.max_value }
            };

            y = integrator.step(t, &y, f)?;
            t += dt_actual;
            results.push((t, y.clone()));
        }

        Ok(results)
    }
}

/// Exponential Euler integrator for stiff systems.
///
/// This method is particularly useful for equations of the form:
/// dy/dt = (y_inf - y) / tau
///
/// which appears in the gating variable equations of the Hodgkin-Huxley model.
#[derive(Debug, Clone)]
pub struct ExponentialEulerIntegrator {
    /// Integration time step (ms)
    pub dt: f64,
}

impl ExponentialEulerIntegrator {
    /// Create a new exponential Euler integrator.
    pub fn new(dt: f64) -> Result<Self> {
        if dt <= 0.0 || !dt.is_finite() {
            return Err(HHError::InvalidParameter {
                parameter: "dt".to_string(),
                value: dt,
                reason: "Time step must be positive and finite".to_string(),
            });
        }

        Ok(Self { dt })
    }

    /// Integrate a single gating variable using exponential Euler.
    ///
    /// For dy/dt = (y_inf - y) / tau, the exact solution over time step dt is:
    /// y(t + dt) = y_inf + (y(t) - y_inf) * exp(-dt / tau)
    ///
    /// # Arguments
    ///
    /// * `y` - Current value of gating variable
    /// * `y_inf` - Steady-state value
    /// * `tau` - Time constant (ms)
    ///
    /// # Returns
    ///
    /// New value of gating variable
    pub fn step_gating(&self, y: f64, y_inf: f64, tau: f64) -> f64 {
        y_inf + (y - y_inf) * (-self.dt / tau).exp()
    }

    /// Alternative formulation using alpha and beta rate constants.
    ///
    /// For dy/dt = alpha * (1 - y) - beta * y
    pub fn step_alpha_beta(&self, y: f64, alpha: f64, beta: f64) -> f64 {
        let y_inf = alpha / (alpha + beta);
        let tau = 1.0 / (alpha + beta);
        self.step_gating(y, y_inf, tau)
    }
}

/// Adaptive step size controller using embedded Runge-Kutta methods.
///
/// This provides error estimation and automatic step size adjustment.
#[derive(Debug, Clone)]
pub struct AdaptiveIntegrator {
    /// Initial time step (ms)
    pub dt_init: f64,
    /// Minimum allowed time step (ms)
    pub dt_min: f64,
    /// Maximum allowed time step (ms)
    pub dt_max: f64,
    /// Error tolerance
    pub tol: f64,
}

impl AdaptiveIntegrator {
    /// Create a new adaptive integrator.
    pub fn new(dt_init: f64, dt_min: f64, dt_max: f64, tol: f64) -> Result<Self> {
        if dt_min <= 0.0 || dt_init < dt_min || dt_max < dt_init {
            return Err(HHError::InvalidParameter {
                parameter: "dt".to_string(),
                value: dt_init,
                reason: "Invalid time step bounds".to_string(),
            });
        }

        Ok(Self {
            dt_init,
            dt_min,
            dt_max,
            tol,
        })
    }

    /// Standard adaptive integrator for HH equations.
    pub fn standard() -> Result<Self> {
        Self::new(0.01, 0.001, 0.1, 1e-6)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_rk4_creation() {
        let integrator = RK4Integrator::new(0.01).unwrap();
        assert_eq!(integrator.dt, 0.01);

        // Invalid time steps
        assert!(RK4Integrator::new(0.0).is_err());
        assert!(RK4Integrator::new(-0.01).is_err());
        assert!(RK4Integrator::new(2.0).is_err());
    }

    #[test]
    fn test_rk4_exponential_decay() {
        // Test with dy/dt = -y, exact solution y(t) = y0 * exp(-t)
        let integrator = RK4Integrator::new(0.01).unwrap();
        let y0 = DVector::from_vec(vec![1.0]);

        let f = |_t: f64, y: &DVector<f64>| -> DVector<f64> {
            DVector::from_vec(vec![-y[0]])
        };

        let results = integrator.integrate(0.0, 1.0, &y0, &f).unwrap();
        let y_final = &results.last().unwrap().1;

        // Compare with exact solution
        let exact = (-1.0_f64).exp();
        assert_relative_eq!(y_final[0], exact, epsilon = 1e-4);
    }

    #[test]
    fn test_rk4_harmonic_oscillator() {
        // Test with simple harmonic oscillator
        // dy1/dt = y2, dy2/dt = -y1
        // Solution: y1 = cos(t), y2 = -sin(t)
        let integrator = RK4Integrator::new(0.01).unwrap();
        let y0 = DVector::from_vec(vec![1.0, 0.0]);

        let f = |_t: f64, y: &DVector<f64>| -> DVector<f64> {
            DVector::from_vec(vec![y[1], -y[0]])
        };

        let results = integrator.integrate(0.0, 2.0 * std::f64::consts::PI, &y0, &f).unwrap();
        let y_final = &results.last().unwrap().1;

        // After one period, should return to initial conditions
        assert_relative_eq!(y_final[0], 1.0, epsilon = 1e-3);
        assert_relative_eq!(y_final[1], 0.0, epsilon = 1e-3);
    }

    #[test]
    fn test_exponential_euler_gating() {
        let integrator = ExponentialEulerIntegrator::new(0.01).unwrap();

        // Test convergence to steady state
        let mut y = 0.0;
        let y_inf = 0.8;
        let tau = 5.0;

        // Simulate for 50 ms (10 time constants)
        for _ in 0..5000 {
            y = integrator.step_gating(y, y_inf, tau);
        }

        // Should be very close to steady state
        assert_relative_eq!(y, y_inf, epsilon = 1e-3);
    }

    #[test]
    fn test_exponential_euler_alpha_beta() {
        let integrator = ExponentialEulerIntegrator::new(0.01).unwrap();

        let alpha = 0.1;
        let beta = 0.05;
        let expected_inf = alpha / (alpha + beta);

        let mut y = 0.0;
        for _ in 0..10000 {
            y = integrator.step_alpha_beta(y, alpha, beta);
        }

        assert_relative_eq!(y, expected_inf, epsilon = 1e-3);
    }

    #[test]
    fn test_adaptive_integrator_creation() {
        let integrator = AdaptiveIntegrator::standard().unwrap();
        assert!(integrator.dt_init > integrator.dt_min);
        assert!(integrator.dt_max > integrator.dt_init);

        // Invalid bounds
        assert!(AdaptiveIntegrator::new(0.1, 0.2, 0.3, 1e-6).is_err());
    }
}
