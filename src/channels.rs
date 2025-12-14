//! Ion channel implementations for the Hodgkin-Huxley model.
//!
//! This module implements the voltage-gated ion channels described in the
//! original Hodgkin & Huxley (1952) paper, including sodium, potassium,
//! and leak channels, along with additional calcium-activated potassium channels.

use crate::constants::PhysicalConstants;

/// Sodium channel gating variables (m and h).
///
/// The sodium current is modeled as I_Na = g_Na * m³ * h * (V - E_Na)
#[derive(Debug, Clone, Copy)]
pub struct SodiumChannel {
    /// Maximum conductance (mS/cm²)
    pub g_max: f64,
}

impl SodiumChannel {
    /// Create a new sodium channel with specified conductance.
    pub fn new(g_max: f64) -> Self {
        Self { g_max }
    }

    /// Standard HH sodium channel (120 mS/cm²).
    pub fn standard() -> Self {
        Self::new(120.0)
    }

    /// Alpha function for m gate (activation).
    ///
    /// From Hodgkin & Huxley 1952, equation 12.
    /// Units: ms⁻¹
    pub fn alpha_m(v: f64) -> f64 {
        let v_shifted = v + 40.0; // Shift from resting potential
        if v_shifted.abs() < 1e-4 {
            // Avoid division by zero using L'Hôpital's rule
            1.0
        } else {
            0.1 * v_shifted / (1.0 - (-v_shifted / 10.0).exp())
        }
    }

    /// Beta function for m gate (activation).
    ///
    /// From Hodgkin & Huxley 1952, equation 13.
    /// Units: ms⁻¹
    pub fn beta_m(v: f64) -> f64 {
        4.0 * (-(v + 65.0) / 18.0).exp()
    }

    /// Alpha function for h gate (inactivation).
    ///
    /// From Hodgkin & Huxley 1952, equation 14.
    /// Units: ms⁻¹
    pub fn alpha_h(v: f64) -> f64 {
        0.07 * (-(v + 65.0) / 20.0).exp()
    }

    /// Beta function for h gate (inactivation).
    ///
    /// From Hodgkin & Huxley 1952, equation 15.
    /// Units: ms⁻¹
    pub fn beta_h(v: f64) -> f64 {
        1.0 / (1.0 + (-(v + 35.0) / 10.0).exp())
    }

    /// Steady-state value for m gate.
    pub fn m_inf(v: f64) -> f64 {
        let alpha = Self::alpha_m(v);
        let beta = Self::beta_m(v);
        alpha / (alpha + beta)
    }

    /// Time constant for m gate (ms).
    pub fn tau_m(v: f64) -> f64 {
        let alpha = Self::alpha_m(v);
        let beta = Self::beta_m(v);
        1.0 / (alpha + beta)
    }

    /// Steady-state value for h gate.
    pub fn h_inf(v: f64) -> f64 {
        let alpha = Self::alpha_h(v);
        let beta = Self::beta_h(v);
        alpha / (alpha + beta)
    }

    /// Time constant for h gate (ms).
    pub fn tau_h(v: f64) -> f64 {
        let alpha = Self::alpha_h(v);
        let beta = Self::beta_h(v);
        1.0 / (alpha + beta)
    }

    /// Calculate sodium current.
    ///
    /// # Arguments
    ///
    /// * `v` - Membrane voltage (mV)
    /// * `m` - m gate variable
    /// * `h` - h gate variable
    /// * `e_na` - Sodium reversal potential (mV)
    ///
    /// # Returns
    ///
    /// Current in µA/cm²
    pub fn current(&self, v: f64, m: f64, h: f64, e_na: f64) -> f64 {
        self.g_max * m.powi(3) * h * (v - e_na)
    }

    /// Calculate conductance.
    pub fn conductance(&self, m: f64, h: f64) -> f64 {
        self.g_max * m.powi(3) * h
    }
}

/// Potassium delayed rectifier channel (n gate).
///
/// The potassium current is modeled as I_K = g_K * n⁴ * (V - E_K)
#[derive(Debug, Clone, Copy)]
pub struct PotassiumChannel {
    /// Maximum conductance (mS/cm²)
    pub g_max: f64,
}

impl PotassiumChannel {
    /// Create a new potassium channel with specified conductance.
    pub fn new(g_max: f64) -> Self {
        Self { g_max }
    }

    /// Standard HH potassium channel (36 mS/cm²).
    pub fn standard() -> Self {
        Self::new(36.0)
    }

    /// Alpha function for n gate.
    ///
    /// From Hodgkin & Huxley 1952, equation 16.
    /// Units: ms⁻¹
    pub fn alpha_n(v: f64) -> f64 {
        let v_shifted = v + 55.0;
        if v_shifted.abs() < 1e-4 {
            // Avoid division by zero
            0.1
        } else {
            0.01 * v_shifted / (1.0 - (-v_shifted / 10.0).exp())
        }
    }

    /// Beta function for n gate.
    ///
    /// From Hodgkin & Huxley 1952, equation 17.
    /// Units: ms⁻¹
    pub fn beta_n(v: f64) -> f64 {
        0.125 * (-(v + 65.0) / 80.0).exp()
    }

    /// Steady-state value for n gate.
    pub fn n_inf(v: f64) -> f64 {
        let alpha = Self::alpha_n(v);
        let beta = Self::beta_n(v);
        alpha / (alpha + beta)
    }

    /// Time constant for n gate (ms).
    pub fn tau_n(v: f64) -> f64 {
        let alpha = Self::alpha_n(v);
        let beta = Self::beta_n(v);
        1.0 / (alpha + beta)
    }

    /// Calculate potassium current.
    ///
    /// # Arguments
    ///
    /// * `v` - Membrane voltage (mV)
    /// * `n` - n gate variable
    /// * `e_k` - Potassium reversal potential (mV)
    ///
    /// # Returns
    ///
    /// Current in µA/cm²
    pub fn current(&self, v: f64, n: f64, e_k: f64) -> f64 {
        self.g_max * n.powi(4) * (v - e_k)
    }

    /// Calculate conductance.
    pub fn conductance(&self, n: f64) -> f64 {
        self.g_max * n.powi(4)
    }
}

/// Calcium-activated potassium channel (a and b gates).
///
/// This slow afterhyperpolarization (AHP) channel contributes to spike
/// frequency adaptation and burst firing patterns.
#[derive(Debug, Clone, Copy)]
pub struct CalciumPotassiumChannel {
    /// Maximum conductance (mS/cm²)
    pub g_max: f64,
}

impl CalciumPotassiumChannel {
    /// Create a new calcium-activated potassium channel.
    pub fn new(g_max: f64) -> Self {
        Self { g_max }
    }

    /// Standard calcium-activated K+ channel (0.3 mS/cm²).
    pub fn standard() -> Self {
        Self::new(0.3)
    }

    /// Alpha function for a gate (activation).
    ///
    /// Based on models from Connor & Stevens (1971).
    /// Units: ms⁻¹
    pub fn alpha_a(v: f64) -> f64 {
        let v_shifted = v + 50.0;
        if v_shifted.abs() < 1e-4 {
            0.02
        } else {
            0.02 * v_shifted / (1.0 - (-v_shifted / 11.0).exp())
        }
    }

    /// Beta function for a gate (activation).
    pub fn beta_a(v: f64) -> f64 {
        0.175 * (-(v + 75.0) / 11.0).exp()
    }

    /// Alpha function for b gate (inactivation).
    pub fn alpha_b(v: f64) -> f64 {
        0.016 * (-(v + 60.0) / 8.0).exp()
    }

    /// Beta function for b gate (inactivation).
    pub fn beta_b(v: f64) -> f64 {
        0.5 / (1.0 + (-(v + 29.0) / 5.0).exp())
    }

    /// Steady-state value for a gate.
    pub fn a_inf(v: f64) -> f64 {
        let alpha = Self::alpha_a(v);
        let beta = Self::beta_a(v);
        alpha / (alpha + beta)
    }

    /// Time constant for a gate (ms).
    pub fn tau_a(v: f64) -> f64 {
        let alpha = Self::alpha_a(v);
        let beta = Self::beta_a(v);
        1.0 / (alpha + beta)
    }

    /// Steady-state value for b gate.
    pub fn b_inf(v: f64) -> f64 {
        let alpha = Self::alpha_b(v);
        let beta = Self::beta_b(v);
        alpha / (alpha + beta)
    }

    /// Time constant for b gate (ms).
    pub fn tau_b(v: f64) -> f64 {
        let alpha = Self::alpha_b(v);
        let beta = Self::beta_b(v);
        1.0 / (alpha + beta)
    }

    /// Calculate calcium-activated potassium current.
    ///
    /// # Arguments
    ///
    /// * `v` - Membrane voltage (mV)
    /// * `a` - a gate variable
    /// * `b` - b gate variable
    /// * `e_k` - Potassium reversal potential (mV)
    pub fn current(&self, v: f64, a: f64, b: f64, e_k: f64) -> f64 {
        self.g_max * a * b * (v - e_k)
    }

    /// Calculate conductance.
    pub fn conductance(&self, a: f64, b: f64) -> f64 {
        self.g_max * a * b
    }
}

/// Leak channel (passive membrane conductance).
///
/// Represents the combined effect of all non-voltage-gated channels.
#[derive(Debug, Clone, Copy)]
pub struct LeakChannel {
    /// Leak conductance (mS/cm²)
    pub g_leak: f64,
    /// Leak reversal potential (mV)
    pub e_leak: f64,
}

impl LeakChannel {
    /// Create a new leak channel.
    pub fn new(g_leak: f64, e_leak: f64) -> Self {
        Self { g_leak, e_leak }
    }

    /// Standard HH leak channel (0.3 mS/cm², -54.387 mV).
    pub fn standard() -> Self {
        Self::new(0.3, -54.387)
    }

    /// Calculate leak current.
    ///
    /// # Arguments
    ///
    /// * `v` - Membrane voltage (mV)
    ///
    /// # Returns
    ///
    /// Current in µA/cm²
    pub fn current(&self, v: f64) -> f64 {
        self.g_leak * (v - self.e_leak)
    }
}

/// Apply Q10 temperature correction to rate constants.
///
/// # Arguments
///
/// * `rate` - Rate constant at reference temperature
/// * `constants` - Physical constants (includes current temperature)
/// * `q10` - Q10 value for the process
/// * `reference_temp` - Reference temperature in Celsius
///
/// # Returns
///
/// Temperature-corrected rate constant
pub fn apply_q10(rate: f64, constants: &PhysicalConstants, q10: f64, reference_temp: f64) -> f64 {
    rate * constants.q10_factor(q10, reference_temp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sodium_channel() {
        let na = SodiumChannel::standard();
        assert_eq!(na.g_max, 120.0);

        // Test at resting potential
        let v = -65.0;
        let m = SodiumChannel::m_inf(v);
        let h = SodiumChannel::h_inf(v);

        // At rest, m should be small (closed) and h should be large (not inactivated)
        assert!(m < 0.1);
        assert!(h > 0.5);
    }

    #[test]
    fn test_potassium_channel() {
        let k = PotassiumChannel::standard();
        assert_eq!(k.g_max, 36.0);

        // Test at resting potential
        let v = -65.0;
        let n = PotassiumChannel::n_inf(v);

        // At rest, n should be small (closed)
        assert!(n < 0.4);
    }

    #[test]
    fn test_gating_variable_bounds() {
        // Test that gating variables stay in [0, 1] range
        for v in (-100..=50).step_by(10) {
            let v = v as f64;

            let m = SodiumChannel::m_inf(v);
            let h = SodiumChannel::h_inf(v);
            let n = PotassiumChannel::n_inf(v);

            assert!(m >= 0.0 && m <= 1.0, "m out of bounds at V={}", v);
            assert!(h >= 0.0 && h <= 1.0, "h out of bounds at V={}", v);
            assert!(n >= 0.0 && n <= 1.0, "n out of bounds at V={}", v);
        }
    }

    #[test]
    fn test_time_constants() {
        // Time constants should be positive
        for v in (-100..=50).step_by(10) {
            let v = v as f64;

            let tau_m = SodiumChannel::tau_m(v);
            let tau_h = SodiumChannel::tau_h(v);
            let tau_n = PotassiumChannel::tau_n(v);

            assert!(tau_m > 0.0, "tau_m negative at V={}", v);
            assert!(tau_h > 0.0, "tau_h negative at V={}", v);
            assert!(tau_n > 0.0, "tau_n negative at V={}", v);
        }
    }

    #[test]
    fn test_leak_channel() {
        let leak = LeakChannel::standard();
        let current = leak.current(-65.0);

        // At resting potential near E_leak, current should be small
        assert!(current.abs() < 5.0);
    }
}
