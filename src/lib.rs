//! # Hodgkin-Huxley Neuron Library
//!
//! A production-ready implementation of the Hodgkin-Huxley neuron model with exact
//! biophysical equations from the seminal 1952 paper.
//!
//! ## Overview
//!
//! This library provides a complete implementation of the Hodgkin-Huxley equations
//! for modeling action potentials in neurons. It includes:
//!
//! - **Accurate biophysics**: Exact equations from Hodgkin & Huxley (1952)
//! - **Multiple neuron types**: Regular spiking, fast spiking, intrinsically bursting, etc.
//! - **Numerical solvers**: RK4 and exponential Euler integrators
//! - **Temperature effects**: Q10 scaling for realistic temperature dependence
//! - **Ion channels**: Na⁺, K⁺, Ca²⁺-activated K⁺, and leak channels
//!
//! ## Mathematical Model
//!
//! The membrane voltage is governed by:
//!
//! ```text
//! C_m * dV/dt = -I_Na - I_K - I_K(Ca) - I_leak + I_ext
//! ```
//!
//! Where:
//! - `I_Na = g_Na * m³ * h * (V - E_Na)` - Fast sodium current
//! - `I_K = g_K * n⁴ * (V - E_K)` - Delayed rectifier potassium current
//! - `I_K(Ca) = g_K(Ca) * a * b * (V - E_K)` - Calcium-activated potassium current
//! - `I_leak = g_leak * (V - E_leak)` - Leak current
//!
//! Gating variables evolve according to:
//!
//! ```text
//! dm/dt = α_m(V) * (1 - m) - β_m(V) * m
//! dh/dt = α_h(V) * (1 - h) - β_h(V) * h
//! dn/dt = α_n(V) * (1 - n) - β_n(V) * n
//! da/dt = α_a(V) * (1 - a) - β_a(V) * a
//! db/dt = α_b(V) * (1 - b) - β_b(V) * b
//! ```
//!
//! ## Quick Start
//!
//! ```
//! use hodgkin_huxley::{HodgkinHuxleyNeuron, neuron_types::NeuronConfig};
//! use nalgebra::DVector;
//!
//! // Create a regular spiking cortical neuron
//! let config = NeuronConfig::regular_spiking();
//! let mut neuron = HodgkinHuxleyNeuron::new(config).unwrap();
//!
//! // Initialize at resting state
//! neuron.initialize_rest();
//!
//! // Apply a current pulse and simulate
//! let i_ext = 10.0; // µA/cm²
//! let dt = 0.01;    // ms
//! let duration = 50.0; // ms
//!
//! for _ in 0..(duration / dt) as usize {
//!     neuron.step(dt, i_ext).unwrap();
//! }
//!
//! // Check if neuron spiked
//! let spikes = neuron.detect_spikes(-20.0);
//! println!("Number of spikes: {}", spikes.len());
//! ```
//!
//! ## Examples
//!
//! ### Simulate different neuron types
//!
//! ```
//! use hodgkin_huxley::{HodgkinHuxleyNeuron, neuron_types::NeuronConfig};
//!
//! // Fast spiking interneuron
//! let fs_neuron = HodgkinHuxleyNeuron::new(NeuronConfig::fast_spiking()).unwrap();
//!
//! // Intrinsically bursting neuron
//! let ib_neuron = HodgkinHuxleyNeuron::new(NeuronConfig::intrinsically_bursting()).unwrap();
//!
//! // Classical squid axon
//! let squid_neuron = HodgkinHuxleyNeuron::new(NeuronConfig::squid_axon()).unwrap();
//! ```
//!
//! ### Record voltage trace
//!
//! ```
//! use hodgkin_huxley::{HodgkinHuxleyNeuron, neuron_types::NeuronConfig};
//!
//! let mut neuron = HodgkinHuxleyNeuron::new(NeuronConfig::regular_spiking()).unwrap();
//! neuron.initialize_rest();
//!
//! let mut voltage_trace = Vec::new();
//! let dt = 0.01;
//!
//! for t in 0..5000 {
//!     let current = if t > 1000 && t < 4000 { 10.0 } else { 0.0 };
//!     neuron.step(dt, current).unwrap();
//!     voltage_trace.push((t as f64 * dt, neuron.voltage()));
//! }
//! ```

pub mod channels;
pub mod constants;
pub mod error;
pub mod neuron_types;
pub mod solvers;

use channels::{CalciumPotassiumChannel, LeakChannel, PotassiumChannel, SodiumChannel};
use constants::{IonConcentrations, PhysicalConstants};
use error::{HHError, Result};
use nalgebra::DVector;
use neuron_types::NeuronConfig;
use solvers::RK4Integrator;

/// Complete Hodgkin-Huxley neuron model with 6 state variables.
///
/// State vector: [V, m, h, n, a, b]
/// - V: Membrane potential (mV)
/// - m: Sodium activation gate
/// - h: Sodium inactivation gate
/// - n: Potassium activation gate
/// - a: Calcium-activated potassium activation gate
/// - b: Calcium-activated potassium inactivation gate
#[derive(Debug, Clone)]
pub struct HodgkinHuxleyNeuron {
    /// Current state vector
    state: DVector<f64>,
    /// Membrane capacitance (µF/cm²)
    c_m: f64,
    /// Sodium channel
    na_channel: SodiumChannel,
    /// Potassium channel
    k_channel: PotassiumChannel,
    /// Calcium-activated potassium channel
    k_ca_channel: CalciumPotassiumChannel,
    /// Leak channel
    leak_channel: LeakChannel,
    /// Ion concentrations
    concentrations: IonConcentrations,
    /// Physical constants
    constants: PhysicalConstants,
    /// Reversal potentials (cached)
    e_na: f64,
    e_k: f64,
    /// Current simulation time (ms)
    time: f64,
    /// Voltage history for spike detection
    voltage_history: Vec<(f64, f64)>,
    /// Integration method
    integrator: RK4Integrator,
}

impl HodgkinHuxleyNeuron {
    /// Create a new Hodgkin-Huxley neuron from configuration.
    ///
    /// # Arguments
    ///
    /// * `config` - Neuron configuration
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::{HodgkinHuxleyNeuron, neuron_types::NeuronConfig};
    ///
    /// let neuron = HodgkinHuxleyNeuron::new(NeuronConfig::squid_axon()).unwrap();
    /// ```
    pub fn new(config: NeuronConfig) -> Result<Self> {
        if config.c_m <= 0.0 {
            return Err(HHError::InvalidParameter {
                parameter: "c_m".to_string(),
                value: config.c_m,
                reason: "Capacitance must be positive".to_string(),
            });
        }

        let e_na = config.concentrations.e_na(&config.constants);
        let e_k = config.concentrations.e_k(&config.constants);

        // Initialize state at resting values
        let v_init = -65.0;
        let state = DVector::from_vec(vec![
            v_init,
            SodiumChannel::m_inf(v_init),
            SodiumChannel::h_inf(v_init),
            PotassiumChannel::n_inf(v_init),
            CalciumPotassiumChannel::a_inf(v_init),
            CalciumPotassiumChannel::b_inf(v_init),
        ]);

        Ok(Self {
            state,
            c_m: config.c_m,
            na_channel: config.na_channel,
            k_channel: config.k_channel,
            k_ca_channel: config.k_ca_channel,
            leak_channel: config.leak_channel,
            concentrations: config.concentrations,
            constants: config.constants,
            e_na,
            e_k,
            time: 0.0,
            voltage_history: Vec::new(),
            integrator: RK4Integrator::new(0.01)?,
        })
    }

    /// Create a standard squid axon neuron.
    pub fn squid_axon() -> Result<Self> {
        Self::new(NeuronConfig::squid_axon())
    }

    /// Create a regular spiking neuron.
    pub fn regular_spiking() -> Result<Self> {
        Self::new(NeuronConfig::regular_spiking())
    }

    /// Create a fast spiking neuron.
    pub fn fast_spiking() -> Result<Self> {
        Self::new(NeuronConfig::fast_spiking())
    }

    /// Create an intrinsically bursting neuron.
    pub fn intrinsically_bursting() -> Result<Self> {
        Self::new(NeuronConfig::intrinsically_bursting())
    }

    /// Initialize neuron at resting state.
    ///
    /// Sets all gating variables to their steady-state values at the
    /// resting potential.
    pub fn initialize_rest(&mut self) {
        let v_rest = self.leak_channel.e_leak;
        self.state[0] = v_rest;
        self.state[1] = SodiumChannel::m_inf(v_rest);
        self.state[2] = SodiumChannel::h_inf(v_rest);
        self.state[3] = PotassiumChannel::n_inf(v_rest);
        self.state[4] = CalciumPotassiumChannel::a_inf(v_rest);
        self.state[5] = CalciumPotassiumChannel::b_inf(v_rest);
        self.time = 0.0;
        self.voltage_history.clear();
    }

    /// Set membrane voltage.
    pub fn set_voltage(&mut self, v: f64) -> Result<()> {
        if !v.is_finite() {
            return Err(HHError::NonFiniteValue {
                location: "voltage".to_string(),
                value: v,
            });
        }
        self.state[0] = v;
        Ok(())
    }

    /// Set gating variables.
    pub fn set_gates(&mut self, m: f64, h: f64, n: f64, a: f64, b: f64) -> Result<()> {
        let gates = [m, h, n, a, b];
        let names = ["m", "h", "n", "a", "b"];

        for (i, (&gate, &name)) in gates.iter().zip(names.iter()).enumerate() {
            if !gate.is_finite() || gate < 0.0 || gate > 1.0 {
                return Err(HHError::InvalidParameter {
                    parameter: name.to_string(),
                    value: gate,
                    reason: "Gating variable must be in [0, 1]".to_string(),
                });
            }
            self.state[i + 1] = gate;
        }
        Ok(())
    }

    /// Get current membrane voltage (mV).
    pub fn voltage(&self) -> f64 {
        self.state[0]
    }

    /// Get gating variables (m, h, n, a, b).
    pub fn gates(&self) -> (f64, f64, f64, f64, f64) {
        (self.state[1], self.state[2], self.state[3], self.state[4], self.state[5])
    }

    /// Get current simulation time (ms).
    pub fn time(&self) -> f64 {
        self.time
    }

    /// Get complete state vector.
    pub fn state(&self) -> &DVector<f64> {
        &self.state
    }

    /// Set integrator time step.
    pub fn set_timestep(&mut self, dt: f64) -> Result<()> {
        self.integrator = RK4Integrator::new(dt)?;
        Ok(())
    }

    /// Compute all ionic currents at current state.
    ///
    /// # Returns
    ///
    /// Tuple of (I_Na, I_K, I_K_Ca, I_leak) in µA/cm²
    pub fn currents(&self) -> (f64, f64, f64, f64) {
        let v = self.state[0];
        let (m, h, n, a, b) = self.gates();

        let i_na = self.na_channel.current(v, m, h, self.e_na);
        let i_k = self.k_channel.current(v, n, self.e_k);
        let i_k_ca = self.k_ca_channel.current(v, a, b, self.e_k);
        let i_leak = self.leak_channel.current(v);

        (i_na, i_k, i_k_ca, i_leak)
    }

    /// Compute total membrane current.
    ///
    /// # Arguments
    ///
    /// * `i_ext` - External applied current (µA/cm²)
    ///
    /// # Returns
    ///
    /// Total current in µA/cm²
    pub fn total_current(&self, i_ext: f64) -> f64 {
        let (i_na, i_k, i_k_ca, i_leak) = self.currents();
        i_na + i_k + i_k_ca + i_leak - i_ext
    }

    /// Compute all conductances at current state.
    ///
    /// # Returns
    ///
    /// Tuple of (g_Na, g_K, g_K_Ca, g_leak) in mS/cm²
    pub fn conductances(&self) -> (f64, f64, f64, f64) {
        let (m, h, n, a, b) = self.gates();

        let g_na = self.na_channel.conductance(m, h);
        let g_k = self.k_channel.conductance(n);
        let g_k_ca = self.k_ca_channel.conductance(a, b);
        let g_leak = self.leak_channel.g_leak;

        (g_na, g_k, g_k_ca, g_leak)
    }

    /// Compute derivatives of state variables (dY/dt).
    ///
    /// # Arguments
    ///
    /// * `i_ext` - External applied current (µA/cm²)
    ///
    /// # Returns
    ///
    /// Vector of derivatives [dV/dt, dm/dt, dh/dt, dn/dt, da/dt, db/dt]
    fn derivatives(&self, i_ext: f64) -> DVector<f64> {
        let v = self.state[0];
        let (m, h, n, a, b) = self.gates();

        // Ionic currents
        let (i_na, i_k, i_k_ca, i_leak) = self.currents();

        // dV/dt = (I_ext - I_Na - I_K - I_K_Ca - I_leak) / C_m
        let dv_dt = (i_ext - i_na - i_k - i_k_ca - i_leak) / self.c_m;

        // Gating variable rate constants
        let alpha_m = SodiumChannel::alpha_m(v);
        let beta_m = SodiumChannel::beta_m(v);
        let alpha_h = SodiumChannel::alpha_h(v);
        let beta_h = SodiumChannel::beta_h(v);
        let alpha_n = PotassiumChannel::alpha_n(v);
        let beta_n = PotassiumChannel::beta_n(v);
        let alpha_a = CalciumPotassiumChannel::alpha_a(v);
        let beta_a = CalciumPotassiumChannel::beta_a(v);
        let alpha_b = CalciumPotassiumChannel::alpha_b(v);
        let beta_b = CalciumPotassiumChannel::beta_b(v);

        // Gating variable derivatives
        let dm_dt = alpha_m * (1.0 - m) - beta_m * m;
        let dh_dt = alpha_h * (1.0 - h) - beta_h * h;
        let dn_dt = alpha_n * (1.0 - n) - beta_n * n;
        let da_dt = alpha_a * (1.0 - a) - beta_a * a;
        let db_dt = alpha_b * (1.0 - b) - beta_b * b;

        DVector::from_vec(vec![dv_dt, dm_dt, dh_dt, dn_dt, da_dt, db_dt])
    }

    /// Advance simulation by one time step using RK4 integration.
    ///
    /// # Arguments
    ///
    /// * `dt` - Time step (ms)
    /// * `i_ext` - External applied current (µA/cm²)
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::{HodgkinHuxleyNeuron, neuron_types::NeuronConfig};
    ///
    /// let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
    /// neuron.initialize_rest();
    ///
    /// // Apply 10 µA/cm² current for 1 ms
    /// neuron.step(0.01, 10.0).unwrap();
    /// ```
    pub fn step(&mut self, dt: f64, i_ext: f64) -> Result<()> {
        // Update integrator if needed
        if (dt - self.integrator.dt).abs() > 1e-10 {
            self.integrator = RK4Integrator::new(dt)?;
        }

        // Save neuron parameters for closure
        let c_m = self.c_m;
        let na_channel = self.na_channel;
        let k_channel = self.k_channel;
        let k_ca_channel = self.k_ca_channel;
        let leak_channel = self.leak_channel;
        let concentrations = self.concentrations;
        let constants = self.constants;
        let e_na = self.e_na;
        let e_k = self.e_k;

        // Define ODE function with captured parameters
        let f = move |_t: f64, y: &DVector<f64>| -> DVector<f64> {
            let temp_neuron = HodgkinHuxleyNeuron {
                state: y.clone(),
                c_m,
                na_channel,
                k_channel,
                k_ca_channel,
                leak_channel,
                concentrations,
                constants,
                e_na,
                e_k,
                time: 0.0,
                voltage_history: Vec::new(),
                integrator: RK4Integrator::new(0.01).unwrap(),
            };
            temp_neuron.derivatives(i_ext)
        };

        // Integrate
        let current_time = self.time;
        let current_state = self.state.clone();
        self.state = self.integrator.step(current_time, &current_state, &f)?;
        self.time += dt;

        // Record voltage for spike detection
        self.voltage_history.push((self.time, self.voltage()));

        // Keep only recent history (last 100 ms)
        if let Some(&(first_time, _)) = self.voltage_history.first() {
            if self.time - first_time > 100.0 {
                let cutoff_time = self.time - 100.0;
                self.voltage_history.retain(|(t, _)| *t >= cutoff_time);
            }
        }

        Ok(())
    }

    /// Simulate for a duration with constant external current.
    ///
    /// # Arguments
    ///
    /// * `duration` - Simulation duration (ms)
    /// * `dt` - Time step (ms)
    /// * `i_ext` - External applied current (µA/cm²)
    ///
    /// # Returns
    ///
    /// Vector of (time, voltage) pairs
    pub fn simulate(&mut self, duration: f64, dt: f64, i_ext: f64) -> Result<Vec<(f64, f64)>> {
        let n_steps = (duration / dt).ceil() as usize;
        let mut trace = Vec::with_capacity(n_steps);

        for _ in 0..n_steps {
            self.step(dt, i_ext)?;
            trace.push((self.time, self.voltage()));
        }

        Ok(trace)
    }

    /// Detect action potential spikes in voltage history.
    ///
    /// # Arguments
    ///
    /// * `threshold` - Voltage threshold for spike detection (mV), typically -20 mV
    ///
    /// # Returns
    ///
    /// Vector of spike times (ms)
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::{HodgkinHuxleyNeuron, neuron_types::NeuronConfig};
    ///
    /// let mut neuron = HodgkinHuxleyNeuron::regular_spiking().unwrap();
    /// neuron.initialize_rest();
    /// neuron.simulate(100.0, 0.01, 10.0).unwrap();
    ///
    /// let spikes = neuron.detect_spikes(-20.0);
    /// println!("Detected {} spikes", spikes.len());
    /// ```
    pub fn detect_spikes(&self, threshold: f64) -> Vec<f64> {
        let mut spikes = Vec::new();
        let mut above_threshold = false;

        for window in self.voltage_history.windows(2) {
            let (t1, v1) = window[0];
            let (_, v2) = window[1];

            // Detect upward threshold crossing
            if !above_threshold && v1 <= threshold && v2 > threshold {
                spikes.push(t1);
                above_threshold = true;
            }

            // Reset when falling below threshold
            if above_threshold && v2 <= threshold {
                above_threshold = false;
            }
        }

        spikes
    }

    /// Calculate instantaneous firing rate from spike times.
    ///
    /// # Arguments
    ///
    /// * `spikes` - Vector of spike times (ms)
    ///
    /// # Returns
    ///
    /// Average firing rate in Hz
    pub fn firing_rate(spikes: &[f64]) -> f64 {
        if spikes.len() < 2 {
            return 0.0;
        }

        let duration = spikes[spikes.len() - 1] - spikes[0];
        if duration <= 0.0 {
            return 0.0;
        }

        ((spikes.len() - 1) as f64 / duration) * 1000.0 // Convert ms to s
    }

    /// Calculate interspike intervals.
    ///
    /// # Arguments
    ///
    /// * `spikes` - Vector of spike times (ms)
    ///
    /// # Returns
    ///
    /// Vector of interspike intervals (ms)
    pub fn interspike_intervals(spikes: &[f64]) -> Vec<f64> {
        if spikes.len() < 2 {
            return Vec::new();
        }

        spikes.windows(2)
            .map(|w| w[1] - w[0])
            .collect()
    }

    /// Clear voltage history.
    pub fn clear_history(&mut self) {
        self.voltage_history.clear();
    }

    /// Get voltage history.
    pub fn voltage_history(&self) -> &[(f64, f64)] {
        &self.voltage_history
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_neuron_creation() {
        let neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        assert_eq!(neuron.c_m, 1.0);
        assert_eq!(neuron.state.len(), 6);
    }

    #[test]
    fn test_initialize_rest() {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.initialize_rest();

        let v = neuron.voltage();
        let (m, h, n, _, _) = neuron.gates();

        // At rest, voltage should be near leak reversal
        assert!((v - neuron.leak_channel.e_leak).abs() < 1.0);

        // Gating variables should be at steady state
        assert_relative_eq!(m, SodiumChannel::m_inf(v), epsilon = 1e-6);
        assert_relative_eq!(h, SodiumChannel::h_inf(v), epsilon = 1e-6);
        assert_relative_eq!(n, PotassiumChannel::n_inf(v), epsilon = 1e-6);
    }

    #[test]
    fn test_set_voltage() {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.set_voltage(-70.0).unwrap();
        assert_eq!(neuron.voltage(), -70.0);

        // Should reject non-finite values
        assert!(neuron.set_voltage(f64::NAN).is_err());
        assert!(neuron.set_voltage(f64::INFINITY).is_err());
    }

    #[test]
    fn test_set_gates() {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.set_gates(0.1, 0.9, 0.3, 0.5, 0.7).unwrap();
        let (m, h, n, a, b) = neuron.gates();
        assert_eq!(m, 0.1);
        assert_eq!(h, 0.9);
        assert_eq!(n, 0.3);
        assert_eq!(a, 0.5);
        assert_eq!(b, 0.7);

        // Should reject out-of-range values
        assert!(neuron.set_gates(-0.1, 0.5, 0.5, 0.5, 0.5).is_err());
        assert!(neuron.set_gates(1.1, 0.5, 0.5, 0.5, 0.5).is_err());
    }

    #[test]
    fn test_currents_at_rest() {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.initialize_rest();

        let (i_na, i_k, i_k_ca, i_leak) = neuron.currents();

        // At rest with gates at steady state, individual currents may be significant
        // but they balance out. The leak reversal potential determines the exact state.
        // With E_leak = -54.387 mV, there will be steady-state currents.
        assert!(i_na.abs() < 50.0, "i_na = {} exceeded bounds", i_na);
        assert!(i_k.abs() < 50.0, "i_k = {} exceeded bounds", i_k);
        assert!(i_k_ca.abs() < 1.0, "i_k_ca = {} exceeded bounds", i_k_ca);
        assert!(i_leak.abs() < 50.0, "i_leak = {} exceeded bounds", i_leak);

        // All values should be finite
        assert!(i_na.is_finite());
        assert!(i_k.is_finite());
        assert!(i_k_ca.is_finite());
        assert!(i_leak.is_finite());
    }

    #[test]
    fn test_step_integration() {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.initialize_rest();

        let v_before = neuron.voltage();
        let t_before = neuron.time();

        // Apply small current
        neuron.step(0.01, 1.0).unwrap();

        let v_after = neuron.voltage();
        let t_after = neuron.time();

        // Voltage should change slightly
        assert_ne!(v_before, v_after);
        // Time should advance
        assert_eq!(t_after, t_before + 0.01);
    }

    #[test]
    fn test_action_potential() {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.initialize_rest();

        // Apply strong current pulse
        let dt = 0.01;
        let current = 20.0; // µA/cm²

        for _ in 0..1000 {
            neuron.step(dt, current).unwrap();
        }

        // Should generate at least one spike
        let spikes = neuron.detect_spikes(-20.0);
        assert!(spikes.len() > 0, "No spikes detected");
    }

    #[test]
    fn test_spike_detection() {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.initialize_rest();

        // Generate spikes with current injection
        neuron.simulate(50.0, 0.01, 15.0).unwrap();

        let spikes = neuron.detect_spikes(-20.0);
        assert!(spikes.len() > 0);

        // Spikes should be in chronological order
        for window in spikes.windows(2) {
            assert!(window[1] > window[0]);
        }
    }

    #[test]
    fn test_firing_rate() {
        let spikes = vec![10.0, 20.0, 30.0, 40.0]; // 4 spikes over 30 ms
        let rate = HodgkinHuxleyNeuron::firing_rate(&spikes);

        // 3 intervals over 30 ms = 100 Hz
        assert_relative_eq!(rate, 100.0, epsilon = 0.1);
    }

    #[test]
    fn test_interspike_intervals() {
        let spikes = vec![10.0, 20.0, 35.0];
        let isis = HodgkinHuxleyNeuron::interspike_intervals(&spikes);

        assert_eq!(isis.len(), 2);
        assert_eq!(isis[0], 10.0);
        assert_eq!(isis[1], 15.0);
    }

    #[test]
    fn test_different_neuron_types() {
        let neurons = vec![
            HodgkinHuxleyNeuron::squid_axon().unwrap(),
            HodgkinHuxleyNeuron::regular_spiking().unwrap(),
            HodgkinHuxleyNeuron::fast_spiking().unwrap(),
            HodgkinHuxleyNeuron::intrinsically_bursting().unwrap(),
        ];

        for mut neuron in neurons {
            neuron.initialize_rest();
            // Should be able to simulate without errors
            assert!(neuron.simulate(10.0, 0.01, 10.0).is_ok());
        }
    }

    #[test]
    fn test_conductances() {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.initialize_rest();

        let (g_na, g_k, g_k_ca, g_leak) = neuron.conductances();

        // All conductances should be positive
        assert!(g_na >= 0.0);
        assert!(g_k >= 0.0);
        assert!(g_k_ca >= 0.0);
        assert!(g_leak > 0.0);

        // At rest, g_na and g_k should be small
        assert!(g_na < 10.0);
        assert!(g_k < 5.0);
    }

    #[test]
    fn test_voltage_history() {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.initialize_rest();

        neuron.simulate(5.0, 0.01, 0.0).unwrap();
        let history = neuron.voltage_history();

        // Should have approximately 500 samples
        assert!(history.len() > 400 && history.len() < 600);

        // Times should be monotonically increasing
        for window in history.windows(2) {
            assert!(window[1].0 > window[0].0);
        }
    }

    #[test]
    fn test_clear_history() {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.initialize_rest();
        neuron.simulate(5.0, 0.01, 0.0).unwrap();

        assert!(neuron.voltage_history().len() > 0);

        neuron.clear_history();
        assert_eq!(neuron.voltage_history().len(), 0);
    }
}
