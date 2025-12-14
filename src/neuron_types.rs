//! Predefined neuron type configurations.
//!
//! This module provides convenient constructors for different types of neurons
//! with physiologically realistic parameters, including regular spiking (RS),
//! fast spiking (FS), and intrinsically bursting (IB) neurons.

use crate::channels::{CalciumPotassiumChannel, LeakChannel, PotassiumChannel, SodiumChannel};
use crate::constants::{IonConcentrations, PhysicalConstants};

/// Configuration for a specific neuron type.
#[derive(Debug, Clone)]
pub struct NeuronConfig {
    /// Membrane capacitance (µF/cm²)
    pub c_m: f64,
    /// Sodium channel
    pub na_channel: SodiumChannel,
    /// Potassium delayed rectifier channel
    pub k_channel: PotassiumChannel,
    /// Calcium-activated potassium channel
    pub k_ca_channel: CalciumPotassiumChannel,
    /// Leak channel
    pub leak_channel: LeakChannel,
    /// Ion concentrations
    pub concentrations: IonConcentrations,
    /// Physical constants
    pub constants: PhysicalConstants,
    /// Neuron type description
    pub description: String,
}

impl NeuronConfig {
    /// Create a custom neuron configuration.
    pub fn new(
        c_m: f64,
        na_channel: SodiumChannel,
        k_channel: PotassiumChannel,
        k_ca_channel: CalciumPotassiumChannel,
        leak_channel: LeakChannel,
        concentrations: IonConcentrations,
        constants: PhysicalConstants,
        description: String,
    ) -> Self {
        Self {
            c_m,
            na_channel,
            k_channel,
            k_ca_channel,
            leak_channel,
            concentrations,
            constants,
            description,
        }
    }

    /// Standard Hodgkin-Huxley squid axon model.
    ///
    /// This is the classical model from Hodgkin & Huxley (1952) with
    /// parameters fitted to squid giant axon at 6.3°C.
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::neuron_types::NeuronConfig;
    ///
    /// let config = NeuronConfig::squid_axon();
    /// assert_eq!(config.c_m, 1.0);
    /// ```
    pub fn squid_axon() -> Self {
        Self {
            c_m: 1.0, // µF/cm²
            na_channel: SodiumChannel::new(120.0),
            k_channel: PotassiumChannel::new(36.0),
            k_ca_channel: CalciumPotassiumChannel::new(0.0), // Not in original HH model
            leak_channel: LeakChannel::new(0.3, -54.387),
            concentrations: IonConcentrations::squid_axon(),
            constants: PhysicalConstants::standard(),
            description: "Hodgkin-Huxley squid giant axon (1952)".to_string(),
        }
    }

    /// Regular spiking (RS) cortical pyramidal neuron.
    ///
    /// These neurons show adaptation and relatively slow, regular firing.
    /// Common in cortical layer 2/3 and layer 5.
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::neuron_types::NeuronConfig;
    ///
    /// let config = NeuronConfig::regular_spiking();
    /// assert!(config.k_ca_channel.g_max > 0.0); // Has adaptation
    /// ```
    pub fn regular_spiking() -> Self {
        Self {
            c_m: 1.0,
            na_channel: SodiumChannel::new(100.0),
            k_channel: PotassiumChannel::new(30.0),
            k_ca_channel: CalciumPotassiumChannel::new(0.3), // Moderate adaptation
            leak_channel: LeakChannel::new(0.1, -70.0),
            concentrations: IonConcentrations::mammalian(),
            constants: PhysicalConstants::physiological(),
            description: "Regular spiking cortical pyramidal neuron".to_string(),
        }
    }

    /// Fast spiking (FS) interneuron.
    ///
    /// These neurons fire rapidly with minimal adaptation, characteristic
    /// of parvalbumin-positive interneurons.
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::neuron_types::NeuronConfig;
    ///
    /// let config = NeuronConfig::fast_spiking();
    /// assert!(config.na_channel.g_max > 120.0); // High Na+ for fast spikes
    /// ```
    pub fn fast_spiking() -> Self {
        Self {
            c_m: 0.75, // Lower capacitance for faster response
            na_channel: SodiumChannel::new(150.0), // Higher Na+ conductance
            k_channel: PotassiumChannel::new(50.0), // Higher K+ conductance
            k_ca_channel: CalciumPotassiumChannel::new(0.05), // Minimal adaptation
            leak_channel: LeakChannel::new(0.15, -65.0),
            concentrations: IonConcentrations::mammalian(),
            constants: PhysicalConstants::physiological(),
            description: "Fast spiking interneuron".to_string(),
        }
    }

    /// Intrinsically bursting (IB) neuron.
    ///
    /// These neurons produce bursts of action potentials, common in
    /// layer 5 pyramidal neurons.
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::neuron_types::NeuronConfig;
    ///
    /// let config = NeuronConfig::intrinsically_bursting();
    /// assert!(config.k_ca_channel.g_max > 0.5); // Strong adaptation for bursting
    /// ```
    pub fn intrinsically_bursting() -> Self {
        Self {
            c_m: 1.2,
            na_channel: SodiumChannel::new(120.0),
            k_channel: PotassiumChannel::new(25.0), // Lower K+ for bursting
            k_ca_channel: CalciumPotassiumChannel::new(0.8), // Strong slow AHP
            leak_channel: LeakChannel::new(0.08, -75.0),
            concentrations: IonConcentrations::mammalian(),
            constants: PhysicalConstants::physiological(),
            description: "Intrinsically bursting layer 5 pyramidal neuron".to_string(),
        }
    }

    /// Low-threshold spiking (LTS) interneuron.
    ///
    /// These neurons can fire from hyperpolarized states and show
    /// characteristic rebound bursts.
    pub fn low_threshold_spiking() -> Self {
        Self {
            c_m: 1.0,
            na_channel: SodiumChannel::new(90.0),
            k_channel: PotassiumChannel::new(28.0),
            k_ca_channel: CalciumPotassiumChannel::new(0.2),
            leak_channel: LeakChannel::new(0.12, -72.0),
            concentrations: IonConcentrations::mammalian(),
            constants: PhysicalConstants::physiological(),
            description: "Low-threshold spiking interneuron".to_string(),
        }
    }

    /// Chattering neuron.
    ///
    /// Produces high-frequency bursts followed by silence, found in
    /// superficial cortical layers.
    pub fn chattering() -> Self {
        Self {
            c_m: 0.9,
            na_channel: SodiumChannel::new(140.0), // High Na+ for rapid firing
            k_channel: PotassiumChannel::new(40.0),
            k_ca_channel: CalciumPotassiumChannel::new(1.2), // Very strong adaptation
            leak_channel: LeakChannel::new(0.1, -68.0),
            concentrations: IonConcentrations::mammalian(),
            constants: PhysicalConstants::physiological(),
            description: "Chattering neuron".to_string(),
        }
    }

    /// Calculate equilibrium potentials for this configuration.
    pub fn reversal_potentials(&self) -> ReversalPotentials {
        ReversalPotentials {
            e_na: self.concentrations.e_na(&self.constants),
            e_k: self.concentrations.e_k(&self.constants),
            e_ca: self.concentrations.e_ca(&self.constants),
            e_cl: self.concentrations.e_cl(&self.constants),
        }
    }

    /// Get resting potential (approximate).
    ///
    /// This is computed by finding the voltage where all currents balance
    /// at resting gating variable values.
    pub fn resting_potential(&self) -> f64 {
        // Use leak reversal as initial guess
        let mut v = self.leak_channel.e_leak;

        // Newton-Raphson iteration to find I_total = 0
        for _ in 0..50 {
            let m_inf = SodiumChannel::m_inf(v);
            let h_inf = SodiumChannel::h_inf(v);
            let n_inf = PotassiumChannel::n_inf(v);
            let a_inf = CalciumPotassiumChannel::a_inf(v);
            let b_inf = CalciumPotassiumChannel::b_inf(v);

            let rev = self.reversal_potentials();

            let i_na = self.na_channel.current(v, m_inf, h_inf, rev.e_na);
            let i_k = self.k_channel.current(v, n_inf, rev.e_k);
            let i_k_ca = self.k_ca_channel.current(v, a_inf, b_inf, rev.e_k);
            let i_leak = self.leak_channel.current(v);

            let i_total = i_na + i_k + i_k_ca + i_leak;

            // Check convergence
            if i_total.abs() < 1e-6 {
                break;
            }

            // Numerical derivative
            let dv = 0.1;
            let m_inf_2 = SodiumChannel::m_inf(v + dv);
            let h_inf_2 = SodiumChannel::h_inf(v + dv);
            let n_inf_2 = PotassiumChannel::n_inf(v + dv);
            let a_inf_2 = CalciumPotassiumChannel::a_inf(v + dv);
            let b_inf_2 = CalciumPotassiumChannel::b_inf(v + dv);

            let i_na_2 = self.na_channel.current(v + dv, m_inf_2, h_inf_2, rev.e_na);
            let i_k_2 = self.k_channel.current(v + dv, n_inf_2, rev.e_k);
            let i_k_ca_2 = self.k_ca_channel.current(v + dv, a_inf_2, b_inf_2, rev.e_k);
            let i_leak_2 = self.leak_channel.current(v + dv);

            let i_total_2 = i_na_2 + i_k_2 + i_k_ca_2 + i_leak_2;

            let di_dv = (i_total_2 - i_total) / dv;

            if di_dv.abs() > 1e-10 {
                v -= i_total / di_dv;
            } else {
                break;
            }
        }

        v
    }
}

/// Reversal potentials for different ions.
#[derive(Debug, Clone, Copy)]
pub struct ReversalPotentials {
    /// Sodium reversal potential (mV)
    pub e_na: f64,
    /// Potassium reversal potential (mV)
    pub e_k: f64,
    /// Calcium reversal potential (mV)
    pub e_ca: f64,
    /// Chloride reversal potential (mV)
    pub e_cl: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_squid_axon_config() {
        let config = NeuronConfig::squid_axon();
        assert_eq!(config.c_m, 1.0);
        assert_eq!(config.na_channel.g_max, 120.0);
        assert_eq!(config.k_channel.g_max, 36.0);
    }

    #[test]
    fn test_reversal_potentials() {
        let config = NeuronConfig::squid_axon();
        let rev = config.reversal_potentials();

        // Sodium should be positive
        assert!(rev.e_na > 0.0);
        // Potassium should be negative
        assert!(rev.e_k < 0.0);
    }

    #[test]
    fn test_resting_potential() {
        let config = NeuronConfig::squid_axon();
        let v_rest = config.resting_potential();

        // Should be somewhere between -80 and -50 mV
        assert!(v_rest > -80.0 && v_rest < -50.0);
    }

    #[test]
    fn test_regular_spiking() {
        let config = NeuronConfig::regular_spiking();
        assert!(config.k_ca_channel.g_max > 0.0);
        assert_eq!(config.description, "Regular spiking cortical pyramidal neuron");
    }

    #[test]
    fn test_fast_spiking() {
        let config = NeuronConfig::fast_spiking();
        assert!(config.na_channel.g_max > 120.0);
        assert!(config.c_m < 1.0); // Lower capacitance
    }

    #[test]
    fn test_intrinsically_bursting() {
        let config = NeuronConfig::intrinsically_bursting();
        assert!(config.k_ca_channel.g_max > 0.5);
    }

    #[test]
    fn test_all_neuron_types() {
        // Ensure all neuron types can be created without panic
        let configs = vec![
            NeuronConfig::squid_axon(),
            NeuronConfig::regular_spiking(),
            NeuronConfig::fast_spiking(),
            NeuronConfig::intrinsically_bursting(),
            NeuronConfig::low_threshold_spiking(),
            NeuronConfig::chattering(),
        ];

        for config in configs {
            let v_rest = config.resting_potential();
            assert!(v_rest.is_finite());
            assert!(v_rest > -100.0 && v_rest < 0.0);
        }
    }
}
