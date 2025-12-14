//! Physical and biological constants for the Hodgkin-Huxley model.
//!
//! This module contains fundamental constants used in neuronal modeling,
//! including thermodynamic constants and physiological ion concentrations.

/// Physical constants used in the Hodgkin-Huxley model.
#[derive(Debug, Clone, Copy)]
pub struct PhysicalConstants {
    /// Faraday constant (C/mol)
    pub faraday: f64,
    /// Universal gas constant (J/(mol·K))
    pub gas_constant: f64,
    /// Temperature in Kelvin
    pub temperature: f64,
    /// Temperature in Celsius
    pub temperature_celsius: f64,
    /// RT/F at current temperature (mV)
    pub rt_over_f: f64,
}

impl PhysicalConstants {
    /// Create physical constants at a given temperature.
    ///
    /// # Arguments
    ///
    /// * `temp_celsius` - Temperature in degrees Celsius
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::constants::PhysicalConstants;
    ///
    /// let constants = PhysicalConstants::new(6.3);
    /// assert!((constants.temperature - 279.45).abs() < 0.1);
    /// ```
    pub fn new(temp_celsius: f64) -> Self {
        const FARADAY: f64 = 96485.332; // C/mol
        const R: f64 = 8.314; // J/(mol·K)

        let temp_kelvin = temp_celsius + 273.15;
        let rt_over_f = (R * temp_kelvin) / FARADAY * 1000.0; // Convert to mV

        Self {
            faraday: FARADAY,
            gas_constant: R,
            temperature: temp_kelvin,
            temperature_celsius: temp_celsius,
            rt_over_f,
        }
    }

    /// Standard temperature (6.3°C, as used in original HH experiments).
    pub fn standard() -> Self {
        Self::new(6.3)
    }

    /// Physiological temperature (37°C).
    pub fn physiological() -> Self {
        Self::new(37.0)
    }

    /// Calculate Q10 temperature scaling factor.
    ///
    /// # Arguments
    ///
    /// * `q10` - Q10 value (typically 2-3 for ion channels)
    /// * `reference_temp` - Reference temperature in Celsius
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::constants::PhysicalConstants;
    ///
    /// let constants = PhysicalConstants::new(37.0);
    /// let scale = constants.q10_factor(3.0, 6.3);
    /// assert!(scale > 1.0); // Higher temperature increases rate
    /// ```
    pub fn q10_factor(&self, q10: f64, reference_temp: f64) -> f64 {
        q10.powf((self.temperature_celsius - reference_temp) / 10.0)
    }
}

impl Default for PhysicalConstants {
    fn default() -> Self {
        Self::standard()
    }
}

/// Ion concentrations inside and outside the cell.
///
/// All concentrations in mM (millimolar).
#[derive(Debug, Clone, Copy)]
pub struct IonConcentrations {
    /// Intracellular Na+ concentration (mM)
    pub na_in: f64,
    /// Extracellular Na+ concentration (mM)
    pub na_out: f64,
    /// Intracellular K+ concentration (mM)
    pub k_in: f64,
    /// Extracellular K+ concentration (mM)
    pub k_out: f64,
    /// Intracellular Ca2+ concentration (mM)
    pub ca_in: f64,
    /// Extracellular Ca2+ concentration (mM)
    pub ca_out: f64,
    /// Intracellular Cl- concentration (mM)
    pub cl_in: f64,
    /// Extracellular Cl- concentration (mM)
    pub cl_out: f64,
}

impl IonConcentrations {
    /// Create ion concentrations with custom values.
    pub fn new(
        na_in: f64, na_out: f64,
        k_in: f64, k_out: f64,
        ca_in: f64, ca_out: f64,
        cl_in: f64, cl_out: f64,
    ) -> Self {
        Self {
            na_in, na_out,
            k_in, k_out,
            ca_in, ca_out,
            cl_in, cl_out,
        }
    }

    /// Standard squid axon concentrations (from Hodgkin & Huxley 1952).
    pub fn squid_axon() -> Self {
        Self {
            na_in: 50.0,
            na_out: 437.0,
            k_in: 400.0,
            k_out: 20.0,
            ca_in: 0.0001,
            ca_out: 2.0,
            cl_in: 52.0,
            cl_out: 560.0,
        }
    }

    /// Mammalian neuron concentrations (physiological).
    pub fn mammalian() -> Self {
        Self {
            na_in: 15.0,
            na_out: 145.0,
            k_in: 140.0,
            k_out: 5.0,
            ca_in: 0.0001,
            ca_out: 2.0,
            cl_in: 10.0,
            cl_out: 120.0,
        }
    }

    /// Calculate Nernst potential for an ion.
    ///
    /// # Arguments
    ///
    /// * `conc_in` - Intracellular concentration (mM)
    /// * `conc_out` - Extracellular concentration (mM)
    /// * `valence` - Ion valence (e.g., +1 for Na+, +2 for Ca2+, -1 for Cl-)
    /// * `constants` - Physical constants
    ///
    /// # Returns
    ///
    /// Nernst potential in mV
    ///
    /// # Examples
    ///
    /// ```
    /// use hodgkin_huxley::constants::{IonConcentrations, PhysicalConstants};
    ///
    /// let conc = IonConcentrations::squid_axon();
    /// let phys = PhysicalConstants::standard();
    /// let e_na = IonConcentrations::nernst(conc.na_in, conc.na_out, 1.0, &phys);
    /// assert!(e_na > 40.0 && e_na < 60.0); // Typical Na+ reversal potential
    /// ```
    pub fn nernst(conc_in: f64, conc_out: f64, valence: f64, constants: &PhysicalConstants) -> f64 {
        (constants.rt_over_f / valence) * (conc_out / conc_in).ln()
    }

    /// Get sodium Nernst potential.
    pub fn e_na(&self, constants: &PhysicalConstants) -> f64 {
        Self::nernst(self.na_in, self.na_out, 1.0, constants)
    }

    /// Get potassium Nernst potential.
    pub fn e_k(&self, constants: &PhysicalConstants) -> f64 {
        Self::nernst(self.k_in, self.k_out, 1.0, constants)
    }

    /// Get calcium Nernst potential.
    pub fn e_ca(&self, constants: &PhysicalConstants) -> f64 {
        Self::nernst(self.ca_in, self.ca_out, 2.0, constants)
    }

    /// Get chloride Nernst potential.
    pub fn e_cl(&self, constants: &PhysicalConstants) -> f64 {
        Self::nernst(self.cl_in, self.cl_out, -1.0, constants)
    }
}

impl Default for IonConcentrations {
    fn default() -> Self {
        Self::squid_axon()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_physical_constants() {
        let constants = PhysicalConstants::new(6.3);
        assert!((constants.temperature - 279.45).abs() < 0.1);
        // RT/F at 6.3°C should be around 23.26 mV
        // R = 8.314 J/(mol·K), F = 96485.332 C/mol
        // RT/F = (8.314 * 279.45) / 96485.332 * 1000 ≈ 24.04 mV
        assert!((constants.rt_over_f - 24.04).abs() < 0.5);
    }

    #[test]
    fn test_q10_scaling() {
        let constants = PhysicalConstants::new(37.0);
        let scale = constants.q10_factor(3.0, 6.3);
        // Should be approximately 3^3.07 ≈ 27.6
        assert!(scale > 20.0 && scale < 35.0);
    }

    #[test]
    fn test_nernst_potential() {
        let conc = IonConcentrations::squid_axon();
        let phys = PhysicalConstants::standard();

        let e_na = conc.e_na(&phys);
        let e_k = conc.e_k(&phys);

        // Na+ should be positive (higher outside)
        assert!(e_na > 40.0);
        // K+ should be negative (higher inside)
        assert!(e_k < -60.0);
    }

    #[test]
    fn test_mammalian_concentrations() {
        let conc = IonConcentrations::mammalian();
        let phys = PhysicalConstants::physiological();

        let e_na = conc.e_na(&phys);
        let e_k = conc.e_k(&phys);

        // Typical mammalian values
        assert!(e_na > 50.0 && e_na < 70.0);
        assert!(e_k > -100.0 && e_k < -70.0);
    }
}
