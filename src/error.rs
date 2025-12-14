//! Error types for the Hodgkin-Huxley neuron library.

use thiserror::Error;

/// Errors that can occur during neuron simulation.
#[derive(Error, Debug, Clone, PartialEq)]
pub enum HHError {
    /// Invalid parameter value
    #[error("Invalid parameter: {parameter} with value {value}: {reason}")]
    InvalidParameter {
        parameter: String,
        value: f64,
        reason: String,
    },

    /// Numerical integration error
    #[error("Integration error at t={time}: {reason}")]
    IntegrationError { time: f64, reason: String },

    /// State vector dimension mismatch
    #[error("State dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: usize, actual: usize },

    /// Non-finite value detected
    #[error("Non-finite value detected in {location}: {value}")]
    NonFiniteValue { location: String, value: f64 },

    /// Temperature out of valid range
    #[error("Temperature {temp_celsius}°C is out of valid range (0-50°C)")]
    InvalidTemperature { temp_celsius: f64 },

    /// Concentration out of valid range
    #[error("Ion concentration out of valid range: {ion}")]
    InvalidConcentration { ion: String },
}

/// Result type for Hodgkin-Huxley operations.
pub type Result<T> = std::result::Result<T, HHError>;
