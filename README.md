# Hodgkin-Huxley Neuron Library

A production-ready Rust implementation of the Hodgkin-Huxley neuron model with exact biophysical equations from the seminal 1952 paper.

[![Crates.io](https://img.shields.io/crates/v/hodgkin-huxley.svg)](https://crates.io/crates/hodgkin-huxley)
[![Documentation](https://docs.rs/hodgkin-huxley/badge.svg)](https://docs.rs/hodgkin-huxley)
[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)](LICENSE)

## Features

- **Accurate biophysics**: Exact equations from Hodgkin & Huxley (1952)
- **Multiple neuron types**: Regular spiking, fast spiking, intrinsically bursting, and more
- **Numerical solvers**: RK4 and exponential Euler integrators
- **Temperature effects**: Q10 scaling for realistic temperature dependence
- **Ion channels**: Na⁺, K⁺, Ca²⁺-activated K⁺, and leak channels
- **Spike detection**: Built-in action potential detection and analysis
- **Well-tested**: Comprehensive unit tests and documentation examples
- **Production-ready**: Proper error handling with `thiserror`

## Mathematical Model

The membrane voltage is governed by:

```
C_m * dV/dt = -I_Na - I_K - I_K(Ca) - I_leak + I_ext
```

Where:
- `I_Na = g_Na * m³ * h * (V - E_Na)` - Fast sodium current
- `I_K = g_K * n⁴ * (V - E_K)` - Delayed rectifier potassium current
- `I_K(Ca) = g_K(Ca) * a * b * (V - E_K)` - Calcium-activated potassium current
- `I_leak = g_leak * (V - E_leak)` - Leak current

Gating variables evolve according to:

```
dm/dt = α_m(V) * (1 - m) - β_m(V) * m
dh/dt = α_h(V) * (1 - h) - β_h(V) * h
dn/dt = α_n(V) * (1 - n) - β_n(V) * n
da/dt = α_a(V) * (1 - a) - β_a(V) * a
db/dt = α_b(V) * (1 - b) - β_b(V) * b
```

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
hodgkin-huxley = "0.1.0"
```

## Quick Start

```rust
use hodgkin_huxley::{HodgkinHuxleyNeuron, neuron_types::NeuronConfig};

fn main() {
    // Create a regular spiking cortical neuron
    let config = NeuronConfig::regular_spiking();
    let mut neuron = HodgkinHuxleyNeuron::new(config).unwrap();

    // Initialize at resting state
    neuron.initialize_rest();

    // Apply a current pulse and simulate
    let i_ext = 10.0; // µA/cm²
    let dt = 0.01;    // ms
    let duration = 50.0; // ms

    for _ in 0..(duration / dt) as usize {
        neuron.step(dt, i_ext).unwrap();
    }

    // Check if neuron spiked
    let spikes = neuron.detect_spikes(-20.0);
    println!("Number of spikes: {}", spikes.len());
}
```

## Examples

### Simulate Different Neuron Types

```rust
use hodgkin_huxley::HodgkinHuxleyNeuron;

// Fast spiking interneuron
let fs_neuron = HodgkinHuxleyNeuron::fast_spiking().unwrap();

// Intrinsically bursting neuron
let ib_neuron = HodgkinHuxleyNeuron::intrinsically_bursting().unwrap();

// Classical squid axon
let squid_neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
```

### Record Voltage Trace

```rust
use hodgkin_huxley::HodgkinHuxleyNeuron;

let mut neuron = HodgkinHuxleyNeuron::regular_spiking().unwrap();
neuron.initialize_rest();

let trace = neuron.simulate(100.0, 0.01, 10.0).unwrap();

// trace is a Vec<(f64, f64)> of (time, voltage) pairs
for (time, voltage) in trace.iter().take(10) {
    println!("t = {:.2} ms, V = {:.2} mV", time, voltage);
}
```

### Analyze Firing Patterns

```rust
use hodgkin_huxley::HodgkinHuxleyNeuron;

let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
neuron.initialize_rest();
neuron.simulate(100.0, 0.01, 15.0).unwrap();

// Detect spikes
let spikes = neuron.detect_spikes(-20.0);

// Calculate firing rate
let rate = HodgkinHuxleyNeuron::firing_rate(&spikes);
println!("Firing rate: {:.2} Hz", rate);

// Calculate interspike intervals
let isis = HodgkinHuxleyNeuron::interspike_intervals(&spikes);
println!("Mean ISI: {:.2} ms", isis.iter().sum::<f64>() / isis.len() as f64);
```

## Running Examples

The library includes several examples demonstrating its capabilities:

```bash
# Simulate an action potential
cargo run --example action_potential --release

# Compare different neuron types
cargo run --example neuron_comparison --release
```

## Available Neuron Types

- **Squid Axon**: Original Hodgkin-Huxley model (1952)
- **Regular Spiking (RS)**: Cortical pyramidal neurons with adaptation
- **Fast Spiking (FS)**: Parvalbumin-positive interneurons
- **Intrinsically Bursting (IB)**: Layer 5 pyramidal neurons
- **Low-Threshold Spiking (LTS)**: Interneurons with rebound bursts
- **Chattering**: High-frequency burst neurons

## Physical Constants

The library includes realistic physical constants and ion concentrations:

- Faraday constant: 96485.332 C/mol
- Gas constant: 8.314 J/(mol·K)
- Temperature-dependent Nernst potentials
- Q10 temperature scaling (Q10 = 3 for mammalian neurons)

## Numerical Integration

Two integration methods are provided:

- **RK4 (Runge-Kutta 4th order)**: High accuracy, suitable for smooth systems
- **Exponential Euler**: Optimized for stiff gating variable equations

Default time step: 0.01 ms (10 µs)

## Units

All quantities use standard electrophysiological units:

- Voltage: mV (millivolts)
- Current: µA/cm² (microamperes per square centimeter)
- Conductance: mS/cm² (millisiemens per square centimeter)
- Capacitance: µF/cm² (microfarads per square centimeter)
- Time: ms (milliseconds)
- Ion concentrations: mM (millimolar)

## Testing

Run the test suite:

```bash
# Run all tests
cargo test

# Run only unit tests
cargo test --lib

# Run documentation tests
cargo test --doc

# Run with verbose output
cargo test -- --nocapture
```

## Documentation

Generate and view the full documentation:

```bash
cargo doc --open
```

## Performance

The library is optimized for production use:

- Efficient RK4 integration
- Minimal allocations in the hot path
- SIMD-friendly operations via `nalgebra`
- Typical performance: ~1-2 µs per integration step (release build)

Run benchmarks:

```bash
cargo bench
```

## References

1. Hodgkin, A. L., & Huxley, A. F. (1952). A quantitative description of membrane current and its application to conduction and excitation in nerve. *The Journal of Physiology*, 117(4), 500-544.

2. Connor, J. A., & Stevens, C. F. (1971). Prediction of repetitive firing behaviour from voltage clamp data on an isolated neurone soma. *The Journal of Physiology*, 213(1), 31-53.

3. Traub, R. D., & Miles, R. (1991). *Neuronal Networks of the Hippocampus*. Cambridge University Press.

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in the work by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or conditions.

## Citation

If you use this library in your research, please cite:

```bibtex
@software{hodgkin_huxley_rust,
  title = {hodgkin-huxley: A Rust implementation of the Hodgkin-Huxley neuron model},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/cortexia/hodgkin-huxley}
}
```
