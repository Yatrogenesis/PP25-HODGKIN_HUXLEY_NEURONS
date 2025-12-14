//! Example: Simulating an action potential with the Hodgkin-Huxley model.
//!
//! This example demonstrates how to:
//! 1. Create a squid axon neuron
//! 2. Apply a current pulse
//! 3. Record the voltage trace
//! 4. Detect action potential spikes

use hodgkin_huxley::HodgkinHuxleyNeuron;

fn main() {
    println!("Hodgkin-Huxley Action Potential Simulation");
    println!("===========================================\n");

    // Create a squid axon neuron (original HH model)
    let mut neuron = HodgkinHuxleyNeuron::squid_axon()
        .expect("Failed to create neuron");

    // Initialize at resting state
    neuron.initialize_rest();
    println!("Initial resting voltage: {:.2} mV", neuron.voltage());

    // Simulation parameters
    let dt = 0.01;           // Time step: 10 µs
    let pulse_start = 5.0;   // Start pulse at 5 ms
    let pulse_end = 10.0;    // End pulse at 10 ms
    let pulse_current = 15.0; // Current amplitude: 15 µA/cm²
    let total_time = 50.0;   // Total simulation time: 50 ms

    println!("Applying {:.1} µA/cm² current pulse from {:.1} to {:.1} ms\n",
             pulse_current, pulse_start, pulse_end);

    // Simulate
    let mut voltage_trace = Vec::new();
    let mut time = 0.0;

    while time < total_time {
        // Determine current
        let current = if time >= pulse_start && time < pulse_end {
            pulse_current
        } else {
            0.0
        };

        // Step the simulation
        neuron.step(dt, current)
            .expect("Integration failed");

        voltage_trace.push((time, neuron.voltage()));
        time += dt;
    }

    // Detect spikes
    let spikes = neuron.detect_spikes(-20.0);
    println!("Action potentials detected: {}", spikes.len());

    if !spikes.is_empty() {
        println!("Spike times (ms):");
        for (i, &spike_time) in spikes.iter().enumerate() {
            println!("  Spike {}: {:.2} ms", i + 1, spike_time);
        }

        // Calculate firing rate
        let firing_rate = HodgkinHuxleyNeuron::firing_rate(&spikes);
        println!("\nAverage firing rate: {:.2} Hz", firing_rate);

        // Calculate interspike intervals
        let isis = HodgkinHuxleyNeuron::interspike_intervals(&spikes);
        if !isis.is_empty() {
            println!("Interspike intervals (ms):");
            for (i, &isi) in isis.iter().enumerate() {
                println!("  ISI {}: {:.2} ms", i + 1, isi);
            }
        }
    }

    // Print voltage trace statistics
    let voltages: Vec<f64> = voltage_trace.iter().map(|(_, v)| *v).collect();
    let v_min = voltages.iter().cloned().fold(f64::INFINITY, f64::min);
    let v_max = voltages.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    println!("\nVoltage statistics:");
    println!("  Minimum: {:.2} mV", v_min);
    println!("  Maximum: {:.2} mV", v_max);
    println!("  Range: {:.2} mV", v_max - v_min);

    // Get final gating variables
    let (m, h, n, a, b) = neuron.gates();
    println!("\nFinal gating variables:");
    println!("  m (Na+ activation): {:.4}", m);
    println!("  h (Na+ inactivation): {:.4}", h);
    println!("  n (K+ activation): {:.4}", n);
    println!("  a (Ca-K+ activation): {:.4}", a);
    println!("  b (Ca-K+ inactivation): {:.4}", b);

    // Get final currents
    let (i_na, i_k, i_k_ca, i_leak) = neuron.currents();
    println!("\nFinal ionic currents (µA/cm²):");
    println!("  I_Na: {:.4}", i_na);
    println!("  I_K: {:.4}", i_k);
    println!("  I_K(Ca): {:.4}", i_k_ca);
    println!("  I_leak: {:.4}", i_leak);

    println!("\nSimulation complete!");
}
