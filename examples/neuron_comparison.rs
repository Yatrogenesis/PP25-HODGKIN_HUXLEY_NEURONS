//! Example: Comparing different neuron types.
//!
//! This example demonstrates the different firing patterns of:
//! - Regular spiking (RS) neurons
//! - Fast spiking (FS) neurons
//! - Intrinsically bursting (IB) neurons

use hodgkin_huxley::HodgkinHuxleyNeuron;

fn main() {
    println!("Neuron Type Comparison");
    println!("======================\n");

    let neuron_types = vec![
        ("Squid Axon", HodgkinHuxleyNeuron::squid_axon().unwrap()),
        ("Regular Spiking", HodgkinHuxleyNeuron::regular_spiking().unwrap()),
        ("Fast Spiking", HodgkinHuxleyNeuron::fast_spiking().unwrap()),
        ("Intrinsically Bursting", HodgkinHuxleyNeuron::intrinsically_bursting().unwrap()),
    ];

    let dt = 0.01;
    let duration = 100.0;
    let current = 10.0;

    for (name, mut neuron) in neuron_types {
        neuron.initialize_rest();

        // Simulate
        neuron.simulate(duration, dt, current)
            .expect("Simulation failed");

        // Detect spikes
        let spikes = neuron.detect_spikes(-20.0);
        let firing_rate = HodgkinHuxleyNeuron::firing_rate(&spikes);

        println!("{}", name);
        println!("  Resting voltage: {:.2} mV", neuron.voltage());
        println!("  Number of spikes: {}", spikes.len());
        println!("  Firing rate: {:.2} Hz", firing_rate);

        if spikes.len() > 1 {
            let isis = HodgkinHuxleyNeuron::interspike_intervals(&spikes);
            let mean_isi = isis.iter().sum::<f64>() / isis.len() as f64;
            let std_isi = (isis.iter().map(|x| (x - mean_isi).powi(2)).sum::<f64>() / isis.len() as f64).sqrt();
            println!("  Mean ISI: {:.2} ms (± {:.2})", mean_isi, std_isi);
        }

        println!();
    }
}
