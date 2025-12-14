use criterion::{black_box, criterion_group, criterion_main, Criterion};
use hodgkin_huxley::{HodgkinHuxleyNeuron, neuron_types::NeuronConfig};

fn bench_single_step(c: &mut Criterion) {
    c.bench_function("single_step", |b| {
        let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
        neuron.initialize_rest();
        b.iter(|| {
            neuron.step(black_box(0.01), black_box(10.0)).unwrap();
        });
    });
}

fn bench_action_potential(c: &mut Criterion) {
    c.bench_function("action_potential_100ms", |b| {
        b.iter(|| {
            let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
            neuron.initialize_rest();
            neuron.simulate(black_box(100.0), black_box(0.01), black_box(15.0)).unwrap();
        });
    });
}

fn bench_different_neuron_types(c: &mut Criterion) {
    let mut group = c.benchmark_group("neuron_types");

    group.bench_function("squid_axon", |b| {
        b.iter(|| {
            let mut neuron = HodgkinHuxleyNeuron::squid_axon().unwrap();
            neuron.initialize_rest();
            neuron.simulate(black_box(50.0), black_box(0.01), black_box(10.0)).unwrap();
        });
    });

    group.bench_function("regular_spiking", |b| {
        b.iter(|| {
            let mut neuron = HodgkinHuxleyNeuron::regular_spiking().unwrap();
            neuron.initialize_rest();
            neuron.simulate(black_box(50.0), black_box(0.01), black_box(10.0)).unwrap();
        });
    });

    group.bench_function("fast_spiking", |b| {
        b.iter(|| {
            let mut neuron = HodgkinHuxleyNeuron::fast_spiking().unwrap();
            neuron.initialize_rest();
            neuron.simulate(black_box(50.0), black_box(0.01), black_box(10.0)).unwrap();
        });
    });

    group.finish();
}

criterion_group!(benches, bench_single_step, bench_action_potential, bench_different_neuron_types);
criterion_main!(benches);
