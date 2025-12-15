---
title: 'hodgkin-huxley-rs: A Production-Ready Rust Implementation of Biophysically Accurate Neuron Models'
tags:
  - Rust
  - neuroscience
  - computational biology
  - Hodgkin-Huxley model
  - action potential
  - neural simulation
authors:
  - name: Francisco Molina Burgos
    orcid: 0009-0008-6093-8267
    affiliation: 1
    corresponding: true
affiliations:
  - name: Independent Researcher, Spain
    index: 1
date: 14 December 2025
bibliography: paper.bib
repository: https://github.com/Yatrogenesis/PP25-HODGKIN_HUXLEY_NEURONS
archive_doi: 10.5281/zenodo.17932872
---

# Summary

`hodgkin-huxley-rs` is an open-source Rust library that implements the Hodgkin-Huxley neuron model with exact biophysical equations from the Nobel Prize-winning 1952 paper [@Hodgkin1952]. The library enables researchers and educators to simulate action potential generation in neurons with high accuracy and performance. It provides six physiologically distinct neuron types commonly studied in neuroscience: squid giant axon (the original model), regular spiking cortical pyramidal neurons, fast spiking interneurons, intrinsically bursting neurons, low-threshold spiking interneurons, and chattering neurons. The software includes temperature-dependent kinetics via Q10 scaling, multiple ion channel types (sodium, potassium, calcium-activated potassium, and leak), and two numerical integration schemes optimized for different accuracy-speed tradeoffs.

# Statement of need

The Hodgkin-Huxley model remains the gold standard for understanding how neurons generate electrical signals. Despite its fundamental importance in neuroscience education and research, existing implementations suffer from several limitations: legacy codebases in C or FORTRAN lack memory safety guarantees, many implementations are restricted to a single neuron type, and few provide the temperature dependence necessary for modeling mammalian neurons.

`hodgkin-huxley-rs` addresses these gaps by providing:

1. **Memory safety**: Rust's ownership model eliminates buffer overflows and use-after-free errors that plague legacy scientific code, making large-scale simulations more reliable.

2. **Multiple neuron phenotypes**: Six pre-configured neuron types with parameters derived from experimental literature [@Connor1971; @Traub1991], enabling comparative studies of neuronal diversity.

3. **Temperature dependence**: Q10 scaling allows accurate simulation of neurons at physiological temperatures (37°C) rather than the original squid axon temperature (6.3°C).

4. **Production-ready error handling**: Comprehensive `Result` types enable robust integration into larger simulation frameworks and pipelines.

5. **High performance**: Integration steps complete in 1–2 μs on commodity hardware, enabling real-time simulation of networks with approximately 10,000 neurons on a single CPU core.

The library targets computational neuroscientists building large-scale neural network simulations, educators teaching electrophysiology, and researchers requiring validated implementations of classical neuron models. It complements existing tools like NEURON [@Hines1997] and Brian2 [@Stimberg2019] by providing a lightweight, embeddable alternative with minimal dependencies.

# Implementation

The membrane potential $V$ evolves according to the current-balance equation:

$$C_m \frac{dV}{dt} = -I_{Na} - I_K - I_{K(Ca)} - I_{leak} + I_{ext}$$

where ionic currents follow the standard Hodgkin-Huxley formalism with voltage-dependent gating variables. The library implements both fourth-order Runge-Kutta (RK4) and exponential Euler integrators. The exponential Euler method exploits the linear structure of gating variable equations, providing unconditional stability for time steps up to 0.1 ms.

Table 1 shows the conductance parameters for each implemented neuron type, derived from experimental measurements.

| Type | $g_{Na}$ (mS/cm²) | $g_K$ (mS/cm²) | $g_{K(Ca)}$ (mS/cm²) |
|------|-------------------|----------------|----------------------|
| Squid Axon | 120 | 36 | 0 |
| Regular Spiking | 100 | 30 | 0.3 |
| Fast Spiking | 150 | 50 | 0.05 |
| Intrinsically Bursting | 120 | 25 | 0.8 |
| Low-Threshold Spiking | 90 | 28 | 0.2 |
| Chattering | 140 | 40 | 1.2 |

: Conductance parameters for implemented neuron types. \label{tab:neurons}

# Validation

The implementation reproduces Figure 12 from Hodgkin & Huxley (1952) within 0.5 mV peak deviation and 0.1 ms timing accuracy. Firing rate curves (f-I relationships) match published electrophysiological data for each neuron type.

# Acknowledgements

The author thanks the Rust scientific computing community for foundational work on numerical libraries.

# References
