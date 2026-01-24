# HBT Analysis in UCC PbPb Collisions

[![Issues](https://img.shields.io/github/issues/estevamgd/HBT-in-UCC-PbPb-Collisions)](https://github.com/estevamgd/HBT-in-UCC-PbPb-Collisions/issues)
[![Last Commit](https://img.shields.io/github/last-commit/estevamgd/HBT-in-UCC-PbPb-Collisions)](https://github.com/estevamgd/HBT-in-UCC-PbPb-Collisions/commits/main)

## Overview

This repository provides a comprehensive toolkit for analyzing Hanbury-Brown and Twiss (HBT) correlations in ultra-central PbPb collisions using CMS Open Data and the ROOT framework. The analysis aims to extract spatial and temporal information about the particle-emitting source, contributing to the study of the Quark-Gluon Plasma (QGP).

## Features

- **Centrality Selection**: Tools for selecting and categorizing events by centrality or multiplicity.
- **Signal and Mixed Event Generation**: Scripts for generating signal and mixed-event distributions, including parallelized and optimized versions.
- **HBT Correlation Analysis**: Calculation of correlation functions with support for various fit models (Exponential, Gaussian, Lévy, etc.).
- **Fitting and Ratio Extraction**: Automated fitting routines and extraction of single and double ratios.
- **Validation and Benchmarking**: Utilities for validating results, comparing methods, and benchmarking performance (including timing comparisons).
- **Data Visualization**: Generation of publication-quality plots for correlation functions, source radii, timing comparisons, and more.
- **ROOT Macro Utilities**: Helper functions for histogram management, file handling, and batch processing.

## Getting Started

### Prerequisites

- [ROOT](https://root.cern/) (with multithreading support)
- C++17 or newer
- CMS Open Data (or compatible ROOT files)

### Build and Run

Most analysis scripts are C++ files or ROOT macros. Example for compiling and running:

```bash
g++ -std=c++17 -pthread src/makeSignal.cpp -o makeSignal `root-config --cflags --libs`
./makeSignal
```

Or, for ROOT macros:

```bash
root -l src/hbt_analysis_perpheralPbPb.C
```

### Directory Structure

- `src/` — Main analysis scripts and macros
- `include/` — Header files and utility functions
- `tests/` — Test scripts and validation tools
- `data/` — Input ROOT files
- `imgs/` — Output plots and images
- `benchmarks/` — Benchmarking results

## Example Analyses

- **Signal and Mixing Distributions**: `makeSignal.cpp`, `makeSignalMixOLMParallel.cpp`
- **HBT Fitting**: `fitting_sr.C`, `fitSingleRatio.cpp`
- **Timing Comparisons**: `compareTiming.cpp`, `timing_comparison_plot_TGraph.cpp`
- **Validation**: `runValidation.cpp`, `validation_func.h`

## Contributing

Contributions are welcome! Please open issues or pull requests for bug fixes, new features, or documentation improvements.
