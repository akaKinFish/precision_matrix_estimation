# Precision Matrix Estimation

MATLAB implementation of frequency-domain precision matrix estimation with cross-frequency smoothing.

## Getting Started

1. Run `startup.m` to configure paths
2. Run `main.m` to execute the main pipeline
3. Run `run_all_tests.m` to run all unit tests

## Project Structure

- `src/` - Source code
  - `modules/` - Individual algorithm modules
  - `core/` - Core algorithms
  - `visualization/` - Plotting functions
- `tests/` - Test files
- `data/` - Data storage
- `results/` - Output results
- `docs/` - Documentation

## Modules

1. **Module 1**: Data preprocessing and whitening
2. **Module 2**: E-step
3. **Module 3**: Active set selection
4. **Module 4**: Gradient computation
5. **Module 5**: Proximal updates
6. **Module 6**: Hyperparameter configuration
7. **Module 7**: Simulation data generation
8. **Module 8**: Recoloring
