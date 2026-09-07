# 2026_RealiseD

An R-based research project comparing the performance of three different types of analysis in the context of rare diseases.

## Overview

The three methods compared are:
* Classical Maximum Likelihood with the lme4 package
* Multiple Imputation with mice on the cluster level followed by a closed form cluster-by-cluster estimator (CBC-estimator)
* The closed form CBC-estimator followed by iterative reweighting

## Project Structure

```
│   .gitignore                      — Files to ignore during version control
│   .lintr                          — R linter configuration for code style
│   README.md                       — This file
│   2026_RealiseD.Rproj            — RStudio project file
│
├───data                            — Data storage (contents ignored by gitignore except tests)
│   ├───raw                        — Original, unmodified data files
│   ├───processed                  — Cleaned and processed data
│   └───test                       — Test datasets
│
├───scripts                        — R analysis scripts
│   ├───data_management.R          — Data loading, cleaning, and preprocessing
│   ├───exploration.R              — Exploratory data analysis
│   ├───analysis.R                 — Main statistical analysis
│   ├───helpers.R                  — Utility functions
│   └───run_all.R                  — Master script to run full pipeline
│
└───results                        — Output from analyses
    ├───graphs                     — Figures and plots
    └───tables                     — Summary tables and results
 
```

## Getting Started

### Prerequisites

- R (latest version recommended)
- RStudio (optional, but recommended)
- Required R packages (see scripts for specific dependencies)

### Setup

1. Clone this repository
2. Open `2026_RealiseD.Rproj` in RStudio
3. Install any required packages listed in the scripts
4. Run `scripts/run_all.R` to execute the complete analysis pipeline

### Code Style

This project uses [lintr](https://lintr.r-lib.org/) for static code analysis. Configuration includes:

- Maximum line length: 120 characters
- Assignment operators: `<-`, `->`, and `=`
- Naming convention: `snake_case` (with exceptions for statistical acronyms like SD, SME, etc.)

## Documentation



## Workflow

The standard analysis workflow is:

1. **Data Management** (`data_management.R`) — Load and prepare data
2. **Exploration** (`exploration.R`) — Explore data structure and distributions
3. **Analysis** (`analysis.R`) — Run main statistical analyses
4. **Results** — Outputs saved to `results/` directory

Run all steps at once with `scripts/run_all.R`, or execute individual scripts as needed.

## Contributing

Contributions are welcome. Suggested workflow:

1. Fork the repository
2. Create a feature branch (`feature/your-feature`)
3. Make changes and test thoroughly
4. Open a pull request with a clear description of changes

## License

See LICENSE file for details (if applicable).

## Contact

For questions about this project, please open an issue on GitHub.
