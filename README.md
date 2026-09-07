# 2026_RealiseD

An R-based research project with choice-based conjoint (CBC) analysis and statistical modeling components.

## Overview

This repository contains reproducible research code and analysis for investigating CBC (Choice-Based Conjoint) estimation methods. It integrates LaTeX documentation with R-based data analysis workflows, following a structured project organization.

## Project Structure

```
│   .gitignore                      — Files to ignore during version control
│   .lintr                          — R linter configuration for code style
│   README.md                       — This file
│   2026_RealiseD.Rproj            — RStudio project file
│   CBCEstimator.tex               — LaTeX source documentation
│   CBCEstimator.pdf               — Compiled PDF documentation
│   CBCEstimator.log               — LaTeX compilation log
│   Tijd.xlsx                       — Time tracking spreadsheet
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
├───results                        — Output from analyses
│   ├───graphs                     — Figures and plots
│   └───tables                     — Summary tables and results
│
├───reports                        — Reports and manuscripts
│   ├───drafts                     — Work-in-progress versions
│   └───final                      — Final report/manuscript versions
│
├───research_question              — Research documentation
│   ├───meeting_notes              — Notes from research meetings
│   └───research_papers            — Reference papers and literature
│
└───supplementary_material         — Additional resources
    ├───R_package_manuals          — Documentation for R packages used
    └───statistical_papers        — Papers on statistical methods
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

The main documentation is provided in `CBCEstimator.tex` and compiled to `CBCEstimator.pdf`. This contains detailed information about the CBC estimation methodology and analysis approach.

## Workflow

The standard analysis workflow is:

1. **Data Management** (`data_management.R`) — Load and prepare data
2. **Exploration** (`exploration.R`) — Explore data structure and distributions
3. **Analysis** (`analysis.R`) — Run main statistical analyses
4. **Results** — Outputs saved to `results/` directory

Run all steps at once with `scripts/run_all.R`, or execute individual scripts as needed.

## Time Tracking

Project time is tracked in `Tijd.xlsx` for project management and estimation purposes.

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
