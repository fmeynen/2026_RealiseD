# Template Project — RStudio

A lightweight template for organizing reproducible RStudio projects.

This repository is a starting template I use for new RStudio projects. It includes a recommended folder structure, a starter R project, basic scripts to get going, a lintr file for consistent and clear formatting.

Table of contents

* [Why this template](#why-this-template)
* [Project structure](#project-structure)
* [Using this repository as a template](#using-this-repository-as-a-template)
* [Using the .lintr file](#-Using-the-.lintr-file)
* [Contributing](#contributing)

## Why this template

This template enforces a clear layout and reproducible workflow for data analysis and R package-adjacent projects. It is minimal so you can adapt it quickly to specific analyses or packages.

## Project structure

This template follows the following opinionated structure:

```
│   .gitignore                       — files to be ignored
│   .lintr                           — linter file for consistent formatting
│   README.md                        — this file
│   Template\_Project\_RStudio.Rproj   — R project file, name changed automatically on first push
│   Tijd.xlsx                        — Record time spent on project
├───.github                          — Github actions
├───data                             — All contents, except for optional test data are ignored by gitignore
│   ├───processed
│   ├───raw
│   └───test
├───reports
│   ├───drafts
│   └───final
├───research\_question
│   ├───meeting\_notes
│   └───research\_papers              — Research papers pertaining to research question, provided by researchers
├───results
│   ├───graphs
│   └───tables                        — Current preference on excel files for inclusion in Word
├───scripts
│       analysis.R
│       data\_management.R
│       exploration.R
│       helpers.R
│       run\_all.R
└───supplementary\_material
    ├───R\_package\_manuals
    └───statistical\_papers           — Any papers on statistical material
```

## Using this repository as a template

To create a new project from this template:

1. Use GitHub's "Use this template" button to create a new repo (or fork/clone).
2. Update the README file.
3. Push once to automatically update the name of the R project. .github/workflows/init.yml is now automatically removed

## Using the .lintr file

I use [lintr](https://lintr.r-lib.org/) for static code analysis. Make sure to add it as an addon in RStudio

The linter file currently has the default settings except for:

* line length of maximum 120
* assignment operators <-, -> and =
* snake\_case, except for some prespecified statistical acronyms (e.g., SD, SME,...)

## Contributing

Contributions are welcome. Suggested workflow:

1. Fork the repo.
2. Create a branch (feature/your-feature).
3. Make changes and add tests where applicable.
4. Open a pull request describing the change.
