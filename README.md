[![Build Status](https://github.com/daavid00/py-micp/actions/workflows/CI.yml/badge.svg)](https://github.com/daavid00/py-micp/actions/workflows/CI.yml)
<a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-3.8%20to%203.14-blue.svg"></a>
[![Code style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/444088296.svg)](https://doi.org/10.5281/zenodo.16880620)
<img src="py-micp.gif" width="830" height="250">

# py-micp: An Open-Source Simulation Workflow for Field-Scale Application of Microbially Induced Calcite Precipitation Technology for Leakage Remediation

This repository contains runscripts to simulate microbially induced calcite
precipitation (MICP) treatment in different domains; in addition to CO2 injection
assessment in leakage paths before and after MICP application. We refer to 
[this paper](https://doi.org/10.1016/j.ijggc.2021.103256) for an extended description of the mathematical model.

The numerical examples accompanying this module are:
* full_diagonal_leak,
* full_double_leak,
* full_single_leak,
* quarter_diagonal_leak,
* quarter_double_leak,
* quarter_single_leak,
* top_vs_vertical_variable_injection, and
* towards_more_optimal_injection_strategies.

The scripts in these examples can be easily modified to run different studies
changing input parameter values, spatial domain properties, and implemented
models in OPM Flow.

## Installation
You will first need to install/download
* [Flow](https://opm-project.org) (tested with Flow==2025.10)
* [Python](https://www.python.org/downloads/) (tested with Python==3.14)
* [GNU Octave](https://www.gnu.org/software/octave/download) (tested with GNU Octave==10.3.0)
* [MRST](https://www.sintef.no/projectweb/mrst/download/) (tested with MRST==2025a)

You will also need to install some python packages, see ```requirements.txt``` 
for a complete list. You can install all the required python packages in a 
virtual environment with the following commands:

```bash
# Clone the repo
git clone https://github.com/daavid00/py-micp.git
# Get inside the folder
cd py-micp
# Create virtual environment
python3 -m venv vpy-micp
# Activate virtual environment
source vpy-micp/bin/activate
# Upgrade pip and setuptools
pip install --upgrade pip setuptools wheel
# Install the python requirements
pip install -r requirements.txt
```

See this [_installation_](https://opm.github.io/pyopmspe11/installation.html#opm-flow) for further details on building OPM Flow from the master branches in Linux, Windows (via [_WSL_](https://learn.microsoft.com/en-us/windows/wsl/)), and macOS.

Once you have installed OPM Flow, GNU Octave, and downloaded MRST:

* Edit line 11 of the python scripts (e.g., 'quarter_single_leak.py') with the
full path to your 'flow' executable (e.g., '/Users/dmar/opm/build/opm-simulators/bin/flow') and simulator flags.

* Edit line 2 of the GNU Octave scripts (e.g., 'quarter_single_leak.m') with the
full path to your MRST 'startup.m' file (e.g., '/Users/dmar/mrst-2025a/startup.m').

Tip: See the [CI.yml](https://github.com/daavid00/py-micp/blob/main/.github/workflows/CI.yml) for installation of 
all dependencies in Ubuntu 24.04 and running the test in _py-micp_.

## Running the scripts
* From the terminal, e.g., for the quarter_single_leak system:

`python3 quarter_single_leak.py`

## Citing
* Landa-Marbán, D. 2025. py-micp: Open-source code to perform studies of MICP treatment and CO2 assessment. https://doi.org/10.5281/zenodo.16880621.

## Publication
The following is a manuscript in which _py-micp_ is used:

1. Tveit, S. and Landa-Marbán, D., 2022. Field-scale optimization of injection
strategies for leakage mitigation using microbially induced calcite
precipitation. https://arxiv.org/abs/2201.00669

Refer to the [_v2024.10 release_](https://github.com/daavid00/py-micp/releases/tag/v2024.10) for the code that was used in that paper. This is relevant 
since in the OPM Flow version 2025.04, the MICP implementation was refactored, resulting in changes in the deck keywords (e.g., now we use `BIOFPARA` instead 
of `MICPPARA`); see the corresponding OPM Flow version of the [_OPM Flow manual_](https://opm-project.org/?page_id=955) for details.

## About py-micp
The _py-micp_ package was funded by the [_Eﬃcient models for microbially induced calcite precipitation as a seal for CO2 storage (MICAP) project_](https://gassnova.no/app/uploads/sites/4/2022/02/Sluttrapport-MICAP.pdf) [project number 268390].
Contributions are more than welcome using the fork and pull request approach.
For new features, please request them raising an issue.

## Related
Below are some tools that might be of interest; check ‘em out 🙂.

* ad-micp: A module to study CO2 leakage remediation by microbially induced calcite precipitation (MICP) (https://github.com/daavid00/ad-micp).
* pyopmspe11: A Python framework using OPM Flow for the CSP SPE11 benchmark project (https://github.com/OPM/pyopmspe11).
* pycopm: An open-source tool to tailor OPM Flow geological models (https://github.com/cssr-tools/pycopm).
* plopm: Quick generation of PNGs, GIFs, and VTKs from a OPM Flow type model (https://github.com/cssr-tools/plopm).
* pofff: An open-source image-based history-matching framework for the FluidFlower benchmark study using OPM Flow (https://github.com/cssr-tools/pofff).
* pyopmnearwell: A Python framework to simulate near well dynamics using OPM Flow (https://github.com/cssr-tools/pyopmnearwell).
* expreccs: A Python framework using OPM Flow to simulate regional and site reservoirs for CO2 storage (https://github.com/cssr-tools/expreccs).
* pymm: An open-source image-based framework for CFD in microsystems (https://github.com/cssr-tools/pymm).

## Contact
David Landa-Marbán (dmar@norceresearch.no).
