[![Build Status](https://github.com/daavid00/py-micp/actions/workflows/py-micp.yml/badge.svg)](https://github.com/daavid00/py-micp/actions/workflows/py-micp.yml)
<a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-3.6%20|%203.7%20|%203.8%20|%203.9%20|%203.10-blue.svg"></a>
[![Code style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<img src="py-micp.gif" width="900" height="250">

# py-micp: An open-source framework to perform microbially induced calcite precipitation and CO2 assessment studies at the field scale

This repository contains runscripts to simulate microbially induced calcite
precipitation (MICP) treatment in different domains; in addition to CO2 injection
assessment in leakage paths before and after MICP application. We refer to
paper [A] for an extended description of the background.

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
* [Flow](https://opm-project.org) (tested with Flow==release/2022.04)
* [Python](https://www.python.org/downloads/) (tested with Python==3.10.5)
* [GNU Octave](https://www.gnu.org/software/octave/download) (tested with GNU Octave==7.1.0)
* [MRST](https://www.sintef.no/projectweb/mrst/download/) (tested with MRST==2022a)

You will also need to install some python packages see ```requirements.txt``` 
for a complete list. You can install all the required python packages in a 
virtual environment with the following commands:

```bash
# Clone the repo
git clone https://github.com/daavid00/py-micp.git
# Get inside the folder
cd py-micp
# Create virtual environment
python3 -m venv .venv
# Activate virtual environment
source .venv/bin/activate
# Upgrade pip and setuptools
pip install --upgrade pip setuptools wheel
# Install the python requirements
pip install -r requirements.txt
```

You can install Flow from binary packages on Ubuntu Linux 16.04 or 18.04 and Red
Hat Enterprise Linux 6 or 7. Installing instruction is found 
[here](https://opm-project.org/?page_id=245). The OPM webpage also has 
instruction for installation on other systems. Once you have installed Flow, 
GNU Octave, and downloaded MRST:
* Edit line 11 of the python scripts (e.g., 'quarter_single_leak.py') with the
full path to your 'flow' executable (e.g., '~/Opm/opm-simulators/build-cmake/bin
/flow').
* Edit line 2 of the GNU Octave scripts (e.g., 'quarter_single_leak.m') with the
full path to your MRST 'startup.m' file (e.g., '~/mrst-2022a/startup.m').

## Running the scripts
* From the terminal, e.g., for the quarter_single_leak system:

`python3 quarter_single_leak.py`

## Paper
[A] Tveit, S. and Landa-Marbán, D. Field-scale optimization of injection
strategies for leakage mitigation using microbially induced calcite
precipitation. Submitted. https://doi.org/10.13140/RG.2.2.22042.16324

## Contact
David Landa-Marbán (dmar@norceresearch.no).
