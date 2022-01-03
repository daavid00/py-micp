## Description
This repository contains runscripts to simulate microbially induced calcite
precipitation treatment in different domains; in addition to CO2 injection
assessment in leakage paths before and after MICP application. They have been
successfully run in Linux and Mac OS (not tested in Windows). We refer to
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

## Requirements
* [OPM](https://opm-project.org) (tested with OPM==release/2021.10)
* [Python](https://www.python.org/downloads/) (tested with Python3==3.8.12)
* [GNU Octave](https://www.gnu.org/software/octave/download) (tested with GNU Octave==6.4.0)

## GNU Octave dependency
* [MRST](https://www.sintef.no/projectweb/mrst/download/) (tested with MRST==2021b)

## Python dependencies
* [Mako](https://www.makotemplates.org) (tested with Mako==1.1.4)
* [matplotlib](https://matplotlib.org) (tested with matplotlib==3.3.2)
* [meshio](https://github.com/nschloe/meshio) (tested with meshio==4.3.12)
* [numpy](https://numpy.org) (tested with numpy==1.19.3)
* [os](https://docs.python.org/3/library/os.html) (tested with os==3.8.12)

## Installation
* Follow the instructions in the links to install requirements and dependencies.
* Edit line 11 of the python scripts (e.g., 'quarter_single_leak.py') with the
full path to your 'flow' executable (e.g., '~/Opm/opm-simulators/build-cmake/bin
/flow').
* Edit line 2 of the GNU Octave scripts (e.g., 'quarter_single_leak.m') with the
full path to your MRST 'startup.m' file (e.g., '~/mrst-2021b/startup.m').

## Running the scripts
* From the terminal, e.g., for the quarter_single_leak system:

`python3 quarter_single_leak.py`

## Paper
* [A] Tveit, S. and Landa-Marbán, D. Field-scale optimization of injection
strategies for leakage mitigation using microbially induced calcite
precipitation. Submitted. https://doi.org/10.13140/RG.2.2.22042.16324

## Contact
David Landa-Marbán (dmar@norceresearch.no).
