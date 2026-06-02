#!/usr/bin/env python
# SPDX-FileCopyrightText: 2021-2026 NORCE Research AS
# SPDX-License-Identifier: GPL-3.0
# pylint: disable=C0116,R0913,R0914,R0915,R0917,R0801

"""
Setting up and simulating a 3D flow system with a straight leakage path, where
the horizontal direction of the injected microbial solution is contoled by
adding injection wells on the opposite corner boundaries to account for the
coming fluid from the boundaries. Information about the model and parameters
can be found in:

Landa-Marbán, D., Kumar, K., Tveit, S., Gasda, S.E. Numerical studies of CO2
leakage remediation by micp-based plugging technology. In: Røkke, N.A. and
Knuutila, h.K. (Eds) Short Papers from the 11th International Trondheim CCS
conference, ISBN: 978-82-536-1714-5, 284-290
"""

import shutil
import subprocess
from pathlib import Path
import numpy as np
from mako.template import Template

# Set the full path to the flow executable and flags
flow = "flow --initial-time-step-in-days=0.0001 --enable-opm-rst-file=true"
# The output dir and vtk needs to be defined for the files to be postprocess
flow += " --output-dir=results --enable-vtk-output=true"


def run_command(command):
    subprocess.run(command, shell=True, check=True)


def ensure_clean_directory(path_name):
    path = Path(path_name)
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def render_template_to_file(template_name, output_path, var):
    mytemplate = Template(filename=template_name)
    filled_template = mytemplate.render(**var)
    with open(output_path, "w", encoding="utf-8") as file_handle:
        file_handle.write(filled_template)


def compact_format(v):
    change_idx = np.flatnonzero(np.diff(v, prepend=v[0] - 1))
    counts = np.diff(np.append(change_idx, v.size))
    vals = v[change_idx]

    out = []
    for val, n in zip(vals, counts):
        val_str = str(int(val)) if float(val).is_integer() else str(val)
        out.append(f"{n}*{val_str} " if n > 1 else f"{val_str} ")
    return out


def towards_more_optimal_injection_strategies():
    # Set the model parameters
    rhob = 35.0  # Density (biofilm) [kg/m^3]
    rhoc = 2710.0  # Density (calcite) [kg/m^3]
    rhow = 1045  # Density (water) [kg/m^3]
    kstr = 2.6e-10  # Detachment rate [m/(Pa s)]
    ko = 2e-5  # Half-velocity coefficient (oxygen) [kg/m^3]
    ku = 21.3  # Half-velocity constant (urea) [kg/m^3]
    mu = 4.17e-5  # Maximum specific growth rate [1/s]
    muu = 0.0161  # Maximum rate of urease utilization [1/s]
    ka = 8.51e-7  # Microbial attachment rate [1/s]
    kd = 3.18e-7  # Microbial death rate [1/s]
    f = 0.5  # Oxygen consumption factor [-]
    yld = 0.5  # Yield growth coefficient [-]
    muw = 2.54e-4  # Water viscoscity [Pa s]
    k = 1e-14  # Aquifer permeability [m^2]
    kl = 2e-14  # Leakage permeability [m^2]
    phi = 0.15  # Porosity [-]
    md = 9.87e-16  # From m^2 to milli darcy [mD]
    yurca = 1.67  # Yielc urea to calcite [-]
    detexp = 0.58  # Exponent for the norm of the detahcment rate [-]
    porofac = [0, 0.4, 0.6, 0.9, 1]  # See PERMFACT in the OPM Flow manual
    permfac = [0.005, 0.005, 0.01, 0.1, 1]  # See PERMFACT table in the OPM Flow manual

    # Delete previous inputs and results
    ensure_clean_directory("decks")
    ensure_clean_directory("results")

    # Create the grid
    h = 30  # Height of the domain [m]
    ht = 5  # Height of the top aquifer [m]
    hl = 5  # Height of the lower aquifer [m]
    le = 100  # Length of the domain [m]
    ux = 14  # x-gap between aquifer and leakage [m]
    uy = 14  # y-gap between aquifer and leakage [m]
    wi = 100  # Width of the domain [m]
    z0 = np.linspace(0, ht, ht + 1)  # Discretization of the z-dir z0
    z1 = np.linspace(ht + 1.0, h - hl + 1.0, h - hl - ht + 1)  # ""     z1
    z2 = np.linspace(h - hl + 2.0, h, 4)  # ""                          z2
    wc = 1  # Cells for each MICP well injection
    xf0 = np.linspace(0, -1.5, 16)  # ""                         -le*exp(xf0)
    xf1 = np.linspace(-21, 21, 44)  # Discretization of the x-dir xf1
    xf2 = np.linspace(-1.5, 0, 16)  # ""                          le*exp(xf2)
    yf0 = np.linspace(0, -1.5, 16)  # ""                         -wi*exp(yf0)
    yf1 = np.linspace(-21, 21, 44)  # Discretization of the y-dir yf1
    yf2 = np.linspace(-1.5, 0, 16)  # ""                          wi*exp(yf2)
    z = np.concatenate((z0, z1, z2)).astype(float)
    xf = np.concatenate((-le * np.exp(xf0), xf1, le * np.exp(xf2))).astype(float)
    yf = np.concatenate((-wi * np.exp(yf0), yf1, wi * np.exp(yf2))).astype(float)
    xxf = len(xf0) + len(xf1) + len(xf2) - 1
    yyf = len(yf0) + len(yf1) + len(yf2) - 1
    zz = len(z0) + len(z1) + len(z2) - 1
    xyzf = xxf * yyf * zz
    ix = int(np.abs(xf).argmin()) + 1
    iy = int(np.abs(yf).argmin()) + 1
    it = int(np.count_nonzero(z < ht))
    il = int(np.count_nonzero(z > h - hl))
    run_command(
        f"octave --eval 'full_single_leak({h}, {ht}, {hl}, {le}, {ux}, "
        f"{uy}, {wi}, {xf.tolist()}, {yf.tolist()}, {z.tolist()})'"
    )

    # Create variables used in this script
    permf = np.full(xyzf, k / md, dtype=float)
    injestra = []
    poro = np.full(xyzf, phi, dtype=float)

    if not Path("micp.mako").is_file():
        raise FileNotFoundError("Required file not found: micp.mako")

    # Define injection strategy
    # [injection time [h], injection rate [m^3/s], microbial concentration [kg/m^3],
    # oxygen concentration [kg/m^3], and urea concentration [kg/m^3]].
    injestra.append([0.01, 2.0e-02, 0.00, 0.00, 0.0])
    injestra.append([15.0, 2.0e-02, 0.00, 0.00, 60.0])
    injestra.append([10.0, 2.0e-02, 0.00, 0.00, 0.0])

    # Set and run the full and quarter domain MICP simulations
    permf[xxf * yyf * it : xxf * yyf * (zz - il)] = kl / md
    permw = compact_format(permf)
    porow = compact_format(poro)
    var = {
        "porofac": porofac,
        "permfac": permfac,
        "detexp": detexp,
        "yurca": yurca,
        "f": f,
        "il": il,
        "ix": ix,
        "iy": iy,
        "injestra": injestra,
        "ka": ka,
        "kd": kd,
        "ko": ko,
        "kstr": kstr,
        "ku": ku,
        "mu": mu,
        "muu": muu,
        "muw": muw,
        "perm": permw,
        "poro": porow,
        "rhob": rhob,
        "rhoc": rhoc,
        "rhow": rhow,
        "wc": wc,
        "xx": xxf,
        "it": it,
        "yld": yld,
        "yy": yyf,
        "zz": zz,
    }
    render_template_to_file("micp.mako", "decks/DIRECTIONAL.DATA", var)
    run_command(f"{flow} decks/DIRECTIONAL.DATA")


if __name__ == "__main__":
    towards_more_optimal_injection_strategies()
