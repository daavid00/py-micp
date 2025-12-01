# Setting up and simulating a 3D flow system with a straight leakage path, where
# two different injection strategies for MICP treatment are compared, where on
# the first one all solutions are injected only using the top of the well and on
# the second one the injection of the growth solution is through a lower section
# of the well. Information about the model and parameters can be found in:

# Landa-Marbán, D., Kumar, K., Tveit, S., Gasda, S.E. Numerical studies of CO2
# leakage remediation by micp-based plugging technology. In: Røkke, N.A. and
# Knuutila, H.K. (Eds) Short Papers from the 11th International Trondheim CCS
# conference, ISBN: 978-82-536-1714-5, 284-290.

# Set the full path to the flow executable and flags
flow = "flow --initial-time-step-in-days=0.0001 --enable-opm-rst-file=true"
# The output dir and vtk needs to be defined for the files to be postprocess
flow += " --output-dir=results --enable-vtk-output=true"

# Import python dependencies
import os
import numpy as np
from mako.template import Template

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
F = 0.5  # Oxygen consumption factor [-]
Y = 0.5  # Yield growth coefficient [-]
muw = 2.54e-4  # Water viscoscity [Pa s]
Ka = 1e-14  # Aquifer permeability [m^2]
Kl = 2e-14  # Leakage permeability [m^2]
phi = 0.15  # Porosity [-]
md = 9.87e-16  # From m^2 to milli darcy [mD]
yurca = 1.67  # Yielc urea to calcite [-]
detexp = 0.58  # Exponent for the norm of the detahcment rate [-]
porofac = [0, 0.4, 0.6, 0.9, 1]  # See PERMFACT in the OPM Flow manual
permfac = [0.005, 0.005, 0.01, 0.1, 1]  # See PERMFACT table in the OPM Flow manual

# Delete previous inputs and results
os.system("rm -rf decks & rm -rf results & wait")
os.system("mkdir decks & mkdir results & wait")

# Create the grid
H = 30  # Height of the domain [m]
ht = 5  # Height of the top aquifer [m]
hl = 5  # Height of the lower aquifer [m]
L = 100  # Length of the domain [m]
ux = 14  # X-gap between aquifer and leakage [m]
uy = 14  # Y-gap between aquifer and leakage [m]
Wi = 100  # Width of the domain [m]
Z0 = np.linspace(0, ht, ht + 1)  # Discretization of the z-dir Z0
Z1 = np.linspace(ht + 1.0, H - hl + 1.0, H - hl - ht + 1)  # ""     Z1
Z2 = np.linspace(H - hl + 2.0, H, 4)  # ""                          Z2
Wc = 1  # Cells for each MICP well injection
Xf0 = np.linspace(0, -1.5, 16)  # ""                         -L*exp(Xf0)
Xf1 = np.linspace(-21, 21, 44)  # Discretization of the x-dir Xf1
Xf2 = np.linspace(-1.5, 0, 16)  # ""                          L*exp(Xf2)
Yf0 = np.linspace(0, -1.5, 16)  # ""                         -Wi*exp(Yf0)
Yf1 = np.linspace(-21, 21, 44)  # Discretization of the y-dir Yf1
Yf2 = np.linspace(-1.5, 0, 16)  # ""                          Wi*exp(Yf2)
XXf, YYf, XXq, YYq, ZZ = [], [], [], [], []
for k in range(len(Z0)):
    ZZ.append(float(Z0[k]))
for k in range(len(Z1)):
    ZZ.append(float(Z1[k]))
for k in range(len(Z2)):
    ZZ.append(float(Z2[k]))
for i in range(len(Xf0)):
    XXf.append(float(-L * np.exp(Xf0[i])))
for i in range(len(Xf1)):
    XXf.append(float(Xf1[i]))
for i in range(len(Xf2)):
    XXf.append(float(L * np.exp(Xf2[i])))
for j in range(len(Yf0)):
    YYf.append(float(-Wi * np.exp(Yf0[j])))
for j in range(len(Yf1)):
    YYf.append(float(Yf1[j]))
for j in range(len(Yf2)):
    YYf.append(float(Wi * np.exp(Yf2[j])))
xxf = len(Xf0) + len(Xf1) + len(Xf2) - 1
yyf = len(Yf0) + len(Yf1) + len(Yf2) - 1
zz = len(Z0) + len(Z1) + len(Z2) - 1
xyzf = xxf * yyf * zz
ix = np.abs(np.transpose(XXf)).argmin() + 1
iy = np.abs(np.transpose(YYf)).argmin() + 1
it = len(np.transpose(ZZ)[tuple([np.transpose(ZZ) < ht])])
il = len(np.transpose(ZZ)[tuple([np.transpose(ZZ) > H - hl])])
os.system(
    f"octave --eval 'full_single_leak({H}, {ht}, {hl}, {L}, {ux}, "
    + f"{uy}, {Wi}, {XXf}, {YYf}, {ZZ})' "
)

# Create variables used in this script
permf = [Ka / md] * xyzf
Lg, t0, t1, injestra = [], [], [], []
poro = phi * np.ones(xyzf, dtype=float)


def compact_format(values):
    """
    Use the 'n*x' notation to write repited values to save storage

    Args:
        values (list): List with the variable values

    Returns:
        values (list): List with the compacted variable values

    """
    n, value0, tmp = 0, float(values[0]), []
    for value in values:
        if value0 != float(value) or len(values) == 1:
            if value0 == 0:
                tmp.append(f"{n}*0 " if n > 1 else "0 ")
            elif value0.is_integer():
                tmp.append(f"{n}*{int(value0)} " if n > 1 else f"{int(value0)} ")
            else:
                tmp.append(f"{n}*{value0} " if n > 1 else f"{value0} ")
            n = 1
            value0 = float(value)
        else:
            n += 1
    if value0 == float(values[-1]) and len(values) > 1:
        if value0 == 0:
            tmp.append(f"{n}*0 " if n > 1 else "0 ")
        elif value0.is_integer():
            tmp.append(f"{n}*{int(value0)} " if n > 1 else f"{int(value0)} ")
        else:
            tmp.append(f"{n}*{value0} " if n > 1 else f"{value0} ")
    return tmp


# Define injection strategy
# [injection time [h], injection rate [m^3/s], microbial concentration [kg/m^3],
# oxygen concentration [kg/m^3], and urea concentration [kg/m^3]]. Before each
# injection of the solutions in the phases we inject only water for a very short
# time [0.01 hour] to ease the convergence of the solution.
injestra.append([0.01, 8e-2, 0, 0, 0])
injestra.append([15.0, 8e-2, 0.01, 0, 0])
injestra.append([2.0, 8e-2, 0, 0, 0])
injestra.append([100.0, 0, 0, 0, 0])
injestra.append([0.01, 8e-2, 0, 0, 0])
injestra.append([30.0, 8e-2, 0, 0.04, 0])
injestra.append([2.0, 8e-2, 0, 0, 0])
injestra.append([30.0, 0, 0, 0, 0])
injestra.append([0.01, 8e-2, 0, 0, 0])
injestra.append([30.0, 8e-2, 0, 0, 60.0])
injestra.append([2.0, 8e-2, 0, 0, 0])
injestra.append([30.0, 0, 0, 0, 0])

# Set and run the full and quarter domain MICP simulations
permf[xxf * yyf * it : xxf * yyf * (zz - il)] = [Kl / md] * (xxf * yyf * (zz - il - it))
permw = compact_format("".join(f"{val} " for val in permf).split())
porow = compact_format("".join(f"{val} " for val in poro).split())
var = {
    "porofac": porofac,
    "permfac": permfac,
    "detexp": detexp,
    "yurca": yurca,
    "F": F,
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
    "Wc": Wc,
    "xx": xxf,
    "Y": Y,
    "yy": yyf,
    "zz": zz,
}
mytemplate = Template(filename="top.mako")
FilledTemplate = mytemplate.render(**var)
with open("decks/TOP.DATA", "w") as f:
    f.write(FilledTemplate)
mytemplate = Template(filename="variable.mako")
FilledTemplate = mytemplate.render(**var)
with open("decks/VARIABLE.DATA", "w") as f:
    f.write(FilledTemplate)
os.system(f"{flow} decks/TOP.DATA")

# {
# Copyright 2021-2025, NORCE Research AS, Computational
# Geosciences and Modeling.

# This file is part of the py-micp module.

# py-micp is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# py-micp is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
# }
