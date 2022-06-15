# Setting up and simulating a 3D flow system with a straight leakage path, where
# the horizontal direction of the injected microbial solution is contoled by
# adding injection wells on the opposite corner boundaries to account for the
# coming fluid from the boundaries. Information about the model and parameters
# can be found in:

# Landa-Marbán, D., Kumar, K., Tveit, S., Gasda, S.E. Numerical studies of CO2
# leakage remediation by micp-based plugging technology. In: Røkke, N.A. and
# Knuutila, H.K. (Eds) Short Papers from the 11th International Trondheim CCS
# conference, ISBN: 978-82-536-1714-5, 284-290.

# Set the full path to the flow executable
flow = "~/Opm/master/build/opm-simulators/bin/flow"

# Import python dependencies
import os
import meshio
import numpy as np
from mako.template import Template
from matplotlib import pyplot as plt

# Set the model parameters
rhob = 35.0  # Density (biofilm) [kg/m^3]
rhoc = 2710.0  # Density (calcite) [kg/m^3]
rhow = 1045  # Density (water) [kg/m^3]
kstr = 2.6e-10  # Detachment rate [m/(Pa s)]
crit = 0.1  # Critical porosity [-]
eta = 3.0  # Fitting factor [-]
ko = 2e-5  # Half-velocity coefficient (oxygen) [kg/m^3]
ku = 21.3  # Half-velocity constant (urea) [kg/m^3]
mu = 4.17e-5  # Maximum specific growth rate [1/s]
muu = 0.0161  # Maximum rate of urease utilization [1/s]
ka = 8.51e-7  # Microbial attachment rate [1/s]
kd = 3.18e-7  # Microbial death rate [1/s]
kmin = 1e-20  # Minimum permeability [m^2]
F = 0.5  # Oxygen consumption factor [-]
Y = 0.5  # Yield growth coefficient [-]
muw = 2.54e-4  # Water viscoscity [Pa s]
Ka = 1e-14  # Aquifer permeability [m^2]
Kl = 2e-14  # Leakage permeability [m^2]
phi = 0.15  # Porosity [-]
md = 9.87e-16  # Milli darcy [mD]
tolclg = 1e-4  # Tolerance before clogging to stop the simulation

# Delete previous inputs and results
os.system("rm -rf inputs & rm -rf results & wait")
os.system("mkdir inputs & mkdir results & wait")

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
XXf = []
YYf = []
XXq = []
YYq = []
ZZ = []
for k in range(len(Z0)):
    ZZ.append(Z0[k])
for k in range(len(Z1)):
    ZZ.append(Z1[k])
for k in range(len(Z2)):
    ZZ.append(Z2[k])
for i in range(len(Xf0)):
    XXf.append(-L * np.exp(Xf0[i]))
for i in range(len(Xf1)):
    XXf.append(Xf1[i])
for i in range(len(Xf2)):
    XXf.append(L * np.exp(Xf2[i]))
for j in range(len(Yf0)):
    YYf.append(-Wi * np.exp(Yf0[j]))
for j in range(len(Yf1)):
    YYf.append(Yf1[j])
for j in range(len(Yf2)):
    YYf.append(Wi * np.exp(Yf2[j]))
xxf = len(Xf0) + len(Xf1) + len(Xf2) - 1
yyf = len(Yf0) + len(Yf1) + len(Yf2) - 1
zz = len(Z0) + len(Z1) + len(Z2) - 1
xyzf = xxf * yyf * zz
ix = abs(np.transpose(XXf)).argmin() + 1
iy = abs(np.transpose(YYf)).argmin() + 1
it = len(np.transpose(ZZ)[tuple([np.transpose(ZZ) < ht])])
il = len(np.transpose(ZZ)[tuple([np.transpose(ZZ) > H - hl])])
os.system(
    f"octave --eval 'full_single_leak({H}, {ht}, {hl}, {L}, {ux}, "
    + f"{uy}, {Wi}, {XXf}, {YYf}, {ZZ})' "
)

# Create variables used in this script
permf = [Ka / md] * xyzf
Lg = []
t0 = []
t1 = []
injestra = []

# Define injection strategy
# [injection time [h], injection rate [m^3/s], microbial concentration [kg/m^3],
# oxygen concentration [kg/m^3], and urea concentration [kg/m^3]].
injestra.append([0.01, 2.0e-02, 0.00, 0.00, 0.0])
injestra.append([15.0, 2.0e-02, 0.00, 0.00, 60.0])
injestra.append([10.0, 2.0e-02, 0.00, 0.00, 0.0])

# We store the maximum injected oxygen and urea concentrations. They are use
# so that the computed oxygen and urea values are within these during the
# solution step
comax = 0
cumax = 0
for i in range(len(injestra)):
    comax = max(injestra[i][3], comax)
    cumax = max(injestra[i][4], cumax)

# Set and run the full and quarter domain MICP simulations
permf[xxf * yyf * it : xxf * yyf * (zz - il)] = [Kl / md] * (xxf * yyf * (zz - il - it))
var = {
    "crit": crit,
    "comax": comax,
    "cumax": cumax,
    "eta": eta,
    "F": F,
    "il": il,
    "ix": ix,
    "iy": iy,
    "injestra": injestra,
    "ka": ka,
    "kd": kd,
    "kmin": kmin,
    "ko": ko,
    "kstr": kstr,
    "ku": ku,
    "mu": mu,
    "muu": muu,
    "muw": muw,
    "perm": permf,
    "phi": phi,
    "rhob": rhob,
    "rhoc": rhoc,
    "rhow": rhow,
    "tolclg": tolclg,
    "Wc": Wc,
    "xx": xxf,
    "it": it,
    "Y": Y,
    "yy": yyf,
    "zz": zz,
}
mytemplate = Template(filename="micp.mako")
FilledTemplate = mytemplate.render(**var)
with open("inputs/directional.data", "w") as f:
    f.write(FilledTemplate)
os.system(
    f"{flow} inputs/directional.data --output-dir=results"
    + " --initial-time-step-in-days=0.0001 --solver-max-time-step-in-days=0.01"
    + " --enable-vtk-output=true & wait\n"
)

# {
# Copyright 2021-2022, NORCE Norwegian Research Centre AS, Computational
# Geosciences and Modeling.

# This file is part of the py-micp module.

# py-micp is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ad-wa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
# }
