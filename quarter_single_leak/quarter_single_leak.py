# Setting up and simulating a 3D flow system with a straight leakage path, where
# CO2 is assessed before and after MICP treatment. Information about the model
# and parameters can be found in:

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
import meshio
import numpy as np
from mako.template import Template
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

# Set the model parameters
rhob = 35  # Density (biofilm) [kg/m^3]
rhoc = 2710  # Density (calcite) [kg/m^3]
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
day = 86400  # Seconds in a day [s]
QCO2 = 1e-4  # Injection rate (CO2) [m^3/s]
dtCO2 = 10  # Time step to print CO2 results [d]
NCO2 = 40  # Number of time steps to run the simulation [-]
detexp = 0.58  # Exponent for the norm of the detahcment rate [-]
yurca = 1.67  # Yielc urea to calcite [-]

# Delete previous decks and results
os.system("rm -rf decks & rm -rf results & wait")
os.system("mkdir decks & mkdir results & wait")

# Define the porosity-permeability functions (aquifer and leakage)
# OLD formulation (release v2024.10)
# def perma(bio, cal):
#     return (
#         (Ka * ((phi - bio - cal - crit) / (phi - crit)) ** eta + kmin)
#         * Ka
#         / (md * (Ka + kmin))
#     )
# def perml(bio, cal):
#     return (
#         (Kl * ((phi - bio - cal - crit) / (phi - crit)) ** eta + kmin)
#         * Kl
#         / (md * (Kl + kmin))
#     )

# NOW (after 2025.04) we use the PERMFACT table, see the OPM Flow manual
porofac = [0, 0.4, 0.6, 0.9, 1]
permfac = [0.005, 0.005, 0.01, 0.1, 1]


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


# Create the grid
H = 30  # Height of the domain [m]
ht = 5  # Height of the top aquifer [m]
hl = 5  # Height of the lower aquifer [m]
L = 100  # Length of the domain [m]
ux = 13  # X-gap between aquifer and leakage [m]
uy = 14  # Y-gap between aquifer and leakage [m]
Wi = 100  # Width of the domain [m]
X0 = np.linspace(0, 21, 22)  # Discretization of the x-dir X0
X1 = np.linspace(-1.5, 0, 16)  # ""                          L*exp(X1)
Y0 = np.linspace(0, 21, 22)  # Discretization of the y-dir Y0
Y1 = np.linspace(-1.5, 0, 16)  # ""                          Wi*exp(Y1)
Z0 = np.linspace(0, ht, ht + 1)  # Discretization of the z-dir Z0
Z1 = np.linspace(ht + 1.0, H - hl + 1.0, H - hl - ht + 1)  # ""     Z1
Z2 = np.linspace(H - hl + 2.0, H, 4)  # ""                          Z2
Wc = 1  # Cells for each MICP well injection
XX, YY, ZZ = [], [], []
for i in range(len(X0)):
    XX.append(float(X0[i]))
for i in range(len(X1)):
    XX.append(float(L * np.exp(X1[i])))
for j in range(len(Y0)):
    YY.append(float(Y0[j]))
for j in range(len(Y1)):
    YY.append(float(Wi * np.exp(Y1[j])))
for k in range(len(Z0)):
    ZZ.append(float(Z0[k]))
for k in range(len(Z1)):
    ZZ.append(float(Z1[k]))
for k in range(len(Z2)):
    ZZ.append(float(Z2[k]))
xx = len(X0) + len(X1) - 1
yy = len(Y0) + len(Y1) - 1
zz = len(Z0) + len(Z1) + len(Z2) - 1
xyz = xx * yy * zz
ix = np.abs(np.transpose(XX)).argmin() + 1
iy = np.abs(np.transpose(YY)).argmin() + 1
it = len(np.transpose(ZZ)[tuple([np.transpose(ZZ) < ht])])
il = len(np.transpose(ZZ)[tuple([np.transpose(ZZ) > H - hl])])
os.system(
    f"octave --eval 'quarter_single_leak({H}, {ht}, {hl}, {L}, {ux}, "
    + f"{uy}, {Wi}, {XX}, {YY}, {ZZ})' "
)

# Create variables used in this script
poro = phi * np.ones(xyz, dtype=float)
biof = [0] * xyz
calc = [0] * xyz
perm = [Ka / md] * xyz
injestra, Lg, t0, t1, ACTNUM, VOLUM = [], [], [], [], [], []

# Set and run the CO2 simulation before MICP treatment
perm[xx * yy * it : xx * yy * (zz - il)] = [Kl / md] * (xx * yy * (zz - il - it))
permw = compact_format("".join(f"{val} " for val in perm).split())
porow = compact_format("".join(f"{val} " for val in poro).split())
mytemplate = Template(filename="co2.mako")
var = {
    "dtCO2": dtCO2,
    "il": il,
    "it": it,
    "NCO2": NCO2,
    "perm": permw,
    "poro": porow,
    "QCO2": QCO2,
    "xx": xx,
    "yy": yy,
    "zz": zz,
}
FilledTemplate = mytemplate.render(**var)
with open("decks/CO2.DATA", "w") as f:
    f.write(FilledTemplate)
os.system(f"{flow} decks/CO2.DATA")

# Define injection strategy
# [injection time [h], injection rate [m^3/s], microbial concentration [kg/m^3],
# oxygen concentration [kg/m^3], and urea concentration [kg/m^3]]. Before each
# injection of the solutions in the phases we inject only water for a very short
# time [0.01 hour] to ease the convergence of the solution. These are the
# optimized times as described in the paper.
# NOW (after 2025.04) we multiply the injection length by two since by updating
# the decks with more realistic pressures, the wells are change to BHP control
injestra.append([2 * 14.9298732684, 2.0e-02, 0.01, 0, 0])
injestra.append([2 * 0.4180189097, 2.0e-02, 0, 0, 0])
injestra.append([2 * 2.5208976414, 0, 0, 0, 0])
injestra.append([0.0100000000, 2.0e-02, 0, 0, 0])
injestra.append([2 * 15.3602915202, 2.0e-02, 0, 0.04, 0])
injestra.append([2 * 1.7326877982, 2.0e-02, 0, 0, 0])
injestra.append([2 * 1.4059078334, 0, 0, 0, 0])
injestra.append([0.0100000000, 2.0e-02, 0, 0, 0])
injestra.append([2 * 0.3929312540, 2.0e-02, 0, 0, 60])
injestra.append([2 * 1.1262173564, 2.0e-02, 0, 0, 0])
injestra.append([2 * 6.0153179270, 0, 0, 0, 0])
injestra.append([0.0100000000, 2.0e-02, 0, 0, 0])
injestra.append([2 * 111.1739194383, 2.0e-02, 0, 0.04, 60])
injestra.append([2 * 0.9334330084, 2.0e-02, 0, 0, 0])
injestra.append([2 * 1.8380402576, 0, 0, 0, 0])
injestra.append([0.0100000000, 2.0e-02, 0, 0, 0])
injestra.append([2 * 1.2220241699, 2.0e-02, 0, 0.04, 60])
injestra.append([2 * 1.1649655408, 2.0e-02, 0, 0, 0])
injestra.append([2 * 1.2496253760, 0, 0, 0, 0])

# Read the input file
mytemplate = Template(filename="micp.mako")

# Set and run the MICP simulation
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
    "xx": xx,
    "Y": Y,
    "yy": yy,
    "zz": zz,
}
FilledTemplate = mytemplate.render(**var)
with open("decks/MICP.DATA", "w") as f:
    f.write(FilledTemplate)
os.system(f"{flow} decks/MICP.DATA")

# Update the porosity and permeability after MICP treatment
permred = interp1d(porofac, permfac, fill_value="extrapolate")
with open("decks/GRID.INC", "r") as f:
    list_of_lines = f.readlines()
act = list_of_lines[-int(xyz / 20.0) - 1 :]
for k in range(len(act)):
    for i in range(len(act[k])):
        if act[k][i] == "0" or act[k][i] == "1":
            ACTNUM.append(int(act[k][i]))
with open("results/MICP.pvd", "r") as f:
    list_of_lines = f.readlines()
ss = list_of_lines[-3]
mesh = meshio.read(f"results/MICP-{int(ss[-13:-8]):05d}.vtu")
for row in mesh.cell_data["calcite volume fraction"]:
    c = row.flatten()
for (
    key
) in (
    mesh.cell_data.keys()
):  # Flow 2025.10 prints biofilm voulme, this has been fixed in Flow master
    if "biofilm" in key:
        biokey = key
        break
for row in mesh.cell_data[biokey]:
    b = row.flatten()
biof[: xx * yy * it] = b[: xx * yy * it]
calc[: xx * yy * it] = c[: xx * yy * it]
perm[: xx * yy * it] = Ka * permred(phi - b[: xx * yy * it] - c[: xx * yy * it])
j = 0
for k in range(it, zz - il):
    for l in range(0, xx * yy):
        if ACTNUM[xx * yy * k + l] == 1:
            biof[xx * yy * k + l] = b[xx * yy * it + j]
            calc[xx * yy * k + l] = c[xx * yy * it + j]
            perm[xx * yy * k + l] = Kl * permred(
                phi - b[xx * yy * it + j] - c[xx * yy * it + j]
            )
            j = j + 1
i = xx * yy * zz - xx * yy * il
I = len(b) - xx * yy * il
biof[i:] = b[I:]
calc[i:] = c[I:]
perm[i:] = Ka * permred(phi - b[I:] - c[I:])

# Set and run the CO2 simulation after MICP treatment
permw = compact_format("".join(f"{float(val)} " for val in perm).split())
poro = [por - bio - cal for por, bio, cal in zip(poro, biof, calc)]
porow = compact_format("".join(f"{float(val)} " for val in poro).split())
mytemplate = Template(filename="co2.mako")
var = {
    "dtCO2": dtCO2,
    "il": il,
    "it": it,
    "NCO2": NCO2,
    "perm": permw,
    "poro": porow,
    "QCO2": QCO2,
    "xx": xx,
    "yy": yy,
    "zz": zz,
}
FilledTemplate = mytemplate.render(**var)
with open("decks/CO2MICP.DATA", "w") as f:
    f.write(FilledTemplate)
os.system(f"{flow} decks/CO2MICP.DATA")

# Obtain the times for the .vtu results and save them in 't' in days
with open("results/CO2.pvd", "r") as f:
    list_of_lines = f.readlines()
nn, k, T = 0, 0, []
for i in range(6, len(list_of_lines)):
    a = []
    for j in range(22, len(list_of_lines[i])):
        a.append(list_of_lines[i][j])
        if list_of_lines[i][j + 1] == '"':
            T.append(a)
            break
for k in range(len(T)):
    aaa = str(T[k]).lstrip("[").rstrip("]")
    bbb = aaa.replace(", ", "")
    xxx = 0
    n = 0
    for i in range(len(bbb)):
        if bbb[i] != "'":
            if bbb[len(bbb) - i - 1] == ".":
                xxx = xxx * pow(10, -n)
                n = 0
            else:
                xxx = xxx + pow(10, n) * float(bbb[len(bbb) - i - 1])
                n = n + 1
    t0.append(xxx / day)
with open("results/CO2MICP.pvd", "r") as f:
    list_of_lines = f.readlines()
nn, k, T = 0, 0, []
for i in range(6, len(list_of_lines)):
    a = []
    for j in range(22, len(list_of_lines[i])):
        a.append(list_of_lines[i][j])
        if list_of_lines[i][j + 1] == '"':
            T.append(a)
            break
for k in range(len(T)):
    aaa = str(T[k]).lstrip("[").rstrip("]")
    bbb = aaa.replace(", ", "")
    xxx, n = 0, 0
    for i in range(len(bbb)):
        if bbb[i] != "'":
            if bbb[len(bbb) - i - 1] == ".":
                xxx = xxx * pow(10, -n)
                n = 0
            else:
                xxx = xxx + pow(10, n) * float(bbb[len(bbb) - i - 1])
                n = n + 1
    t1.append(xxx / day)

# Compute the mass percentage of leaked CO2 to the upper aquifer before and after
# MICP treatment
for k in range(it):
    z = ZZ[k + 1] - ZZ[k]
    for j in range(yy):
        y = YY[j + 1] - YY[j]
        for i in range(xx):
            VOLUM.append([(XX[i + 1] - XX[i]) * y * z])
a = []
II = xx * yy * it - 1
for j in range(len(t0)):
    mesh = meshio.read(f"results/CO2-{j:05d}.vtu")
    for row in mesh.cell_data["saturation_gas"]:
        CO2 = row
    for row in mesh.cell_data["porosity"]:
        p = row
    CO2M = QCO2 * t0[j] * day
    if t0[j] > 0:
        a.append(100.0 * np.sum(CO2[:II] * p[:II] * VOLUM[:II]) / CO2M)
    else:
        a.append(np.sum(CO2[:II] * p[:II] * VOLUM[:II]))
Lg.append(a)
a = []
for j in range(len(t1)):
    mesh = meshio.read(f"results/CO2MICP-{j:05d}.vtu")
    for row in mesh.cell_data["saturation_gas"]:
        CO2 = row
    for row in mesh.cell_data["porosity"]:
        p = row
    CO2M = QCO2 * t1[j] * day
    if t1[j] > 0:
        a.append(100.0 * np.sum(CO2[:II] * p[:II] * VOLUM[:II]) / CO2M)
    else:
        a.append(np.sum(CO2[:II] * p[:II] * VOLUM[:II]))
Lg.append(a)
print(f"Percentage of leake CO2:{Lg[0][-1]}")
print(f"Percentage of leake CO2 after MICP treatment:{Lg[1][-1]}")

# Plot the results over time and save them as co2mass_comparison.png
lw = 5
plt.figure(figsize=(5, 4), dpi=256)
plt.rc("font", size=9)
plt.rc("legend", fontsize=9)
axes = plt.subplot(1, 1, 1)
plt.plot(
    t0, Lg[0][:], color=[1, 0.2, 0.2], linewidth=lw, linestyle="-", label="Without MICP"
)
plt.plot(
    t1, Lg[1][:], color=[1, 0.5, 0], linewidth=lw, linestyle="-", label="After MICP"
)
plt.xlabel("t [days]")
plt.ylabel("Leaked $CO_{2}$/injected $CO_{2}$ [%]")
plt.grid()
plt.legend(loc="upper left")
plt.savefig("results/co2mass_comparison.png")
plt.show()

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
