#!/usr/bin/env python
# SPDX-FileCopyrightText: 2021-2026 NORCE Research AS
# SPDX-License-Identifier: GPL-3.0
# pylint: disable=C0116,R0913,R0914,R0915,R0917,R0801

"""Setting up and simulating a 3D flow system with a straight leakage path, where
CO2 is assessed before and after MICP treatment. Information about the model
and parameters can be found in:

Landa-Marbán, D., Kumar, K., Tveit, S., Gasda, S.E. Numerical studies of CO2
leakage remediation by micp-based plugging technology. In: Røkke, N.A. and
Knuutila, h.K. (Eds) Short Papers from the 11th International Trondheim CCS
conference, ISBN: 978-82-536-1714-5, 284-290"""

import shutil
import subprocess
import xml.etree.ElementTree as element_tree
from pathlib import Path
import meshio
import numpy as np
from mako.template import Template
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

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


def read_pvd_entries(pvd_path, day):
    tree = element_tree.parse(pvd_path)
    root = tree.getroot()
    dataset_paths, times = [], []
    base_dir = Path(pvd_path).parent
    for dataset in root.iter():
        if not dataset.tag.endswith("DataSet"):
            continue
        dataset_file = dataset.attrib.get("file")
        dataset_time = dataset.attrib.get("timestep")
        if dataset_file is None or dataset_time is None:
            continue
        dataset_paths.append(str(base_dir / dataset_file))
        times.append(float(dataset_time) / day)
    return dataset_paths, times


def parse_keyword_values(file_path, keyword):
    with open(file_path, "r", encoding="utf-8") as file_handle:
        lines = file_handle.readlines()
    values, reading = [], False
    for line in lines:
        stripped = line.strip()
        if not reading:
            if stripped.startswith(keyword):
                trailing = stripped[len(keyword) :].strip()
                line_tokens = trailing.split() if trailing else []
                reading = True
            else:
                continue
        else:
            line_tokens = stripped.split()
        for line_token in line_tokens:
            token = line_token.rstrip("/")
            if token:
                if "*" in token:
                    count_text, value_text = token.split("*", 1)
                    try:
                        count = int(count_text)
                        value = int(float(value_text))
                        values.extend([value] * count)
                    except ValueError:
                        pass
                else:
                    try:
                        values.append(int(float(token)))
                    except ValueError:
                        pass
            if line_token.endswith("/"):
                return np.asarray(values, dtype=int)
    return np.asarray(values, dtype=int)


def flatten_cell_data(mesh, key):
    data = mesh.cell_data[key]
    if len(data) == 1:
        return np.asarray(data[0]).ravel()
    return np.concatenate([np.asarray(row).ravel() for row in data])


def compact_format(v):
    change_idx = np.flatnonzero(np.diff(v, prepend=v[0] - 1))
    counts = np.diff(np.append(change_idx, v.size))
    vals = v[change_idx]

    out = []
    for val, n in zip(vals, counts):
        val_str = str(int(val)) if float(val).is_integer() else str(val)
        out.append(f"{n}*{val_str} " if n > 1 else f"{val_str} ")
    return out


def compute_leakage_series(dataset_paths, times, qco2, day, volume, top_cells):
    a = []
    active_volumes = volume[:top_cells]
    for dataset_index, dataset_path in enumerate(dataset_paths):
        mesh = meshio.read(dataset_path)
        co2 = flatten_cell_data(mesh, "saturation_gas")
        p = flatten_cell_data(mesh, "porosity")
        leaked_mass = np.sum(co2[:top_cells] * p[:top_cells] * active_volumes)
        co2m = qco2 * times[dataset_index] * day
        if times[dataset_index] > 0:
            a.append(100.0 * leaked_mass / co2m)
        else:
            a.append(leaked_mass)
    return a


def full_single_leak():
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
    day = 86400  # Seconds in a day [s]
    qco2 = 4e-4  # Injection rate (CO2) [m^3/s]
    dtco2 = 10  # Time step to print CO2 results [d]
    nco2 = 40  # Number of time steps to run the simulation [-]
    detexp = 0.58  # Exponent for the norm of the detahcment rate [-]
    yurca = 1.67  # Yielc urea to calcite [-]

    for required_file in ("co2.mako", "micp.mako"):
        if not Path(required_file).is_file():
            raise FileNotFoundError(f"Required file not found: {required_file}")

    # Delete previous inputs and results
    ensure_clean_directory("decks")
    ensure_clean_directory("results")

    # Define the porosity-permeability functions (aquifer and leakage)
    # OLD formulation (release v2024.10)
    # def perma(bio, cal):
    #     return (
    #         (k * ((phi - bio - cal - crit) / (phi - crit)) ** eta + kmin)
    #         * k
    #         / (md * (k + kmin))
    #     )
    # def perml(bio, cal):
    #     return (
    #         (kl * ((phi - bio - cal - crit) / (phi - crit)) ** eta + kmin)
    #         * kl
    #         / (md * (kl + kmin))
    #     )

    # NOW (after 2025.04) we use the PERMFACT table, see the OPM Flow manual
    porofac = [0, 0.4, 0.6, 0.9, 1]
    permfac = [0.005, 0.005, 0.01, 0.1, 1]

    # Create the grid
    h = 30  # Height of the domain [m]
    ht = 5  # Height of the top aquifer [m]
    hl = 5  # Height of the lower aquifer [m]
    le = 100  # Length of the domain [m]
    ux = 13  # x-gap between aquifer and leakage [m]
    uy = 14  # y-gap between aquifer and leakage [m]
    wi = 100  # Width of the domain [m]
    x0 = np.linspace(0, -1.5, 16)  # ""                         -le*exp(x0)
    x1 = np.linspace(-21, 21, 44)  # Discretization of the x-dir x1
    x2 = np.linspace(-1.5, 0, 16)  # ""                          le*exp(x2)
    y0 = np.linspace(0, -1.5, 16)  # ""                         -wi*exp(y0)
    y1 = np.linspace(-21, 21, 44)  # Discretization of the y-dir y1
    y2 = np.linspace(-1.5, 0, 16)  # ""                          wi*exp(y2)
    z0 = np.linspace(0, ht, ht + 1)  # Discretization of the z-dir z0
    z1 = np.linspace(ht + 1.0, h - hl + 1.0, h - hl - ht + 1)  # ""     z1
    z2 = np.linspace(h - hl + 2.0, h, 4)  # ""                          z2
    wc = 1  # Cells for each MICP well injection
    x = np.concatenate((-le * np.exp(x0), x1, le * np.exp(x2))).astype(float)
    y = np.concatenate((-wi * np.exp(y0), y1, wi * np.exp(y2))).astype(float)
    z = np.concatenate((z0, z1, z2)).astype(float)
    xx = len(x0) + len(x1) + len(x2) - 1
    yy = len(y0) + len(y1) + len(y2) - 1
    zz = len(z0) + len(z1) + len(z2) - 1
    xyz = xx * yy * zz
    ix = int(np.abs(x).argmin()) + 1
    iy = int(np.abs(y).argmin()) + 1
    it = int(np.count_nonzero(z < ht))
    il = int(np.count_nonzero(z > h - hl))
    run_command(
        f"octave --eval 'full_single_leak({h}, {ht}, {hl}, {le}, {ux}, "
        f"{uy}, {wi}, {x.tolist()}, {y.tolist()}, {z.tolist()})'"
    )

    # Create variables used in this script
    poro = np.full(xyz, phi, dtype=float)
    biof = np.zeros(xyz, dtype=float)
    calc = np.zeros(xyz, dtype=float)
    injestra = []
    perm = np.full(xyz, k / md, dtype=float)
    lg = []
    t0 = []
    t1 = []
    actnum = np.asarray([], dtype=int)
    volume = np.asarray([], dtype=float)

    # Set and run the CO2 simulation before MICP treatment
    perm[xx * yy * it : xx * yy * (zz - il)] = kl / md
    permw = compact_format(perm)
    porow = compact_format(poro)
    var = {
        "biof": biof.tolist(),
        "calc": calc.tolist(),
        "dtco2": dtco2,
        "il": il,
        "ix": ix,
        "it": it,
        "iy": iy,
        "nco2": nco2,
        "perm": permw,
        "poro": porow,
        "qco2": qco2,
        "xx": xx,
        "yy": yy,
        "zz": zz,
    }
    render_template_to_file("co2.mako", "decks/CO2.DATA", var)
    run_command(f"{flow} decks/CO2.DATA")

    # Define injection strategy
    # [injection time [h], injection rate [m^3/s], microbial concentration [kg/m^3],
    # oxygen concentration [kg/m^3], and urea concentration [kg/m^3]]. Before each
    # injection of the solutions in the phases we inject only water for a very short
    # time [0.01 hour] to ease the convergence of the solution.
    # This injection strategy results in the full remediation of the leakage path
    # (i.e., the percentage of leake CO2 is zero) and it is identical to the ad-hoc
    # one in the quarter_single_leak example (the injection rate is scaled up by a
    # factor of four.
    # NOW (after 2025.04) this strategy seems not good, the task is to find one
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
    injestra.append([0.01, 8e-2, 0, 0, 0])
    injestra.append([30.0, 8e-2, 0, 0.04, 60.0])
    injestra.append([2.0, 8e-2, 0, 0, 0])
    injestra.append([30.0, 0, 0, 0, 0])
    injestra.append([0.01, 8e-2, 0, 0, 0])
    injestra.append([30.0, 8e-2, 0, 0.04, 60.0])
    injestra.append([2.0, 8e-2, 0, 0, 0])
    injestra.append([30.0, 0, 0, 0, 0])
    injestra.append([0.01, 8e-2, 0, 0, 0])
    injestra.append([30.0, 8e-2, 0, 0.04, 60.0])
    injestra.append([2.0, 8e-2, 0, 0, 0])
    injestra.append([30.0, 0, 0, 0, 0])

    # Read the input file

    # Set and run the MICP simulation
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
        "xx": xx,
        "yld": yld,
        "yy": yy,
        "zz": zz,
    }
    render_template_to_file("micp.mako", "decks/MICP.DATA", var)
    run_command(f"{flow} decks/MICP.DATA")

    # Update the porosity and permeability after MICP treatment
    permred = interp1d(porofac, permfac, fill_value="extrapolate")
    actnum = parse_keyword_values("decks/GRID.INC", "ACTNUM")
    if len(actnum) < xyz:
        raise ValueError(
            f"ACTNUM contains {len(actnum)} values, expected at least {xyz}"
        )
    micp_dataset_paths, _ = read_pvd_entries("results/MICP.pvd", day)
    if not micp_dataset_paths:
        raise ValueError("No VTU datasets found in results/MICP.pvd")
    mesh = meshio.read(micp_dataset_paths[-1])
    c = flatten_cell_data(mesh, "calcite volume fraction")
    biokey = next((key for key in mesh.cell_data.keys() if "biofilm" in key), None)
    if biokey is None:
        raise KeyError("No biofilm cell_data entry found in MICP VTU file")
    b = flatten_cell_data(mesh, biokey)
    top_cells = xx * yy * it
    mid_start = top_cells
    mid_stop = xx * yy * (zz - il)
    bottom_count = xx * yy * il
    bottom_start = xyz - bottom_count
    biof[:top_cells] = b[:top_cells]
    calc[:top_cells] = c[:top_cells]
    perm[:top_cells] = k * permred(phi - b[:top_cells] - c[:top_cells])
    active_mask = actnum[mid_start:mid_stop] == 1
    active_indices = np.flatnonzero(active_mask) + mid_start
    active_count = active_indices.size
    source_slice = slice(mid_start, mid_start + active_count)
    biof[active_indices] = b[source_slice]
    calc[active_indices] = c[source_slice]
    perm[active_indices] = kl * permred(phi - b[source_slice] - c[source_slice])
    biof[bottom_start:] = b[-bottom_count:]
    calc[bottom_start:] = c[-bottom_count:]
    perm[bottom_start:] = k * permred(phi - b[-bottom_count:] - c[-bottom_count:])

    # Set and run the CO2 simulation after MICP treatment
    permw = compact_format(perm)
    poro = poro - biof - calc
    porow = compact_format(poro)
    var = {
        "biof": biof.tolist(),
        "calc": calc.tolist(),
        "dtco2": dtco2,
        "il": il,
        "ix": ix,
        "it": it,
        "iy": iy,
        "nco2": nco2,
        "perm": permw,
        "poro": porow,
        "qco2": qco2,
        "xx": xx,
        "yy": yy,
        "zz": zz,
    }
    render_template_to_file("co2.mako", "decks/CO2MICP.DATA", var)
    run_command(f"{flow} decks/CO2MICP.DATA")

    # Obtain the times for the .vtu results and save them in 't' in days
    co2_dataset_paths, t0 = read_pvd_entries("results/CO2.pvd", day)
    co2micp_dataset_paths, t1 = read_pvd_entries("results/CO2MICP.pvd", day)

    # Compute the mass percentage of leaked CO2 to the upper aquifer before and after
    # MICP treatment
    x_widths = np.diff(x[: xx + 1])
    y_widths = np.diff(y[: yy + 1])
    z_widths = np.diff(z[: it + 1])
    volume = (
        z_widths[:, None, None] * y_widths[None, :, None] * x_widths[None, None, :]
    ).reshape(-1)
    lg.append(
        compute_leakage_series(co2_dataset_paths, t0, qco2, day, volume, top_cells)
    )
    lg.append(
        compute_leakage_series(co2micp_dataset_paths, t1, qco2, day, volume, top_cells)
    )
    print(f"Percentage of leake CO2:{lg[0][-1]}")
    print(f"Percentage of leake CO2 after MICP treatment:{lg[1][-1]}")

    # Plot the results over time and save them as co2mass_comparison.png
    lw = 5
    plt.figure(figsize=(5, 4), dpi=256)
    plt.rc("font", size=9)
    plt.rc("legend", fontsize=9)
    plt.plot(
        t0,
        lg[0][:],
        color=[1, 0.2, 0.2],
        linewidth=lw,
        linestyle="-",
        label="Without MICP",
    )
    plt.plot(
        t1, lg[1][:], color=[1, 0.5, 0], linewidth=lw, linestyle="-", label="After MICP"
    )
    plt.xlabel("t [days]")
    plt.ylabel("Leaked $CO_{2}$/injected $CO_{2}$ [%]")
    plt.grid()
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig("results/co2mass_comparison.png")
    plt.show()


if __name__ == "__main__":
    full_single_leak()
