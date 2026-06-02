# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: GPL-3.0
# pylint: disable=R0801,R0914

"""Test the quarter_double_leak model"""

import pathlib
import shutil
import subprocess
from pathlib import Path
import numpy as np
from opm.io.ecl import EGrid as OpmGrid
from opm.io.ecl import EclFile as OpmFile

EPS_PORV = 1e-1
EPS = 1e-6

mainpth = pathlib.Path(__file__).parents[1]


def test_quarter_double_leak(tmp_path, monkeypatch, mrst_startup):
    """See quarter_double_leak"""
    monkeypatch.chdir(tmp_path)
    shutil.copytree(mainpth / "quarter_double_leak", tmp_path, dirs_exist_ok=True)
    for value in [
        "nco2 = 40|nco2 = 1",
        "injestra,|[[1,2e-2,0.01,0,0]],",
        "plt.sh|#plt.sh",
    ]:
        subprocess.run(
            ["sed", "-i.bak", f"s|{value}|g", "quarter_double_leak.py"],
            check=True,
        )
    Path("quarter_double_leak.py.bak").unlink(missing_ok=True)
    default = "/Users/dmar/Github/py-micp/MRST/startup.m"
    subprocess.run(
        ["sed", "-i.bak", f"s|{default}|{mrst_startup}|g", "quarter_double_leak.m"],
        check=True,
    )
    Path("quarter_double_leak.m.bak").unlink(missing_ok=True)

    subprocess.run(
        [
            "python3",
            "quarter_double_leak.py",
        ],
        check=True,
    )

    for case, inj, pro in zip(["CO2", "MICP", "CO2MICP"], [3, 16, 3], [9, 3, 9]):
        name = tmp_path / "decks" / f"{case}.DATA"
        with open(name, "r", encoding="utf8") as f:
            lines = f.readlines()
        content = "".join(lines)
        assert content.count("TSTEP") == 1
        assert content.count("INJE0") == inj
        assert content.count("PROD0") == pro
        egrid = OpmGrid(f"{tmp_path}/results/{case.upper()}.EGRID")
        nx, ny, nz = egrid.dimension
        assert nx == 37
        assert ny == 37
        assert nz == 30
        init = OpmFile(f"{tmp_path}/results/{case.upper()}.INIT")
        assert abs(np.sum(init["PORV"]) - 15006.0) < EPS_PORV
        assert abs(min(init["DX"]) - 1) < EPS
        assert abs(max(init["DX"]) - 9.516258) < EPS
        assert abs(min(init["DY"]) - 1) < EPS
        assert abs(max(init["DY"]) - 9.516258) < EPS
        assert abs(min(init["DZ"]) - 1) < EPS
        assert abs(max(init["DZ"]) - 1) < EPS
        name = tmp_path / "results" / f"{case}-00001.vtu"
        with open(name, "r", encoding="utf8") as f:
            lines = f.readlines()
        assert len(lines) == 47951 if case == "MICP" else 42216

    assert (tmp_path / "results" / "co2mass_comparison.png").exists()
