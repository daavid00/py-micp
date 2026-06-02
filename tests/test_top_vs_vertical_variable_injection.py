# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: GPL-3.0
# pylint: disable=R0801,R0914

"""Test the top_vs_vertical_variable_injection model"""

import pathlib
import shutil
import subprocess
from pathlib import Path

EPS = 1e-6

mainpth = pathlib.Path(__file__).parents[1]


def test_top_vs_vertical_variable_injection(tmp_path, monkeypatch, mrst_startup):
    """See top_vs_vertical_variable_injection"""
    monkeypatch.chdir(tmp_path)
    shutil.copytree(
        mainpth / "top_vs_vertical_variable_injection", tmp_path, dirs_exist_ok=True
    )
    subprocess.run(
        [
            "sed",
            "-i.bak",
            "s|injestra,|[[1,2e-2,0.01,0.01,0.01]],|g",
            "top_vs_vertical_variable_injection.py",
        ],
        check=True,
    )
    Path("top_vs_vertical_variable_injection.py.bak").unlink(missing_ok=True)
    default = "/Users/dmar/Github/py-micp/MRST/startup.m"
    subprocess.run(
        ["sed", "-i.bak", f"s|{default}|{mrst_startup}|g", "full_single_leak.m"],
        check=True,
    )
    Path("full_single_leak.m.bak").unlink(missing_ok=True)

    subprocess.run(
        [
            "python3",
            "top_vs_vertical_variable_injection.py",
        ],
        check=True,
    )

    for case, inj in zip(["TOP", "VARIABLE"], [0, 0.01]):
        name = tmp_path / "decks" / f"{case}.DATA"
        with open(name, "r", encoding="utf8") as f:
            lines = f.readlines()
        content = "".join(lines)
        assert f"INJE03 0 {inj}" in content
        name = tmp_path / "results" / f"{case}-00001.vtu"
        with open(name, "r", encoding="utf8") as f:
            lines = f.readlines()
        assert len(lines) == 195638
