# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: GPL-3.0
# pylint: disable=R0801,R0914

"""Test the towards_more_optimal_injection_strategies model"""

import pathlib
import shutil
import subprocess
from pathlib import Path

EPS = 1e-6

mainpth = pathlib.Path(__file__).parents[1]


def test_towards_more_optimal_injection_strategies(tmp_path, monkeypatch, mrst_startup):
    """See towards_more_optimal_injection_strategies"""
    monkeypatch.chdir(tmp_path)
    shutil.copytree(
        mainpth / "towards_more_optimal_injection_strategies",
        tmp_path,
        dirs_exist_ok=True,
    )
    subprocess.run(
        [
            "sed",
            "-i.bak",
            "s|injestra,|[[1,2e-2,0.01,0,0]],|g",
            "towards_more_optimal_injection_strategies.py",
        ],
        check=True,
    )
    Path("towards_more_optimal_injection_strategies.py.bak").unlink(missing_ok=True)
    default = "/Users/dmar/Github/py-micp/MRST/startup.m"
    subprocess.run(
        ["sed", "-i.bak", f"s|{default}|{mrst_startup}|g", "full_single_leak.m"],
        check=True,
    )
    Path("full_single_leak.m.bak").unlink(missing_ok=True)

    subprocess.run(
        [
            "python3",
            "towards_more_optimal_injection_strategies.py",
        ],
        check=True,
    )

    name = tmp_path / "decks" / "DIRECTIONAL.DATA"
    with open(name, "r", encoding="utf8") as f:
        lines = f.readlines()
    content = "".join(lines)
    assert content.count("INJE0") == 28
    assert content.count("PROD0") == 12
    name = tmp_path / "results" / "DIRECTIONAL-00001.vtu"
    with open(name, "r", encoding="utf8") as f:
        lines = f.readlines()
    assert len(lines) == 195638
