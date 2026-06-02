# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: GPL-3.0
# pylint: disable=C0116

"""Command flag for the path of the MRST 'startup.m'"""

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--mrst_startup",
        action="store",
        default="/Users/dmar/Github/py-micp/MRST/startup.m",
    )


@pytest.fixture(scope="session")
def mrst_startup(request):
    return request.config.option.mrst_startup
