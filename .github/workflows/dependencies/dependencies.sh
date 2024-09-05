#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends\
    build-essential \
    g++ gfortran    \
    libopenmpi-dev  \
    openmpi-bin     \
    python3         \
    python3-pip

python3 -m pip install -U pip
python3 -m pip install -U build packaging setuptools wheel
