#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update

sudo apt-get install -y --no-install-recommends \
    build-essential    \
    g++-10 gfortran-10 \
    libopenmpi-dev     \
    openmpi-bin        \
    python3            \
    python3-pip

python3 -m pip install -U pip setuptools wheel
