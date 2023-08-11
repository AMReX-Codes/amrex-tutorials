#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y  \
    build-essential      \
    clang gfortran       \
    python3              \
    python3-pip

python3 -m pip install -U pip setuptools wheel
