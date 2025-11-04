#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

brew update
brew install gcc@15 || true
brew install libomp || true
brew install --cc=gcc-15 open-mpi || true
brew install ccache || true

