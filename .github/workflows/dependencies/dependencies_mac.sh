#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

brew update
brew reinstall gcc || true
brew install libomp || true
brew install open-mpi || true
brew install ccache || true
