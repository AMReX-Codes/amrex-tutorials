#!/usr/bin/env bash

set -eu -o pipefail

sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 2
sudo update-alternatives --set python /usr/bin/python3
python -m pip install --upgrade pip
python -m pip install --upgrade wheel
python -m pip install --upgrade cmake matplotlib==3.2.2 numpy yt
