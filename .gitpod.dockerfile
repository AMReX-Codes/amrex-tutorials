FROM gitpod/workspace-full-vnc

RUN sudo apt-get update  && sudo apt-get install -y --no-install-recommends build-essential g++ libopenmpi-dev openmpi-bin paraview python3-paraview ffmpeg python3 python3-pip python3-setuptools && sudo rm -rf /var/lib/apt/lists/*