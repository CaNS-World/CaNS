#!/bin/sh
#
# installs NVHPC SDK for building with the NVIDIA compiler
# see: https://developer.nvidia.com/hpc-sdk-downloads
#
# get packages
#
sudo apt-get install libfftw3-dev
#
NVHPC_VERSION_A=25.3
NVHPC_VERSION_B="$(echo $NVHPC_VERSION_A | sed 's/\./-/g')"
curl https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list
sudo apt-get update -y
sudo apt-get install -y nvhpc-${NVHPC_VERSION_B}
#sudo apt-get install -y nvhpc-${NVHPC_VERSION_B}-cuda-multi
#
# save NVHPC version in GITHUB_ENV for later use when setting the NVHPC environment
#
echo "NVHPC_VERSION_A=$NVHPC_VERSION_A" >> $GITHUB_ENV
echo "NVHPC_VERSION_B=$NVHPC_VERSION_B" >> $GITHUB_ENV
