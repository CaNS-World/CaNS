#!/bin/sh
#
# installs NVHPC SDK for building with the NVIDIA compiler
# see: https://developer.nvidia.com/hpc-sdk-downloads
#
# optional environment variables:
#   NVHPC_VERSION=latest|26.3
#
# get packages
#
set -eu
#
sudo apt-get install -y --no-install-recommends libfftw3-dev
#
NVHPC_VERSION_A="${NVHPC_VERSION_A:-${NVHPC_VERSION:-latest}}"
#
curl -fsSL https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list
sudo apt-get update -y -o Acquire::Retries=3
#
if [ "$NVHPC_VERSION_A" = "latest" ]; then
	sudo apt-get install -y --no-install-recommends nvhpc
	#sudo apt-get install -y --no-install-recommends nvhpc-cuda-multi
	NVHPC_VERSION_A="$(dpkg-query -W -f='${Version}\n' nvhpc | sed 's/-0$//')"
else
	sudo apt-get install -y --no-install-recommends "nvhpc-$(echo "$NVHPC_VERSION_A" | sed 's/\./-/g')"
	#sudo apt-get install -y --no-install-recommends "nvhpc-$(echo "$NVHPC_VERSION_A" | sed 's/\./-/g')-cuda-multi"
fi
#
NVHPC_VERSION_B="$(echo "$NVHPC_VERSION_A" | sed 's/\./-/g')"
#
# save NVHPC version in GITHUB_ENV for later use when setting the NVHPC environment
#
echo "NVHPC_VERSION_A=$NVHPC_VERSION_A" >> $GITHUB_ENV
echo "NVHPC_VERSION_B=$NVHPC_VERSION_B" >> $GITHUB_ENV
