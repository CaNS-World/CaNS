#!/bin/sh
#
# installs Intel compiler (`ifx`) and MPI library, and sets up the build environment
# see: https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2023-0/apt.html
#
# get packages
#
sudo apt-get install libfftw3-dev
#
# download the key to system keyring
#
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
| gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt-get update
sudo apt-get install intel-hpckit
sudo apt-get install intel-fortran-essentials
