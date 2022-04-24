# 2DECOMP&FFT

2DECOMP&FFT is a library for 2D pencil decomposition and highly scalable distributed 3D Fast Fourier Transforms

#### Table of Contents

- [Overview](doc/overview.md)
- Software
  - [Download](doc/download.md)
  - [Installation](doc/installation.md)
- [Domain decomposition strategies](doc/decomposition.md)
- [Fast Fourier Transform (FFT) review](doc/fft.md)
- APIs
  - [2D pencil decomposition APIs](doc/api_decomposition.md)
  - [FFT APIs](doc/api_fft.md)
  - [Halo cell support](doc/api_halo.md)
  - [Parallel I/O](doc/api_io.md)
  - [Non-blocking communication](doc/api_nonblocking.md)
- Performance benchmarks
  - [2DECOMP&FFT vs. P3DFFT](doc/p3dfft.md)
  - [HECToR](doc/hector.md)
  - [JUGENE](doc/jugene.md)
- Applications and case studies
  - [Sample applications](doc/samples.md)
  - [Case study - Vortex generation using FFT](doc/vortex.md)
  - [Incompact3D - a CFD application for turbulence research](doc/incompact3d.md)
  - [DSTAR - a CFD application for studies of turbulence, aeroacoustics, combustion and multiphase flow](doc/dstar.md)
- Miscellaneous technical subjects
  - [Interactive decomposition map](https://monet.nag.co.uk/2decomp/decomp_map.php) 
  - [Using the 1D slab decompostion mode](doc/1d_decomp.md)
  - [Shared-memory optimisation](doc/shared_memory.md)
  - [Process grid](doc/pgrid.md)
  - [Padded all-to-all optimisation](doc/padded_alltoall.md)
  - [Precision guidelines](doc/precision.md)
  - [Memory comsumption](doc/memory.md)

#### Software License

Copyright &copy; 2011-2021, The Numerical Algorithms Group (NAG)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of the copyright owner nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.