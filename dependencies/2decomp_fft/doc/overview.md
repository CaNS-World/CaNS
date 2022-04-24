## Overview

### Introduction

The 2DECOMP&FFT library is a software framework in Fortran to build large-scale parallel applications. It is designed for applications using three-dimensional structured mesh and spatially implicit numerical algorithms. At the foundation it implements a general-purpose 2D pencil decomposition for data distribution on distributed-memory platforms. On top, it provides a highly scalable and efficient interface to perform three-dimensional distributed FFTs. The library is optimised for supercomputers and scales well to hundreds of thousands of cores. It relies on MPI but provides a user-friendly programming interface that hides communication details from application developers.

### Features

Here is a list of 2DECOMP&FFT's main features:

* General-purpose 2D pencil decomposition module to support building large-scale parallel applications on distributed memory systems.
* Highly scalable and efficient distributed Fast Fourier Transform module, supporting three dimensional FFTs (both complex-to-complex and real-to-complex/complex-to-real).
* Halo-cell support allowing explicit message passing between neighbouring blocks.
* Parallel I/O module to support the handling of large data sets.
* Shared-memory optimisation on the communication code for multi-core systems.

2DECOMP&FFT distinguishes itself from many other popular distributed FFT libraries by exposing its communication APIs upon which many other parallel algorithms can be built.

2DECOMP&FFT is designed to be:

* **Scalable** - The library and applications built upon it are known to scale to o(10^5) cores on major supercomputers.
* **Flexible** - Software framework to support building higher-level libraries and many types of applications.
* **User-friendly** - Black-box implementation and very clean application programming interface hiding most communication details from applications.
* **Portable** - Code tested on many major supercomputing architectures. The FFT library interfaces with almost every popular external FFT implementations.

### History

This software package was originally derived from several projects funded under the HECToR Distributed Computational Science and Engineering (dCSE) programme operated by NAG Ltd. HECToR - a UK Research Councils' high end computing service - served as the UK's national supercomputer for open science between 2008 and 2014.

The active development of this library completed in 2012. It has been in production use in many research applications since then. The code quality appears to be very good with almost no major bugs reported over the years. Its performance remains very competitive as reported by a [recent study](https://www.icl.utk.edu/files/publications/2021/icl-utk-1490-2021.pdf). 

Since August 2021, this project is hosted in NAG's official GitHub account to facilitate future development and maintenance.

### Citation

If you wish to cite this work, you are recommended to use the following paper:

* N. Li and S. Laizet, "2DECOMP&FFT – A highly scalable 2D decomposition library and FFT interface", Cray User Group 2010 conference, Edinburgh, 2010.