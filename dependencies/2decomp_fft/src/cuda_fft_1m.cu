//=======================================================================
// This is part of the 2DECOMP&FFT library
// 
// 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
// decomposition. It also implements a highly scalable distributed
// three-dimensional Fast Fourier Transform (FFT).
//
// Copyright (C) 2009-2021 Ning Li, the Numerical Algorithms Group (NAG)
//
//=======================================================================

// This contains CUDA code that compute multiple 1D FFTs on NVidia GPU

#ifdef DOUBLE_PREC
#define CUFFT_REAL_TYPE cufftDoubleReal
#define CUFFT_COMPLEX_TYPE cufftDoubleComplex
#define CUFFT_PLAN_TYPE_C2C CUFFT_Z2Z
#define CUFFT_PLAN_TYPE_R2C CUFFT_D2Z
#define CUFFT_PLAN_TYPE_C2R CUFFT_Z2D
#define CUFFT_EXEC_TYPE_C2C cufftExecZ2Z
#define CUFFT_EXEC_TYPE_R2C cufftExecD2Z
#define CUFFT_EXEC_TYPE_C2R cufftExecZ2D
#else
#define CUFFT_REAL_TYPE cufftReal
#define CUFFT_COMPLEX_TYPE cufftComplex
#define CUFFT_PLAN_TYPE_C2C CUFFT_C2C
#define CUFFT_PLAN_TYPE_R2C CUFFT_R2C
#define CUFFT_PLAN_TYPE_C2R CUFFT_C2R
#define CUFFT_EXEC_TYPE_C2C cufftExecC2C
#define CUFFT_EXEC_TYPE_R2C cufftExecR2C
#define CUFFT_EXEC_TYPE_C2R cufftExecC2R
#endif

#include <stdio.h>
#include <stdlib.h>
#include "cufft.h"
#include "cuda.h"

extern "C" void fft_1m_r2c_(int *nx, int *m, CUFFT_REAL_TYPE *h_a, CUFFT_COMPLEX_TYPE *h_b)
{
  unsigned long size1 = sizeof(CUFFT_REAL_TYPE) * (*nx) * (*m);
  unsigned long size2 = sizeof(CUFFT_COMPLEX_TYPE) * (*nx/2+1) * (*m);
  CUFFT_REAL_TYPE *d_ic = NULL;
  CUFFT_COMPLEX_TYPE *d_oc = NULL; 
  cufftHandle plan;
  cudaMalloc((void **)&d_ic, size1);
  cudaMalloc((void **)&d_oc, size2);
  cudaMemcpy(d_ic, h_a, size1, cudaMemcpyHostToDevice);
  int dims[1] = {*nx};
  cufftPlanMany(&plan,1,dims,NULL,1,0,NULL,1,0,CUFFT_PLAN_TYPE_R2C,*m);
  CUFFT_EXEC_TYPE_R2C(plan, d_ic, d_oc);
  cudaMemcpy(h_b, d_oc, size2, cudaMemcpyDeviceToHost);
  cudaFree(d_ic);
  cudaFree(d_oc);
  cufftDestroy(plan);
}


extern "C" void fft_1m_c2r_(int *nx, int *m, CUFFT_COMPLEX_TYPE *h_a, CUFFT_REAL_TYPE *h_b)
{
  unsigned long size1 = sizeof(CUFFT_COMPLEX_TYPE) * (*nx/2+1)*(*m);
  unsigned long size2 = sizeof(CUFFT_REAL_TYPE) * (*nx)*(*m);
  CUFFT_COMPLEX_TYPE *d_ic = NULL;
  CUFFT_REAL_TYPE *d_oc = NULL; 
  cufftHandle plan;
  cudaMalloc((void **)&d_ic, size1);
  cudaMalloc((void **)&d_oc, size2);
  cudaMemcpy(d_ic, h_a, size1, cudaMemcpyHostToDevice);
  int dims[1] = {*nx};
  cufftPlanMany(&plan,1,dims,NULL,1,0,NULL,1,0,CUFFT_PLAN_TYPE_C2R,*m);
  CUFFT_EXEC_TYPE_C2R(plan, d_ic, d_oc);
  cudaMemcpy(h_b, d_oc, size2, cudaMemcpyDeviceToHost);
  cudaFree(d_ic);
  cudaFree(d_oc);
  cufftDestroy(plan);
}


extern "C" void fft_1m_c2c_(int *nx, int *m, CUFFT_COMPLEX_TYPE *h_a, CUFFT_COMPLEX_TYPE *h_b, int *sign)
{
  unsigned long size1 = sizeof(CUFFT_COMPLEX_TYPE) * (*nx) * (*m);
  CUFFT_COMPLEX_TYPE *d_ic = NULL;
  CUFFT_COMPLEX_TYPE *d_oc = NULL; 
  cufftHandle plan;
  cudaMalloc((void **)&d_ic, size1);
  cudaMalloc((void **)&d_oc, size1);
  cudaMemcpy(d_ic, h_a, size1, cudaMemcpyHostToDevice);
  int dims[1] = {*nx};
  cufftPlanMany(&plan,1,dims,NULL,1,0,NULL,1,0,CUFFT_PLAN_TYPE_C2C,*m);
  CUFFT_EXEC_TYPE_C2C(plan, d_ic, d_oc, *sign);
  cudaMemcpy(h_b, d_oc, size1, cudaMemcpyDeviceToHost);
  cudaFree(d_ic);
  cudaFree(d_oc);
  cufftDestroy(plan);
}
