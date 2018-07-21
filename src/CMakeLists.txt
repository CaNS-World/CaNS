SET ( SRC
  param.f90
  common_mpi.f90 
  bound.f90 
  chkdiv.f90 
  chkdt.f90 
  correc.f90 
  debug.f90 
  fftw.f90
  fft.f90
  fillps.f90
  2decomp/decomp_2d.f90
  initflow.f90 
  initgrid.f90 
  initmpi.f90
  initsolver.f90 
  2decomp/io.f90
  load.f90 
  output.f90
  mom.f90
  momd.f90
  moms.f90
  rk.f90
  solver.f90
  sanity.f90
  main.f90
  )

ADD_EXECUTABLE ( cans
  ${DECOMP_SRC}
  ${SRC}
  )

TARGET_LINK_LIBRARIES ( cans
  fftw3
  fftw3_threads
  lapack
  )
