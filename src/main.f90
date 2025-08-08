! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
!
!        CCCCCCCCCCCCC                    NNNNNNNN        NNNNNNNN    SSSSSSSSSSSSSSS
!     CCC::::::::::::C                    N:::::::N       N::::::N  SS:::::::::::::::S
!   CC:::::::::::::::C                    N::::::::N      N::::::N S:::::SSSSSS::::::S
!  C:::::CCCCCCCC::::C                    N:::::::::N     N::::::N S:::::S     SSSSSSS
! C:::::C       CCCCCC   aaaaaaaaaaaaa    N::::::::::N    N::::::N S:::::S
!C:::::C                 a::::::::::::a   N:::::::::::N   N::::::N S:::::S
!C:::::C                 aaaaaaaaa:::::a  N:::::::N::::N  N::::::N  S::::SSSS
!C:::::C                          a::::a  N::::::N N::::N N::::::N   SS::::::SSSSS
!C:::::C                   aaaaaaa:::::a  N::::::N  N::::N:::::::N     SSS::::::::SS
!C:::::C                 aa::::::::::::a  N::::::N   N:::::::::::N        SSSSSS::::S
!C:::::C                a::::aaaa::::::a  N::::::N    N::::::::::N             S:::::S
! C:::::C       CCCCCC a::::a    a:::::a  N::::::N     N:::::::::N             S:::::S
!  C:::::CCCCCCCC::::C a::::a    a:::::a  N::::::N      N::::::::N SSSSSSS     S:::::S
!   CC:::::::::::::::C a:::::aaaa::::::a  N::::::N       N:::::::N S::::::SSSSSS:::::S
!     CCC::::::::::::C  a::::::::::aa:::a N::::::N        N::::::N S:::::::::::::::SS
!        CCCCCCCCCCCCC   aaaaaaaaaa  aaaa NNNNNNNN         NNNNNNN  SSSSSSSSSSSSSSS
!-------------------------------------------------------------------------------------
! CaNS -- Canonical Navier-Stokes Solver
!-------------------------------------------------------------------------------------
program cans
  use, intrinsic :: iso_fortran_env, only: compiler_version,compiler_options
  use, intrinsic :: iso_c_binding  , only: C_PTR
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan
  use mpi
  use decomp_2d
  use mod_bound          , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv         , only: chkdiv
  use mod_chkdt          , only: chkdt
  use mod_common_mpi     , only: myid,ierr,dinfo_ptdma
  use mod_correc         , only: correc
  use mod_fft            , only: fftini,fftend
  use mod_fillps         , only: fillps
  use mod_initflow       , only: initflow,initscal
  use mod_initgrid       , only: initgrid
  use mod_initmpi        , only: initmpi
  use mod_initsolver     , only: initsolver
  use mod_solve_helmholtz, only: solve_helmholtz,rhs_bound
#if _USE_HDF5
  use mod_load_hdf5      , only: load_one
#else
  use mod_load           , only: load_one
#endif
  use mod_mom            , only: bulk_forcing
  use mod_rk             , only: rk,rk_scal
  use mod_output         , only: out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d
  use mod_param          , only: ng,l,dl,dli, &
                                 gtype,gr, &
                                 cfl,dtmax,dt_f, &
                                 visc,alpha_max, &
                                 inivel,is_wallturb, &
                                 nstep,time_max,tw_max,stop_type, &
                                 restart,is_overwrite_save,nsaves_max, &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                                 cbcvel,bcvel,cbcpre,bcpre, &
                                 is_forced,bforce,velf, &
                                 gacc,nscal,beta, &
                                 dims, &
                                 nb,is_bound, &
                                 rkcoeff,small, &
                                 datadir, &
                                 read_input, &
                                 is_debug,is_debug_poisson, &
                                 is_timing, &
                                 is_impdiff,is_impdiff_1d, &
                                 is_poisson_pcr_tdma, &
                                 is_mask_divergence_check
  use mod_sanity         , only: test_sanity_input,test_sanity_solver
  use mod_scal           , only: scalar,initialize_scalars,bulk_forcing_s
#if !defined(_OPENACC)
  use mod_solver         , only: solver
#else
  use mod_solver_gpu     , only: solver => solver_gpu
  use mod_workspaces     , only: init_wspace_arrays,set_cufft_wspace,cudecomp_finalize
  use mod_common_cudecomp, only: istream_acc_queue_1,ap_z_ptdma
#endif
  use mod_timer          , only: timer_tic,timer_toc,timer_print
  use mod_updatep        , only: updatep
  use mod_utils          , only: bulk_mean
#if defined(_OPENACC)
  use mod_utils          , only: device_memory_footprint, arr_ptr
#endif
  use mod_types
  use omp_lib
  implicit none
  integer , dimension(3) :: lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
  real(rp), target, allocatable, dimension(:,:,:) :: u,v,w,p,pp
  real(rp), allocatable, dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
  real(rp), dimension(0:1,3) :: tauxo,tauyo,tauzo
  real(rp), dimension(3) :: f
#if !defined(_OPENACC)
  type(C_PTR), dimension(2,2) :: arrplanp
#else
  integer    , dimension(2,2) :: arrplanp
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:) :: ap,bp,cp
  integer , dimension(3) :: n_z_d
  real(rp), allocatable, dimension(:,:,:) :: ap_d,cp_d
  logical :: is_ptdma_update_p
  real(rp) :: normfftp
  type(rhs_bound) :: rhsbp
  real(rp) :: alpha
#if !defined(_OPENACC)
  type(C_PTR), dimension(2,2) :: arrplanu,arrplanv,arrplanw
#else
  integer    , dimension(2,2) :: arrplanu,arrplanv,arrplanw
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyu,lambdaxyv,lambdaxyw
  real(rp), allocatable, dimension(:) :: au,av,aw,bu,bv,bw,cu,cv,cw
  real(rp) :: normfftu,normfftv,normfftw
  type(rhs_bound) :: rhsbu,rhsbv,rhsbw
  !
  real(rp) :: dt,dti,dt_cfl,time,dtrk,dtrki,divtot,divmax
  integer :: irk,istep
  real(rp), allocatable, dimension(:) :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi, &
                                         dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g, &
                                         grid_vol_ratio_c,grid_vol_ratio_f
  real(rp) :: meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl
  !
  type(scalar), target, allocatable, dimension(:) :: scalars
  type(scalar), pointer :: s
  real(rp) :: meanscal
  real(rp), allocatable, dimension(:) :: fs
  integer :: iscal,is
  !
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(42) :: var
  !
  real(rp) :: dt12,dt12av,dt12min,dt12max
  real(rp) :: twi,tw
  !
  integer  :: savecounter
  character(len=7  ) :: fldnum
  character(len=3  ) :: scalnum
  character(len=4  ) :: chkptnum
  character(len=100) :: filename
  type(arr_ptr)    , allocatable, dimension(:) ::   io_vars
  character(len=10), allocatable, dimension(:) :: c_io_vars
  integer :: k,kk
  logical :: is_done,kill
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  !
  ! read parameter file
  !
  call read_input(myid)
  !
  ! initialize MPI/OpenMP
  !
  !$ call omp_set_num_threads(omp_get_max_threads())
  call initmpi(ng,dims,cbcpre,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,nb,is_bound)
  twi = MPI_WTIME()
  savecounter = 0
  !
  ! allocate variables
  !
  allocate(u( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           v( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           w( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           p( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)))
  allocate(lambdaxyp(n_z(1),n_z(2)))
  allocate(ap(n_z(3)),bp(n_z(3)),cp(n_z(3)))
  if(is_poisson_pcr_tdma) then
#if defined(_OPENACC)
    n_z_d(:) = ap_z_ptdma%shape(:)
#else
    n_z_d(:) = dinfo_ptdma%zsz(:)
#endif
    allocate(ap_d(n_z_d(1),n_z_d(2),n_z_d(3)), &
             cp_d(n_z_d(1),n_z_d(2),n_z_d(3)))
  end if
  allocate(dzc( 0:n(3)+1), &
           dzf( 0:n(3)+1), &
           zc(  0:n(3)+1), &
           zf(  0:n(3)+1), &
           dzci(0:n(3)+1), &
           dzfi(0:n(3)+1))
  allocate(dzc_g( 0:ng(3)+1), &
           dzf_g( 0:ng(3)+1), &
           zc_g(  0:ng(3)+1), &
           zf_g(  0:ng(3)+1), &
           dzci_g(0:ng(3)+1), &
           dzfi_g(0:ng(3)+1))
  allocate(grid_vol_ratio_c,mold=dzc)
  allocate(grid_vol_ratio_f,mold=dzf)
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
  if(is_impdiff) then
    allocate(lambdaxyu(n_z(1),n_z(2)), &
             lambdaxyv(n_z(1),n_z(2)), &
             lambdaxyw(n_z(1),n_z(2)))
    allocate(au(n_z(3)),bu(n_z(3)),cu(n_z(3)), &
             av(n_z(3)),bv(n_z(3)),cv(n_z(3)), &
             aw(n_z(3)),bw(n_z(3)),cw(n_z(3)))
    allocate(rhsbu%x(n(2),n(3),0:1), &
             rhsbu%y(n(1),n(3),0:1), &
             rhsbu%z(n(1),n(2),0:1), &
             rhsbv%x(n(2),n(3),0:1), &
             rhsbv%y(n(1),n(3),0:1), &
             rhsbv%z(n(1),n(2),0:1), &
             rhsbw%x(n(2),n(3),0:1), &
             rhsbw%y(n(1),n(3),0:1), &
             rhsbw%z(n(1),n(2),0:1))
  end if
  !
  allocate(scalars(nscal))
  call initialize_scalars(scalars,nscal,n,n_z)
  allocate(fs(nscal))
  !$acc enter data copyin(scalars(:))
  !
  if(is_debug) then
    if(myid == 0) print*, 'This executable of CaNS was built with compiler: ', compiler_version()
    if(myid == 0) print*, 'Using the options: ', compiler_options()
    block
      character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version
      integer :: ilen
      call MPI_GET_LIBRARY_VERSION(mpi_version,ilen,ierr)
      if(myid == 0) print*, 'MPI Version: ', trim(mpi_version)
    end block
    if(myid == 0) print*, ''
  end if
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, '*** Beginning of simulation ***'
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, ''
  call initgrid(gtype,ng(3),gr,l(3),dzc_g,dzf_g,zc_g,zf_g)
  if(myid == 0) then
    open(99,file=trim(datadir)//'grid.bin',action='write',form='unformatted',access='stream',status='replace')
    write(99) dzc_g(1:ng(3)),dzf_g(1:ng(3)),zc_g(1:ng(3)),zf_g(1:ng(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=0,ng(3)+1
      write(99,*) 0.,zf_g(kk),zc_g(kk),dzf_g(kk),dzc_g(kk)
    end do
    close(99)
    open(99,file=trim(datadir)//'geometry.out')
      write(99,*) ng(1),ng(2),ng(3)
      write(99,*) l(1),l(2),l(3)
    close(99)
  end if
  !$acc enter data copyin(lo,hi,n) async
  !$acc enter data copyin(bforce,dl,dli,l) async
  !$acc enter data copyin(gacc) async
  !$acc enter data copyin(zc_g,zf_g,dzc_g,dzf_g) async
  !$acc enter data create(zc,zf,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g) async
  !
  !$acc parallel loop default(present) private(k) async
  do kk=lo(3)-1,hi(3)+1
    k = kk-(lo(3)-1)
    zc( k) = zc_g(kk)
    zf( k) = zf_g(kk)
    dzc(k) = dzc_g(kk)
    dzf(k) = dzf_g(kk)
    dzci(k) = dzc(k)**(-1)
    dzfi(k) = dzf(k)**(-1)
  end do
  !$acc data copy(ng) async
  !$acc parallel loop default(present) async
  do k=0,ng(3)+1
    dzci_g(k) = dzc_g(k)**(-1)
    dzfi_g(k) = dzf_g(k)**(-1)
  end do
  !$acc end data
  !$acc enter data create(grid_vol_ratio_c,grid_vol_ratio_f) async
  !$acc parallel loop default(present) async
  do k=0,n(3)+1
    grid_vol_ratio_c(k) = dl(1)*dl(2)*dzc(k)/(l(1)*l(2)*l(3))
    grid_vol_ratio_f(k) = dl(1)*dl(2)*dzf(k)/(l(1)*l(2)*l(3))
  end do
  !$acc update self(zc,zf,dzc,dzf,dzci,dzfi) async
  !$acc exit data copyout(zc_g,zf_g,dzc_g,dzf_g,dzci_g,dzfi_g) async ! not needed on the device
  !$acc wait
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity_input(ng,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_forced)
  !
  ! initialize Poisson solver
  !
  call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:), &
                  lambdaxyp,['c','c','c'],ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
  !$acc enter data copyin(lambdaxyp,ap,bp,cp) async
  !$acc enter data copyin(rhsbp,rhsbp%x,rhsbp%y,rhsbp%z) async
  if(is_poisson_pcr_tdma) then
    !$acc enter data create(ap_d,cp_d) async
  end if
  !$acc wait
  if(is_impdiff) then
    call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,1),bcvel(:,:,1), &
                    lambdaxyu,['f','c','c'],au,bu,cu,arrplanu,normfftu,rhsbu%x,rhsbu%y,rhsbu%z)
    call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,2),bcvel(:,:,2), &
                    lambdaxyv,['c','f','c'],av,bv,cv,arrplanv,normfftv,rhsbv%x,rhsbv%y,rhsbv%z)
    call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,3),bcvel(:,:,3), &
                    lambdaxyw,['c','c','f'],aw,bw,cw,arrplanw,normfftw,rhsbw%x,rhsbw%y,rhsbw%z)
    do iscal=1,nscal
      s => scalars(iscal)
      call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,s%cbc,s%bc, &
                      s%lambdaxy,['c','c','c'],s%a,s%b,s%c,s%arrplan,s%normfft, &
                      s%rhsb%x,s%rhsb%y,s%rhsb%z)
    end do
    if(is_impdiff_1d) then
      deallocate(lambdaxyu,lambdaxyv,lambdaxyw)
      call fftend(arrplanu)
      call fftend(arrplanv)
      call fftend(arrplanw)
      deallocate(rhsbu%x,rhsbu%y,rhsbv%x,rhsbv%y,rhsbw%x,rhsbw%y)
      do iscal=1,nscal
        s => scalars(iscal)
        deallocate(s%rhsb%x,s%rhsb%y)
        deallocate(s%lambdaxy)
        call fftend(s%arrplan)
      end do
    end if
    !$acc enter data copyin(lambdaxyu,au,bu,cu,lambdaxyv,av,bv,cv,lambdaxyw,aw,bw,cw) async
    !$acc enter data copyin(rhsbu,rhsbu%x,rhsbu%y,rhsbu%z) async
    !$acc enter data copyin(rhsbv,rhsbv%x,rhsbv%y,rhsbv%z) async
    !$acc enter data copyin(rhsbw,rhsbw%x,rhsbw%y,rhsbw%z) async
    do iscal=1,nscal
      s => scalars(iscal)
      !$acc enter data copyin(s%lambdaxy,s%a,s%b,s%c) async
      !$acc enter data copyin(s%rhsb,s%rhsb%x,s%rhsb%y,s%rhsb%z) async
    end do
    !$acc wait
  end if
#if defined(_OPENACC)
  !
  ! determine workspace sizes and allocate the memory
  !
  call init_wspace_arrays()
  call set_cufft_wspace(pack(arrplanp,.true.),istream_acc_queue_1)
  if(is_impdiff .and. .not.is_impdiff_1d) then
    call set_cufft_wspace(pack(arrplanu,.true.),istream_acc_queue_1)
    call set_cufft_wspace(pack(arrplanv,.true.),istream_acc_queue_1)
    call set_cufft_wspace(pack(arrplanw,.true.),istream_acc_queue_1)
  end if
  if(myid == 0) print*,'*** Device memory footprint (Gb): ', &
                  device_memory_footprint(n,n_z,nscal)/(1._sp*1024**3), ' ***'
#endif
  if(is_debug_poisson) then
    call test_sanity_solver(ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,dli,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g, &
                            nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
  end if
  !
  allocate(io_vars(4+nscal),c_io_vars(4+nscal)) ! u,v,w,p,scalars
  io_vars(1)%arr => u; c_io_vars(1) = '_u'
  io_vars(2)%arr => v; c_io_vars(2) = '_v'
  io_vars(3)%arr => w; c_io_vars(3) = '_w'
  io_vars(4)%arr => p; c_io_vars(4) = '_p'
  do iscal=1,nscal
    write(scalnum,'(i3.3)') iscal
    io_vars(4+iscal)%arr => scalars(iscal)%val; c_io_vars(4+iscal) = '_s_'//scalnum
  end do
  !
  is_ptdma_update_p = .true.
  !
  if(.not.restart) then
    istep = 0
    time = 0.
    call initflow(inivel,bcvel,ng,lo,l,dl,zc,zf,dzc,dzf,visc,is_forced,velf,bforce,is_wallturb,u,v,w,p)
    do iscal=1,nscal
      s => scalars(iscal)
      call initscal(s%ini,s%bc,ng,lo,l,dl,zc,dzf,s%alpha,s%is_forced,s%scalf,s%val)
    end do
    if(myid == 0) print*, '*** Initial condition succesfully set ***'
  else
#ifdef _USE_HDF5
    call load_one('r',trim(datadir)//'fld.h5',c_io_vars, &
                  MPI_COMM_WORLD,ng,[1,1,1],lo,hi,io_vars,time,istep, 4+nscal)
#else
    do is=1,4+nscal
      call load_one('r',trim(datadir)//'fld'//trim(c_io_vars(is))//'.bin', &
                    MPI_COMM_WORLD,ng,[1,1,1],lo,hi,io_vars(is)%arr,time,istep)
    end do
#endif
    if(myid == 0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  end if
  !$acc enter data copyin(u,v,w,p,dudtrko,dvdtrko,dwdtrko) create(pp)
  call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
  do iscal=1,nscal
    s => scalars(iscal)
    !$acc enter data copyin(s%val,s%dsdtrko) async(1)
    call boundp(s%cbc,n,s%bc,nb,is_bound,dl,dzc,s%val)
  end do
  !$acc wait
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  !$acc wait ! not needed but to prevent possible future issues
  !$acc update self(u,v,w,p)
  do iscal=1,nscal
    !$acc update self(scalars(iscal)%val)
  end do
  if(iout1d > 0.and.mod(istep,max(iout1d,1)) == 0) then
    include 'out1d.h90'
  end if
  if(iout2d > 0.and.mod(istep,max(iout2d,1)) == 0) then
    include 'out2d.h90'
  end if
  if(iout3d > 0.and.mod(istep,max(iout3d,1)) == 0) then
    include 'out3d.h90'
  end if
  !
  call chkdt(n,dl,dzci,dzfi,visc,alpha_max,u,v,w,dt_cfl)
  dt = merge(dt_f,min(cfl*dt_cfl,dtmax),dt_f > 0.)
  if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
  dti = 1./dt
  kill = .false.
  !
  ! main loop
  !
  if(myid == 0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  do while(.not.is_done)
    if(is_timing) then
      !$acc wait(1)
      dt12 = MPI_WTIME()
    end if
    istep = istep + 1
    time = time + dt
    if(myid == 0) print*, 'Time step #', istep, 'Time = ', time
    tauxo(:,:) = 0.; tauyo(:,:) = 0.; tauzo(:,:) = 0.
    dpdl(:)     = 0.
    fs(1:nscal) = 0.
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      dtrki = dtrk**(-1)
      do iscal=1,nscal
        s => scalars(iscal)
        call rk_scal(rkcoeff(:,irk),n,dli,l,dzci,dzfi,grid_vol_ratio_f,s%alpha,dt,is_bound,u,v,w, &
                     s%is_forced,s%scalf,s%source,s%fluxo,s%dsdtrko,s%val,s%f)
        call bulk_forcing_s(n,s%is_forced,s%f,s%val)
        fs(iscal) = fs(iscal) + s%f
        if(is_impdiff) then
          call solve_helmholtz(n,ng,hi,s%arrplan,s%normfft,-0.5*s%alpha*dtrk, &
                               s%lambdaxy,s%a,s%b,s%c,s%rhsb%x,s%rhsb%y,s%rhsb%z,is_bound,s%cbc,['c','c','c'],s%val)
        end if
        call boundp(s%cbc,n,s%bc,nb,is_bound,dl,dzc,s%val)
      end do
      call rk(rkcoeff(:,irk),n,dli,dzci,dzfi,grid_vol_ratio_c,grid_vol_ratio_f,visc,dt,p, &
              is_forced,velf,bforce,gacc,beta,scalars,dudtrko,dvdtrko,dwdtrko,u,v,w,f)
      call bulk_forcing(n,is_forced,f,u,v,w)
      dpdl(:) = dpdl(:) + f(:)
      if(is_impdiff) then
        alpha = -.5*visc*dtrk
        call solve_helmholtz(n,ng,hi,arrplanu,normfftu,alpha, &
                             lambdaxyu,au,bu,cu,rhsbu%x,rhsbu%y,rhsbu%z,is_bound,cbcvel(:,:,1),['f','c','c'],u)
        call solve_helmholtz(n,ng,hi,arrplanv,normfftv,alpha, &
                             lambdaxyv,av,bv,cv,rhsbv%x,rhsbv%y,rhsbv%z,is_bound,cbcvel(:,:,2),['c','f','c'],v)
        call solve_helmholtz(n,ng,hi,arrplanw,normfftw,alpha, &
                             lambdaxyw,aw,bw,cw,rhsbw%x,rhsbw%y,rhsbw%z,is_bound,cbcvel(:,:,3),['c','c','f'],w)
      end if
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
      call fillps(n,dli,dzfi,dtrki,u,v,w,pp)
      call updt_rhs_b(['c','c','c'],cbcpre,n,is_bound,rhsbp%x,rhsbp%y,rhsbp%z,pp)
      call solver(n,ng,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre,['c','c','c'],pp,is_ptdma_update_p,ap_d,cp_d)
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,pp)
      call correc(n,dli,dzci,dtrk,pp,u,v,w)
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.true.,dl,dzc,dzf,u,v,w)
      call updatep(n,dli,dzci,dzfi,alpha,pp,p)
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
    end do
    dpdl(:)     = -dpdl(:)*dti
    fs(1:nscal) = fs(1:nscal)*dti
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep >= nstep   ) is_done = is_done.or..true.
    end if
    if(stop_type(2)) then ! maximum simulation time reached
      if(time  >= time_max) is_done = is_done.or..true.
    end if
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600.
      if(tw    >= tw_max  ) is_done = is_done.or..true.
    end if
    if(icheck > 0.and.mod(istep,max(icheck,1)) == 0) then
      if(myid == 0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,visc,alpha_max,u,v,w,dt_cfl)
      dt = merge(dt_f,min(cfl*dt_cfl,dtmax),dt_f > 0.)
      if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
      if(dt_cfl < small) then
        if(myid == 0) print*, 'ERROR: time step is too small.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
      dti = 1./dt
      call chkdiv(lo,hi,dli,dzfi,u,v,w,divtot,divmax)
      if(myid == 0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
      if(.not.is_mask_divergence_check) then
        if(divmax > small.or.is_nan(divtot)) then
          if(myid == 0) print*, 'ERROR: maximum divergence is too large.'
          if(myid == 0) print*, 'Aborting...'
          is_done = .true.
          kill = .true.
        end if
      end if
    end if
    !
    ! output routines below
    !
    if(iout0d > 0.and.mod(istep,max(iout0d,1)) == 0) then
      !allocate(var(4))
      var(1) = 1.*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      !
      if(any(is_forced(:)).or.any(abs(bforce(:)) > 0.)) then
        meanvelu = 0.
        meanvelv = 0.
        meanvelw = 0.
        if(is_forced(1).or.abs(bforce(1)) > 0.) then
          call bulk_mean(n,grid_vol_ratio_f,u,meanvelu)
        end if
        if(is_forced(2).or.abs(bforce(2)) > 0.) then
          call bulk_mean(n,grid_vol_ratio_f,v,meanvelv)
        end if
        if(is_forced(3).or.abs(bforce(3)) > 0.) then
          call bulk_mean(n,grid_vol_ratio_c,w,meanvelw)
        end if
        if(.not.any(is_forced(:))) dpdl(:) = -bforce(:) ! constant pressure gradient
        var(1)   = time
        var(2:4) = dpdl(1:3)
        var(5:7) = [meanvelu,meanvelv,meanvelw]
        call out0d(trim(datadir)//'forcing.out',7,var)
      end if
      !
      do iscal=1,nscal
        s => scalars(iscal)
        write(scalnum,'(i3.3)') iscal
        if(s%is_forced.or.abs(s%source) > 0.) then
          meanscal = 0.
          call bulk_mean(n,grid_vol_ratio_f,s%val,meanscal)
          if(.not.s%is_forced) fs(:) = s%source
          var(1:3) = [time,fs(iscal),meanscal]
          call out0d(trim(datadir)//'forcing_s_'//scalnum//'.out',3,var)
        end if
      end do
    end if
    write(fldnum,'(i7.7)') istep
    if(iout1d > 0.and.mod(istep,max(iout1d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out1d.h90'
    end if
    if(iout2d > 0.and.mod(istep,max(iout2d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out2d.h90'
    end if
    if(iout3d > 0.and.mod(istep,max(iout3d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out3d.h90'
    end if
    if(isave > 0.and.((mod(istep,max(isave,1)) == 0).or.(is_done.and..not.kill))) then
      if(is_overwrite_save) then
        filename = 'fld'
      else
        filename = 'fld_'//fldnum
        if(nsaves_max > 0) then
          if(savecounter >= nsaves_max) savecounter = 0
          savecounter = savecounter + 1
          write(chkptnum,'(i4.4)') savecounter
          filename = 'fld_'//chkptnum
          var(1) = 1.*istep
          var(2) = time
          var(3) = 1.*savecounter
          call out0d(trim(datadir)//'log_checkpoints.out',3,var)
        end if
      end if
      !$acc wait
      !$acc update self(u,v,w,p)
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
#ifdef _USE_HDF5
      call load_one('w',trim(datadir)//'fld.h5',c_io_vars, &
                  MPI_COMM_WORLD,ng,[1,1,1],lo,hi,io_vars,time,istep, 4+nscal)
#else
      do is=1,4+nscal
        call load_one('w',trim(datadir)//trim(filename)//trim(c_io_vars(is))//'.bin', &
                      MPI_COMM_WORLD,ng,[1,1,1],lo,hi,io_vars(is)%arr,time,istep)
      end do
#endif
      if(.not.is_overwrite_save) then
        !
        ! fld_*.bin -> last checkpoint file (symbolic link)
        !
        do is=1,4+nscal
          call gen_alias(myid,trim(datadir),trim(filename)//trim(c_io_vars(is))//'.bin','fld'//trim(c_io_vars(is))//'.bin')
        end do
      end if
      if(myid == 0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
    end if
    if(is_timing) then
      !$acc wait(1)
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid == 0) print*, 'Avrg, min & max elapsed time: '
      if(myid == 0) print*, dt12av/(1.*product(dims)),dt12min,dt12max
    end if
  end do
  !
  ! clear ffts
  !
  call fftend(arrplanp)
  if(is_impdiff .and. .not.is_impdiff_1d) then
    call fftend(arrplanu)
    call fftend(arrplanv)
    call fftend(arrplanw)
  end if
  if(myid == 0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
#if defined(_OPENACC)
  call cudecomp_finalize
#endif
  call MPI_FINALIZE(ierr)
end program cans
