! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_sanity
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound     , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv    , only: chkdiv
  use mod_common_mpi, only: myid,ierr
  use mod_correc    , only: correc
  use mod_debug     , only: chk_helmholtz
  use mod_fft       , only: fftend
  use mod_fillps    , only: fillps
  use mod_initflow  , only: add_noise
  use mod_initmpi   , only: initmpi
  use mod_initsolver, only: initsolver
  use mod_param     , only: ipencil_axis,is_impdiff,is_impdiff_1d,is_poisson_dtdma,small, &
                            is_bulk_forcing_wall_shear,is_boussinesq_buoyancy,gacc
#if !defined(_OPENACC)
  use mod_solver    , only: solver
#else
  use mod_solver_gpu, only: solver => solver_gpu
#endif
  use mod_types
  implicit none
  private
  public test_sanity_input,test_sanity_solver,test_sanity_halo
  contains
  subroutine test_sanity_input(ng,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_forced)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    integer , intent(in), dimension(3) :: ng
    integer , intent(in), dimension(2) :: dims
    logical , intent(in), dimension(3) :: stop_type
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
    real(rp)        , intent(in), dimension(0:1,3,3) :: bcvel
    real(rp)        , intent(in), dimension(0:1,3)   :: bcpre
    logical         , intent(in), dimension(3)       :: is_forced
    logical :: passed,passed_loc
    !
    call chk_dims(ng,dims,cbcvel,cbcpre,passed);      if(.not.passed) call abortit
    call chk_stop_type(stop_type,passed);             if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed);    if(.not.passed) call abortit
    call chk_forcing(cbcvel,cbcpre,is_forced,passed); if(.not.passed) call abortit
    !
    if(is_poisson_dtdma) then
      passed_loc = ng(3)/dims(2) >= 2
      if(is_impdiff.and.(.not.is_impdiff_1d .or. dims(2) > 1)) then
        if(cbcvel(1,3,3) /= 'P') passed_loc = passed_loc.and.(ng(3)/dims(2)-1 >= 2)
      end if
      if(myid == 0.and.(.not.passed_loc)) &
        print*, 'ERROR: DTDMA requires at least two active points per Z slab.'
      passed = passed.and.passed_loc
      if(.not.passed) call abortit
    end if
    !
    if(is_impdiff_1d.and.(.not.is_impdiff)) then
      if(myid == 0)  print*, 'ERROR: `is_impdiff_1d = T` requires `is_impdiff = T` (forced in `param.f90`).'; call abortit
    end if
    if(is_impdiff_1d .and. .not.(ipencil_axis == 3) .and. .not.is_poisson_dtdma) then
      if(dims(2) > 1) then
        if(myid == 0)  print*, 'Warning: a run with implicit Z diffusion (`is_impdiff_1d = T`) is much more efficient &
                                       & when the flow is not decomposed along the Z direction.'
      end if
    end if
    if(is_poisson_dtdma.and.(ipencil_axis == 3)) then
      if(myid == 0)  print*, 'ERROR: `is_poisson_dtdma = T` requires X/Y-aligned pencils.'; call abortit
    end if
  end subroutine test_sanity_input
  !
  subroutine chk_stop_type(stop_type,passed)
    implicit none
    logical, intent(in), dimension(3) :: stop_type
    logical, intent(out) :: passed
    passed = .true.
    if(.not.any(stop_type(:))) then
      if(myid == 0) print*, 'ERROR: stopping criterion not chosen.'
      passed = .false.
    end if
  end subroutine chk_stop_type
  !
  subroutine chk_dims(ng,dims,cbcvel,cbcpre,passed)
    use mod_param, only: nscal,cbcscal
    implicit none
    integer, intent(in), dimension(3) :: ng
    integer, intent(in), dimension(2) :: dims
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
    logical, intent(out) :: passed
    integer, dimension(2) :: ii
    logical :: passed_loc
    integer :: idir,iscal
    character(len=2) :: bc01v
    passed = .true.
    passed_loc = all(ng(:) >= 1)
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: ng(:) must be positive in every direction.'
    passed = passed.and.passed_loc
    !
    ii = pack([1,2,3],[1,2,3] /= ipencil_axis)
    passed_loc = all(dims(:)<=ng(ii)).and.all(dims(:)>=1)
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: 1 <= dims(:) <= [ng(1),ng(2)], or [ng(1),ng(3)], or [ng(2),ng(3)] depending on the decomposition.'
    passed = passed.and.passed_loc
    !
    if(is_impdiff) then
      passed_loc = .true.
      do idir=1,2
        if(ng(idir) /= 1) cycle
        bc01v = cbcvel(0,idir,idir)//cbcvel(1,idir,idir)
        passed_loc = passed_loc.and.bc01v /= 'NN'
#if !defined(_OPENACC)
        if(.not.is_impdiff_1d) passed_loc = passed_loc.and.(bc01v == 'PP'.or.bc01v == 'ND')
#endif
      end do
      if(myid == 0.and.(.not.passed_loc)) &
        print*, 'ERROR: unsupported one-point X/Y direction for the selected implicit velocity BCs.'
      passed = passed.and.passed_loc
    end if
    !
    if(ng(3) == 1) then
      if(cbcpre(0,3)//cbcpre(1,3) == 'PP') then
        passed_loc = all(cbcvel(:,3,:) == 'P')
        do iscal=1,nscal
          passed_loc = passed_loc.and.all(cbcscal(:,3,iscal) == 'P')
        end do
      else
#if defined(_OPENACC)
        passed_loc = .false.
#else
        passed_loc = .not.is_impdiff
        if(cbcpre(0,3)//cbcpre(1,3) == 'NN') passed_loc = passed_loc.and.any(cbcpre(:,1:2) == 'D')
#endif
      end if
      if(myid == 0.and.(.not.passed_loc)) &
        print*, 'ERROR: unsupported one-point Z configuration; use periodic BCs for all fields.'
      passed = passed.and.passed_loc
    end if
  end subroutine chk_dims
  !
  subroutine chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed)
    use mod_param, only: nscal,cbcscal
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
    real(rp)        , intent(in), dimension(0:1,3,3) :: bcvel
    real(rp)        , intent(in), dimension(0:1,3  ) :: bcpre
    logical         , intent(out) :: passed
    character(len=2) :: bc01v,bc01p
    integer :: ivel,idir,iscal
#if defined(_OPENACC)
    character(len=2) :: bc01s
#endif
    logical :: passed_loc
    passed = .true.
    !
    ! check validity of pressure and velocity BCs
    !
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,3
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        passed_loc = passed_loc.and.( (bc01v == 'PP').or. &
                                      (bc01v == 'ND').or. &
                                      (bc01v == 'DN').or. &
                                      (bc01v == 'NN').or. &
                                      (bc01v == 'DD') )
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity BCs not valid.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01p == 'PP').or. &
                                    (bc01p == 'ND').or. &
                                    (bc01p == 'DN').or. &
                                    (bc01p == 'NN').or. &
                                    (bc01p == 'DD') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: pressure BCs not valid.'
    passed = passed.and.passed_loc
    !
    ! check that all variables have the same periodic directions
    !
    passed_loc = .true.
    do idir=1,3
      do ivel=1,3
        passed_loc = passed_loc.and.all((cbcvel(:,idir,ivel) == 'P').eqv.(cbcpre(:,idir) == 'P'))
      end do
      do iscal=1,nscal
        passed_loc = passed_loc.and.all((cbcscal(:,idir,iscal) == 'P').eqv.(cbcpre(:,idir) == 'P'))
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: velocity and scalar periodic BCs must match pressure BCs in every direction.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      ivel = idir
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01v == 'PP'.and.bc01p == 'PP').or. &
                                    (bc01v == 'ND'.and.bc01p == 'DN').or. &
                                    (bc01v == 'DN'.and.bc01p == 'ND').or. &
                                    (bc01v == 'DD'.and.bc01p == 'NN').or. &
                                    (bc01v == 'NN'.and.bc01p == 'DD') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity and pressure BCs not compatible.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,2
      passed_loc = passed_loc.and.((bcpre(0,idir) == 0.).and.(bcpre(1,idir) == 0.))
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: pressure BCs in directions x and y must be homogeneous (value = 0.).'
    passed = passed.and.passed_loc
    if(is_impdiff .and. .not.is_impdiff_1d) then
      passed_loc = .true.
      do ivel = 1,3
        do idir=1,2
          passed_loc = passed_loc.and.((bcvel(0,idir,ivel) == 0.).and.(bcvel(1,idir,ivel) == 0.))
        end do
      end do
      if(myid == 0.and.(.not.passed_loc)) &
        print*, 'ERROR: velocity BCs with implicit diffusion in directions x and y must be homogeneous (value = 0.).'
      passed = passed.and.passed_loc
    end if
#if defined(_OPENACC)
    if(is_impdiff .and. .not.is_impdiff_1d) then
      passed_loc = .true.
      do iscal=1,nscal
        do idir=1,2
          bc01s = cbcscal(0,idir,iscal)//cbcscal(1,idir,iscal)
          passed_loc = passed_loc.and.( (bc01s == 'PP').or. &
                                        (bc01s == 'ND').or. &
                                        (bc01s == 'DN').or. &
                                        (bc01s == 'NN').or. &
                                        (bc01s == 'DD') )
        end do
      end do
      if(myid == 0.and.(.not.passed_loc)) &
        print*, 'ERROR: scalar BCs along x or y not valid.'
      passed = passed.and.passed_loc
    end if
#endif
  end subroutine chk_bc
  !
  subroutine chk_forcing(cbcvel,cbcpre,is_forced,passed)
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre
    logical         , intent(in), dimension(3) :: is_forced
    logical         , intent(out) :: passed
    logical :: passed_loc
    integer :: idir
    passed = .true.
    !
    ! 1) check for compatibility between pressure BCs and flow forcing
    !
    do idir=1,3
      if(is_forced(idir)) then
        passed = passed.and.(cbcpre(0,idir)//cbcpre(1,idir) == 'PP')
      end if
    end do
    if(myid == 0.and.(.not.passed)) &
    print*, 'ERROR: Flow cannot be forced in a non-periodic direction; check the BCs and is_forced in `input.nml`.'
    !
    ! 2) check the requirements of the bulk forcing that balances the wall shear
    !
    if(is_bulk_forcing_wall_shear) then
      passed_loc = .true.
      do idir=1,3
        if(is_forced(idir)) passed_loc = passed_loc.and..not.any(cbcvel(:,:,idir) == 'N')
      end do
      passed_loc = passed_loc.and..not.(is_impdiff.and..not.is_impdiff_1d.and.any(cbcpre(:,1:2) /= 'P'))
      passed_loc = passed_loc.and..not.(is_boussinesq_buoyancy.and.any(is_forced(:).and.abs(gacc(:)) > 0.))
      if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: `is_bulk_forcing_wall_shear = T` requires no-slip walls for the forced velocity components, &
              &no fully implicit diffusion with walls along x or y, and no gravity along a forced direction with &
              &Boussinesq buoyancy.'
      passed = passed.and.passed_loc
    end if
  end subroutine chk_forcing
  !
  subroutine test_sanity_solver(ng,lo,hi,n,nh_v,nh_p,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,n_z,l,dl,dli,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g, &
                                nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
#if defined(_OPENACC)
    use mod_workspaces     , only: set_cufft_wspace
    use mod_common_cudecomp, only: istream_acc_queue_1
#endif
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
    integer , intent(in)               :: nh_v,nh_p,nh_g
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(3) :: dl,dli
    real(rp), intent(in), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi,dzci_g,dzfi_g
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
    real(rp), intent(in), dimension(0:1,3,3)          :: bcvel
    real(rp), intent(in), dimension(0:1,3)            :: bcpre
    real(rp), allocatable, dimension(:,:,:) :: u,v,w,p
#if !defined(_OPENACC) || defined(_USE_HIP)
    type(C_PTR), dimension(2,2) :: arrplan
#else
    integer    , dimension(2,2) :: arrplan
#endif
    real(rp) :: normfft
    real(rp), allocatable, dimension(:,:) ::  lambdaxy
    real(rp), allocatable, dimension(:) :: a,b,c,bb
    real(rp), allocatable, dimension(:,:,:) :: rhsbx,rhsby,rhsbz
    real(rp) :: dt,dti,alpha
    real(rp) :: divtot,divmax,restot,resmax
    integer :: i,j,k
    logical :: passed,passed_loc
    passed = .true.
    !$acc wait
    allocate(u(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v), &
             v(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v), &
             w(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v), &
             p(1-nh_p:n(1)+nh_p,1-nh_p:n(2)+nh_p,1-nh_p:n(3)+nh_p), &
             lambdaxy(n_z(1),n_z(2)), &
             a(n_z(3)),b(n_z(3)),c(n_z(3)), &
             rhsbx(n(2),n(3),0:1), &
             rhsby(n(1),n(3),0:1), &
             rhsbz(n(1),n(2),0:1))
    !$acc enter data copyin(n,n_z)
    !
    ! initialize velocity below with some random noise
    !
    u(:,:,:) = 0.
    v(:,:,:) = 0.
    w(:,:,:) = 0.
    p(:,:,:) = 0.
    call add_noise(ng,lo,123,.5_rp,u(1:n(1),1:n(2),1:n(3)))
    call add_noise(ng,lo,456,.5_rp,v(1:n(1),1:n(2),1:n(3)))
    call add_noise(ng,lo,789,.5_rp,w(1:n(1),1:n(2),1:n(3)))
    !$acc enter data copyin(u,v,w,p)
    !
    ! test pressure correction
    !
    call initsolver(ng,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:), &
                    lambdaxy,['c','c','c'],a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
    !$acc enter data copyin(lambdaxy,a,b,c,rhsbx,rhsby,rhsbz)
#if defined(_OPENACC)
    call set_cufft_wspace(pack(arrplan,.true.),istream_acc_queue_1)
#endif
    dt  = acos(-1.) ! value is irrelevant
    dti = dt**(-1)
    call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w,.false.,-1)
    call fillps(n,nh_v,nh_p,nh_g,dli,dzfi,dti,u,v,w,p)
    call updt_rhs_b(['c','c','c'],cbcpre,n,nh_p,is_bound,rhsbx,rhsby,rhsbz,p)
    call solver(n,nh_p,ng,arrplan,normfft,lambdaxy,a,b,c,cbcpre,['c','c','c'],p)
    call boundp(cbcpre,n,nh_p,nh_g,bcpre,nb,is_bound,dl,dzc,p)
    call correc(n,nh_v,nh_p,nh_g,dli,dzci,dt,p,u,v,w)
    call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w,.true.,1,0)
    call chkdiv(lo,hi,nh_v,nh_g,l,dli,dzfi,u,v,w,divtot,divmax)
    passed_loc = divmax < small
    if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: Pressure correction: Divergence is too large, with maximum = ', divmax
    passed = passed.and.passed_loc
    call fftend(arrplan)
    if(is_impdiff .and. .not.is_impdiff_1d) then
      allocate(bb(n_z(3)))
      alpha = -acos(-1.) ! negative shift, as in implicit diffusion
      !$acc parallel loop collapse(3) default(present)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=0,n(3)+1
        do j=0,n(2)+1
          do i=0,n(1)+1
            u(i,j,k) = 0.
            v(i,j,k) = 0.
            w(i,j,k) = 0.
          end do
        end do
      end do
      !$acc update self(u,v,w)
      call add_noise(ng,lo,123,.5_rp,u(1:n(1),1:n(2),1:n(3)))
      call add_noise(ng,lo,456,.5_rp,v(1:n(1),1:n(2),1:n(3)))
      call add_noise(ng,lo,789,.5_rp,w(1:n(1),1:n(2),1:n(3)))
      !$acc update device(u,v,w)
      !$acc enter data create(bb)
      call initsolver(ng,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,1),bcvel(:,:,1), &
                      lambdaxy,['f','c','c'],a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
      !$acc update device(lambdaxy,a,b,c,rhsbx,rhsby,rhsbz)
#if defined(_OPENACC)
      call set_cufft_wspace(pack(arrplan,.true.),istream_acc_queue_1)
#endif
      call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w)
      !$acc parallel loop collapse(3) default(present)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=0,n(3)+1
        do j=0,n(2)+1
          do i=0,n(1)+1
            u(i,j,k) = u(i,j,k)*alpha
            p(i,j,k) = u(i,j,k)
          end do
        end do
      end do
      !$acc parallel loop default(present)
      !$OMP parallel do   DEFAULT(shared)
      do k=1,n_z(3)
        bb(k) = b(k) + alpha
      end do
      call updt_rhs_b(['f','c','c'],cbcvel(:,:,1),n,nh_v,is_bound,rhsbx,rhsby,rhsbz,u)
      call solver(n,nh_v,ng,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,:,1),['f','c','c'],u)
      call fftend(arrplan)
      ! actually, we are only interested in the boundary condition in `u`
      call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w)
      call chk_helmholtz(lo,hi,nh_p,nh_v,nh_g,l,dli,dzci,dzfi,alpha,p,u,cbcvel(:,:,1),is_bound,['f','c','c'],restot,resmax)
      passed_loc = resmax < small
      if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: wrong solution of Helmholtz equation in x direction.'
      passed = passed.and.passed_loc
      !
      call initsolver(ng,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,2),bcvel(:,:,2), &
                      lambdaxy,['c','f','c'],a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
      !$acc update device(lambdaxy,a,b,c,rhsbx,rhsby,rhsbz)
#if defined(_OPENACC)
      call set_cufft_wspace(pack(arrplan,.true.),istream_acc_queue_1)
#endif
      call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w)
      !$acc parallel loop collapse(3) default(present)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=0,n(3)+1
        do j=0,n(2)+1
          do i=0,n(1)+1
            v(i,j,k) = v(i,j,k)*alpha
            p(i,j,k) = v(i,j,k)
          end do
        end do
      end do
      !$acc parallel loop default(present)
      !$OMP parallel do   DEFAULT(shared)
      do k=1,n_z(3)
        bb(k) = b(k) + alpha
      end do
      call updt_rhs_b(['c','f','c'],cbcvel(:,:,2),n,nh_v,is_bound,rhsbx,rhsby,rhsbz,v)
      call solver(n,nh_v,ng,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,:,2),['c','f','c'],v)
      call fftend(arrplan)
      ! actually, we are only interested in the boundary condition in `v`
      call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w)
      call chk_helmholtz(lo,hi,nh_p,nh_v,nh_g,l,dli,dzci,dzfi,alpha,p,v,cbcvel(:,:,2),is_bound,['c','f','c'],restot,resmax)
      passed_loc = resmax < small
      if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: wrong solution of Helmholtz equation in y direction.'
      passed = passed.and.passed_loc
      !
      call initsolver(ng,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,3),bcvel(:,:,3), &
                      lambdaxy,['c','c','f'],a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
      !$acc update device(lambdaxy,a,b,c,rhsbx,rhsby,rhsbz)
#if defined(_OPENACC)
      call set_cufft_wspace(pack(arrplan,.true.),istream_acc_queue_1)
#endif
      call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w)
      !$acc parallel loop collapse(3) default(present)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=0,n(3)+1
        do j=0,n(2)+1
          do i=0,n(1)+1
            w(i,j,k) = w(i,j,k)*alpha
            p(i,j,k) = w(i,j,k)
          end do
        end do
      end do
      !$acc parallel loop default(present)
      !$OMP parallel do   DEFAULT(shared)
      do k=1,n_z(3)
        bb(k) = b(k) + alpha
      end do
      call updt_rhs_b(['c','c','f'],cbcvel(:,:,3),n,nh_v,is_bound,rhsbx,rhsby,rhsbz,w)
      call solver(n,nh_v,ng,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,:,3),['c','c','f'],w)
      call fftend(arrplan)
      ! actually, we are only interested in the boundary condition in `w`
      call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w)
      call chk_helmholtz(lo,hi,nh_p,nh_v,nh_g,l,dli,dzci,dzfi,alpha,p,w,cbcvel(:,:,3),is_bound,['c','c','f'],restot,resmax)
      passed_loc = resmax < small
      if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: wrong solution of Helmholtz equation in z direction.'
      passed = passed.and.passed_loc
      !$acc exit data delete(bb)
    end if
    !$acc exit data delete(u,v,w,p,lambdaxy,a,b,c,rhsbx,rhsby,rhsbz)
    if(.not.passed) then
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      error stop
    end if
  end subroutine test_sanity_solver
  !
  subroutine abortit
    implicit none
    if(myid == 0) print*, ''
    if(myid == 0) print*, '*** Simulation aborted due to errors in the input file ***'
    if(myid == 0) print*, '    check `input.nml`.'
    call decomp_2d_finalize
    call MPI_FINALIZE(ierr)
    error stop
  end subroutine abortit
  !
  subroutine test_sanity_halo(ng,lo,n,nh_v,nh_p,nh_s,nh_g,dl,dzc,dzf,zc_g,zf_g,nb,is_bound, &
                              cbcvel,bcvel,cbcpre,bcpre,cbcscal,bcscal)
    !
    ! checks the halo exchanges and boundary conditions on all ghost layers, edges and
    ! corners included: the interior of each field is set to f(i,j,k), a function of the
    ! global indices, the halos are updated, and every point is compared against its
    ! expected value (see `expected_value`)
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,n
    integer , intent(in)               :: nh_v,nh_p,nh_s,nh_g
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(1-nh_g:) :: dzc,dzf,zc_g,zf_g
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    real(rp)        , intent(in), dimension(0:1,3,3) :: bcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
    real(rp)        , intent(in), dimension(0:1,3)   :: bcpre
    character(len=1), intent(in), dimension(0:1,3), optional :: cbcscal
    real(rp)        , intent(in), dimension(0:1,3), optional :: bcscal
    real(rp), allocatable, dimension(:,:,:) :: u,v,w,p
    real(rp), dimension(3) :: errmax
    logical , dimension(3,3) :: is_c
    logical :: passed
    integer :: ivel
    !
    is_c(:,:) = .true.
    do ivel=1,3
      is_c(ivel,ivel) = .false. ! velocity component `ivel` is staggered along `ivel`
    end do
    passed = .true.
    !
    ! velocity
    !
    allocate(u(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v), &
             v(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v), &
             w(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v))
    call set_field(n,nh_v,lo,u)
    call set_field(n,nh_v,lo,v)
    call set_field(n,nh_v,lo,w)
    !$acc enter data copyin(u,v,w)
    call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w)
    !$acc wait
    !$acc exit data copyout(u,v,w)
    errmax(1) = field_error(ng,lo,n,nh_v,nh_g,is_c(:,1),cbcvel(:,:,1),bcvel(:,:,1),dl,zc_g,zf_g,u)
    errmax(2) = field_error(ng,lo,n,nh_v,nh_g,is_c(:,2),cbcvel(:,:,2),bcvel(:,:,2),dl,zc_g,zf_g,v)
    errmax(3) = field_error(ng,lo,n,nh_v,nh_g,is_c(:,3),cbcvel(:,:,3),bcvel(:,:,3),dl,zc_g,zf_g,w)
    call MPI_ALLREDUCE(MPI_IN_PLACE,errmax,3,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(myid == 0.and.any(errmax(:) > 0.)) &
      print*, 'ERROR: wrong velocity halo values; maximum relative errors of u, v, w = ', errmax(:)
    passed = passed.and.all(errmax(:) <= 0.)
    deallocate(u,v,w)
    !
    ! pressure
    !
    allocate(p(1-nh_p:n(1)+nh_p,1-nh_p:n(2)+nh_p,1-nh_p:n(3)+nh_p))
    call set_field(n,nh_p,lo,p)
    !$acc enter data copyin(p)
    call boundp(cbcpre,n,nh_p,nh_g,bcpre,nb,is_bound,dl,dzc,p)
    !$acc wait
    !$acc exit data copyout(p)
    errmax(1) = field_error(ng,lo,n,nh_p,nh_g,[.true.,.true.,.true.],cbcpre,bcpre,dl,zc_g,zf_g,p)
    call MPI_ALLREDUCE(MPI_IN_PLACE,errmax(1),1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(myid == 0.and.errmax(1) > 0.) &
      print*, 'ERROR: wrong pressure halo values; maximum relative error = ', errmax(1)
    passed = passed.and.errmax(1) <= 0.
    deallocate(p)
    !
    ! scalar
    !
    if(present(cbcscal).and.present(bcscal)) then
      allocate(p(1-nh_s:n(1)+nh_s,1-nh_s:n(2)+nh_s,1-nh_s:n(3)+nh_s))
      call set_field(n,nh_s,lo,p)
      !$acc enter data copyin(p)
      call boundp(cbcscal,n,nh_s,nh_g,bcscal,nb,is_bound,dl,dzc,p)
      !$acc wait
      !$acc exit data copyout(p)
      errmax(1) = field_error(ng,lo,n,nh_s,nh_g,[.true.,.true.,.true.],cbcscal,bcscal,dl,zc_g,zf_g,p)
      call MPI_ALLREDUCE(MPI_IN_PLACE,errmax(1),1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid == 0.and.errmax(1) > 0.) &
        print*, 'ERROR: wrong scalar halo values; maximum relative error = ', errmax(1)
      passed = passed.and.errmax(1) <= 0.
      deallocate(p)
    end if
    if(.not.passed) then
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      error stop
    end if
    if(myid == 0) print*, '*** Halo check passed for widths (velocity, pressure, scalars): ', nh_v,nh_p,nh_s, ' ***'
  end subroutine test_sanity_halo
  !
  pure real(rp) function fval(i,j,k)
    !
    ! test field: exact and distinct for less than 64 points per direction, halos included
    !
    implicit none
    integer, intent(in) :: i,j,k
    fval = i + 64._rp*j + 4096._rp*k
  end function fval
  !
  subroutine set_field(n,nh,lo,p)
    implicit none
    integer , intent(in), dimension(3) :: n,lo
    integer , intent(in)               :: nh
    real(rp), intent(out), dimension(1-nh:,1-nh:,1-nh:) :: p
    integer :: i,j,k
    p(:,:,:) = -huge(1._rp)/8 ! stale halo values
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          p(i,j,k) = fval(i+lo(1)-1,j+lo(2)-1,k+lo(3)-1)
        end do
      end do
    end do
  end subroutine set_field
  !
  function field_error(ng,lo,n,nh,nh_g,is_c,cbc,bc,dl,zc_g,zf_g,p) result(errmax)
    !
    ! maximum relative error over all points, halos included (zero if within round-off)
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,n
    integer , intent(in)               :: nh,nh_g
    logical , intent(in), dimension(3) :: is_c
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp), intent(in), dimension(0:1,3) :: bc
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(1-nh_g:) :: zc_g,zf_g
    real(rp), intent(in), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: errmax,val,err
    integer :: i,j,k
    errmax = 0.
    do k=1-nh,n(3)+nh
      do j=1-nh,n(2)+nh
        do i=1-nh,n(1)+nh
          val = expected_value([i,j,k]+lo(:)-1,ng,nh_g,is_c,cbc,bc,dl,zc_g,zf_g)
          err = abs(p(i,j,k)-val)/max(1._rp,abs(val))
          if(err > small) errmax = max(errmax,err)
        end do
      end do
    end do
  end function field_error
  !
  recursive function expected_value(ii,ng,nh_g,is_c,cbc,bc,dl,zc_g,zf_g) result(val)
    !
    ! value of the test field at global index `ii`, inside or outside the domain, after the halo
    ! exchange and the boundary conditions, which are imposed along x, y, and z, in this order:
    ! outside the domain, it is the value at the periodic image or, for a boundary, the value at
    ! the mirror point transformed by the boundary condition (cf. `set_bc` under `bound.f90`),
    ! starting from the last direction of the sequence along which `ii` is outside
    !
    implicit none
    integer , intent(in), dimension(3) :: ii,ng
    integer , intent(in)               :: nh_g
    logical , intent(in), dimension(3) :: is_c
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp), intent(in), dimension(0:1,3) :: bc
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(1-nh_g:) :: zc_g,zf_g
    real(rp) :: val
    integer , dimension(3) :: im
    integer :: idir,i,ib,dh
    real(rp) :: dr
    do idir=3,1,-1
      i = ii(idir)
      im(:) = ii(:)
      if(cbc(0,idir)//cbc(1,idir) == 'PP') then
        if(i >= 1 .and. i <= ng(idir)) cycle
        im(idir) = modulo(i-1,ng(idir))+1
        val = expected_value(im,ng,nh_g,is_c,cbc,bc,dl,zc_g,zf_g)
        return
      end if
      if(is_c(idir)) then
        if(     i < 1       ) then
          ib = 0; dh = -i         ; im(idir) = 1+dh
        else if(i > ng(idir)) then
          ib = 1; dh = i-ng(idir)-1; im(idir) = ng(idir)-dh
        else
          cycle
        end if
      else ! the boundary node is also set by the boundary condition
        if(     i <= 0       ) then
          ib = 0; dh = -i         ; im(idir) = max(dh,1)
        else if(i >= ng(idir)) then
          ib = 1; dh = i-ng(idir) ; im(idir) = ng(idir)-max(dh,1)
        else
          cycle
        end if
      end if
      if(idir == 3) then
        if(is_c(idir)) then
          dr = abs(zc_g(im(idir))-zc_g(i))
        else
          dr = abs(zf_g(im(idir))-zf_g(i))
        end if
      else
        dr = abs(im(idir)-i)*dl(idir)
      end if
      select case(cbc(ib,idir))
      case('D')
        if(.not.is_c(idir) .and. dh == 0) then
          val = bc(ib,idir)
        else
          val = 2.*bc(ib,idir) - expected_value(im,ng,nh_g,is_c,cbc,bc,dl,zc_g,zf_g)
        end if
      case('N')
        val = expected_value(im,ng,nh_g,is_c,cbc,bc,dl,zc_g,zf_g) + (2*ib-1)*dr*bc(ib,idir)
      case default
        val = huge(1._rp)
      end select
      return
    end do
    val = fval(ii(1),ii(2),ii(3))
  end function expected_value
end module mod_sanity
