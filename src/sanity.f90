module mod_sanity
  use iso_c_binding, only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound     , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv    , only: chkdiv
  use mod_common_mpi, only: myid,ierr,n_z
  use mod_correc    , only: correc
  use mod_debug     , only: chk_helmholtz
  use mod_fft       , only: fftend
  use mod_fillps    , only: fillps
  use mod_initflow  , only: add_noise
  use mod_initmpi   , only: initmpi
  use mod_initsolver, only: initsolver
  use mod_param     , only: small
  use mod_solver    , only: solver
  use mod_types
  implicit none
  private
  public test_sanity
  contains
  subroutine test_sanity(ng,n,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_outflow,is_forced, &
                         dli,dzci_g,dzfi_g,dzci,dzfi)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,n
    integer , intent(in), dimension(3) :: dims
    logical , intent(in), dimension(3) :: stop_type
    character(len=1), intent(in), dimension(0:1,3,3)  :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)    :: cbcpre
    real(rp), intent(in), dimension(0:1,3,3)          :: bcvel
    real(rp), intent(in), dimension(0:1,3)            :: bcpre
    logical , intent(in), dimension(0:1,3)            :: is_outflow
    logical , intent(in), dimension(3)                :: is_forced
    real(rp), intent(in), dimension(3)                :: dli
    real(rp), intent(in), dimension(0:)               :: dzci_g,dzfi_g,dzci,dzfi
    logical :: passed
    !
    call chk_dims(ng,dims,passed);                 if(.not.passed) call abortit
    call chk_stop_type(stop_type,passed);          if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed); if(.not.passed) call abortit
    call chk_outflow(cbcpre,is_outflow,passed);    if(.not.passed) call abortit
    call chk_forcing(cbcpre,is_forced  ,passed);   if(.not.passed) call abortit 
    !call chk_solvers(ng,n,dli,dzci_g,dzfi_g,dzci,dzfi,cbcvel,cbcpre,bcvel,bcpre,is_outflow,passed)
    if(.not.passed) call abortit
    return
  end subroutine test_sanity
  !
  subroutine chk_stop_type(stop_type,passed)
  implicit none
  logical, intent(in), dimension(3) :: stop_type
  logical, intent(out) :: passed
  passed = .true.
  if(.not.any(stop_type(:))) then
    if(myid.eq.0) print*, 'ERROR: stopping criterion not chosen.'
    passed = .false.
  endif
  return 
  end subroutine chk_stop_type
  !
  subroutine chk_dims(ng,dims,passed)
    implicit none
    integer, intent(in), dimension(3) :: ng
    integer, intent(in), dimension(3) :: dims
    logical, intent(out) :: passed
    logical :: passed_loc
    passed = .true.
    passed_loc = all(mod(ng(:),2).eq.0)
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: itot, jtot and ktot should be even.'
    passed = passed.and.passed_loc
    passed_loc = all(mod(ng(1:3),dims(1:3)).eq.0)
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: itot, jtot and ktot should be divisable by dims(1), dims(2) and dims(3), respectively.'
    passed = passed.and.passed_loc
    passed_loc = (mod(ng(2),dims(1)).eq.0).and.(mod(ng(3),dims(2)).eq.0)
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: jtot should be divisable by both dims(1) and dims(2), and &
                     &ktot should be divisable by dims(2)'
    passed = passed.and.passed_loc
    return
  end subroutine chk_dims
  !
  subroutine chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed)
  implicit none
  character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
  character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
  real(rp)        , intent(in), dimension(0:1,3,3) :: bcvel
  real(rp)        , intent(in), dimension(0:1,3  ) :: bcpre
  logical         , intent(out) :: passed
  character(len=2) :: bc01v,bc01p
  integer :: ivel,idir
  logical :: passed_loc
  passed = .true.
  !
  ! check validity of pressure and velocity BCs
  !
  passed_loc = .true.
  do ivel = 1,3
    do idir=1,3
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      passed_loc = passed_loc.and.( (bc01v.eq.'PP').or. &
                                    (bc01v.eq.'ND').or. &
                                    (bc01v.eq.'DN').or. &
                                    (bc01v.eq.'NN').or. &
                                    (bc01v.eq.'DD') )
    enddo
  enddo
  if(myid.eq.0.and.(.not.passed_loc)) print*, 'ERROR: velocity BCs not valid.'
  passed = passed.and.passed_loc
  !
  passed_loc = .true.
  do idir=1,3
    bc01p = cbcpre(0,idir)//cbcpre(1,idir)
    passed_loc = passed_loc.and.( (bc01p.eq.'PP').or. &
                                  (bc01p.eq.'ND').or. &
                                  (bc01p.eq.'DN').or. &
                                  (bc01p.eq.'NN').or. &
                                  (bc01p.eq.'DD') )
  enddo
  if(myid.eq.0.and.(.not.passed_loc)) print*, 'ERROR: pressure BCs not valid.' 
  passed = passed.and.passed_loc
  !
  passed_loc = .true.
  do idir=1,3
    ivel = idir
    bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
    bc01p = cbcpre(0,idir)//cbcpre(1,idir)
    passed_loc = passed_loc.and.( (bc01v.eq.'PP'.and.bc01p.eq.'PP').or. &
                                  (bc01v.eq.'ND'.and.bc01p.eq.'DN').or. &
                                  (bc01v.eq.'DN'.and.bc01p.eq.'ND').or. &
                                  (bc01v.eq.'DD'.and.bc01p.eq.'NN').or. &
                                  (bc01v.eq.'NN'.and.bc01p.eq.'DD') )
  enddo
  if(myid.eq.0.and.(.not.passed_loc)) print*, 'ERROR: velocity and pressure BCs not compatible.'
  passed = passed.and.passed_loc
  !
  passed_loc = .true.
  do idir=1,2
    passed_loc = passed_loc.and.((bcpre(0,idir).eq.0.).and.(bcpre(1,idir).eq.0.))
  enddo
  if(myid.eq.0.and.(.not.passed_loc)) &
    print*, 'ERROR: pressure BCs in directions x and y must be homogeneous (value = 0.).'
  passed = passed.and.passed_loc
#ifdef IMPDIFF
  passed_loc = .true.
  do ivel = 1,3
    do idir=1,2
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      passed_loc = passed_loc.and.(bc01v.ne.'NN')
    enddo
  enddo
  if(myid.eq.0.and.(.not.passed_loc)) &
    print*, 'ERROR: Neumann-Neumann velocity BCs with implicit diffusion currently not supported in x and y; only in z.'
  passed = passed.and.passed_loc
  !
  passed_loc = .true.
  do ivel = 1,3
    do idir=1,2
      passed_loc = passed_loc.and.((bcvel(0,idir,ivel).eq.0.).and.(bcvel(1,idir,ivel).eq.0.))
    enddo
  enddo
  if(myid.eq.0.and.(.not.passed_loc)) &
    print*, 'ERROR: velocity BCs with implicit diffusion in directions x and y must be homogeneous (value = 0.).'
  passed = passed.and.passed_loc
#endif
  return 
  end subroutine chk_bc
  !
  subroutine chk_outflow(cbcpre,is_outflow,passed)
  implicit none
  logical         , intent(in), dimension(0:1,3  ) :: is_outflow
  character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
  logical         , intent(out) :: passed
  integer :: idir,ibound
  passed = .true.
  !
  ! 1) check for compatibility between pressure BCs and outflow BC
  !
  do idir=1,3
    do ibound = 0,1
      passed = passed.and. &
               (cbcpre(ibound,idir).eq.'D'.and.(is_outflow(ibound,idir))) .or. &
               (.not.is_outflow(ibound,idir))
    enddo
  enddo
  if(myid.eq.0.and.(.not.passed)) &
    print*, 'ERROR: Dirichlet pressure BC should be an outflow direction; check the BC or is_outflow in dns.in.'
  return 
  end subroutine chk_outflow
  !
  subroutine chk_forcing(cbcpre,is_forced,passed)
  implicit none
  character(len=1), intent(in), dimension(0:1,3) :: cbcpre
  logical         , intent(in), dimension(3) :: is_forced
  logical         , intent(out) :: passed
  integer :: idir
  passed = .true.
  !
  ! 1) check for compatibility between pressure BCs and forcing BC
  !
  do idir=1,3
    if(is_forced(idir)) then
      passed = passed.and.(cbcpre(0,idir)//cbcpre(1,idir).eq.'PP')
    endif
  enddo
  if(myid.eq.0.and.(.not.passed)) &
  print*, 'ERROR: Flow cannot be forced in a non-periodic direction; check the BCs and is_forced in dns.in.'
  return 
  end subroutine chk_forcing
  !
  subroutine chk_solvers(ng,n,dli,dzci_g,dzfi_g,dzci,dzfi,cbcvel,cbcpre,bcvel,bcpre,is_outflow,passed)
  implicit none
  integer , intent(in), dimension(3) :: ng,n
  real(rp), intent(in), dimension(3) :: dli
  real(rp), intent(in), dimension(0:) :: dzci_g,dzfi_g,dzci,dzfi
  character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
  character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
  real(rp), intent(in), dimension(0:1,3,3)          :: bcvel
  real(rp), intent(in), dimension(0:1,3)            :: bcpre
  logical , intent(in), dimension(0:1,3)            :: is_outflow
  logical , intent(out) :: passed
  real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: u,v,w,p,up,vp,wp
  type(C_PTR), dimension(2,2) :: arrplan
  real(rp), dimension(n_z(1),n_z(2)) :: lambdaxy
  real(rp) :: normfft
  real(rp), dimension(n_z(3)) :: a,b,c,bb
  real(rp), dimension(n(2),n(3),0:1) :: rhsbx
  real(rp), dimension(n(1),n(3),0:1) :: rhsby
  real(rp), dimension(n(1),n(2),0:1) :: rhsbz
  logical , dimension(0:1,3)            :: no_outflow
  real(rp), dimension(3) :: dl
  real(rp), dimension(0:n(3)+1) :: dzc,dzf
  real(rp) :: dt,dti,alpha
  real(rp) :: divtot,divmax,resmax
  logical :: passed_loc
  passed = .true.
  !
  ! initialize velocity below with some random noise
  !
  up(:,:,:) = 0.
  vp(:,:,:) = 0.
  wp(:,:,:) = 0.
  call add_noise(n,123,.5_rp,up(1:n(1),1:n(2),1:n(3)))
  call add_noise(n,456,.5_rp,vp(1:n(1),1:n(2),1:n(3)))
  call add_noise(n,789,.5_rp,wp(1:n(1),1:n(2),1:n(3)))
  !
  ! test pressure correction
  !
  call initsolver(n,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:),lambdaxy,(/'c','c','c'/),a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
  dl  = dli**(-1)
  dzc = dzci**(-1)
  dzf = dzfi**(-1)
  dt  = acos(-1.) ! value is irrelevant
  dti = dt**(-1)
  no_outflow(:,:) = .false.
  call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp)
  call fillps(n,dli,dzfi,dti,up,vp,wp,p)
  call updt_rhs_b((/'c','c','c'/),cbcpre,n,rhsbx,rhsby,rhsbz,p)
  call solver(n,arrplan,normfft,lambdaxy,a,b,c,cbcpre(:,3),(/'c','c','c'/),p)
  call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
  call correc(n,dli,dzci,dt,p,up,vp,wp,u,v,w)
  call bounduvw(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
  call chkdiv(n,dli,dzfi,u,v,w,divtot,divmax)
  passed_loc = divmax.lt.small
  if(myid.eq.0.and.(.not.passed_loc)) &
  print*, 'ERROR: Pressure correction: Divergence is too large.'
  passed = passed.and.passed_loc
  call fftend(arrplan)
#ifdef IMPDIFF
  alpha = acos(-1.) ! irrelevant
  up(:,:,:) = 0.
  vp(:,:,:) = 0.
  wp(:,:,:) = 0.
  call add_noise(n,123,.5_rp,up(1:n(1),1:n(2),1:n(3)))
  call add_noise(n,456,.5_rp,vp(1:n(1),1:n(2),1:n(3)))
  call add_noise(n,789,.5_rp,wp(1:n(1),1:n(2),1:n(3)))
  call initsolver(n,dli,dzci_g,dzfi_g,cbcvel(:,:,1),bcvel(:,:,1),lambdaxy,(/'f','c','c'/),a,b,c,arrplan,normfft, &
                  rhsbx,rhsby,rhsbz)
  call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp)
  up(:,:,:) = up(:,:,:)*alpha
  u( :,:,:) = up(:,:,:)
  bb(:) = b(:) + alpha
  call updt_rhs_b((/'f','c','c'/),cbcvel(:,:,1),n,rhsbx,rhsby,rhsbz,up)
  call solver(n,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,3,1),(/'f','c','c'/),up)
  call fftend(arrplan)
  call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp) ! actually we are only interested in boundary condition in up
  call chk_helmholtz(n,dli,dzci,dzfi,alpha,u,up,cbcvel(:,:,1),(/'f','c','c'/),resmax)
  passed_loc = resmax.lt.small
  if(myid.eq.0.and.(.not.passed_loc)) &
  print*, 'ERROR: wrong solution of Helmholtz equation in x direction.'
  passed = passed.and.passed_loc
  !
  call initsolver(n,dli,dzci_g,dzfi_g,cbcvel(:,:,2),bcvel(:,:,2),lambdaxy,(/'c','f','c'/),a,b,c,arrplan,normfft, &
                  rhsbx,rhsby,rhsbz)
  call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp)
  vp(:,:,:) = vp(:,:,:)*alpha
  v( :,:,:) = vp(:,:,:)
  bb(:) = b(:) + alpha
  call updt_rhs_b((/'c','f','c'/),cbcvel(:,:,2),n,rhsbx,rhsby,rhsbz,vp)
  call solver(n,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,3,2),(/'c','f','c'/),vp)
  call fftend(arrplan)
  call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp) ! actually we are only interested in boundary condition in up
  call chk_helmholtz(n,dli,dzci,dzfi,alpha,v,vp,cbcvel(:,:,2),(/'c','f','c'/),resmax)
  passed_loc = resmax.lt.small
  if(myid.eq.0.and.(.not.passed_loc)) &
  print*, 'ERROR: wrong solution of Helmholtz equation in y direction.'
  passed = passed.and.passed_loc
  !
  call initsolver(n,dli,dzci_g,dzfi_g,cbcvel(:,:,3),bcvel(:,:,3),lambdaxy,(/'c','c','f'/),a,b,c,arrplan,normfft, &
                  rhsbx,rhsby,rhsbz)
  call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp)
  wp(:,:,:) = wp(:,:,:)*alpha
  w( :,:,:) = wp(:,:,:)
  bb(:) = b(:) + alpha
  call updt_rhs_b((/'c','c','f'/),cbcvel(:,:,3),n,rhsbx,rhsby,rhsbz,wp)
  call solver(n,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,3,3),(/'c','c','f'/),wp)
  call fftend(arrplan)
  call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp) ! actually we are only interested in boundary condition in up
  call chk_helmholtz(n,dli,dzci,dzfi,alpha,w,wp,cbcvel(:,:,3),(/'c','c','f'/),resmax)
  passed_loc = resmax.lt.small
  if(myid.eq.0.and.(.not.passed_loc)) &
  print*, 'ERROR: wrong solution of Helmholtz equation in z direction.'
  passed = passed.and.passed_loc
#endif
  return
  end subroutine chk_solvers
  !
  subroutine abortit
      implicit none
      if(myid.eq.0) print*, ''
      if(myid.eq.0) print*, '*** Simulation aborted due to errors in the input file ***'
      if(myid.eq.0) print*, '    check dns.in'
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      call exit
      return
  end subroutine abortit
end module mod_sanity
