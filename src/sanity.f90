module mod_sanity
!use mod_common_mpi, only:myid
implicit none
contains
  subroutine test_sanity

    if(passed) then 
    call chk_n(ng,dims)
    if(.not.passed) exit
    call chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed)
    if(.not.passed) exit
    endif

    if(.not.passed) then
      if(myid.eq.0) print*, '*** Simulation aborted. Check bc.h90 and setup.h90 and restart ***'
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      call exit
    endif
  end subroutine test_sanity
  subroutine chk_n(ng,dims)
    implicit none
    if(myid.eq.0.and.(.not.passed)) print*, 'itot should be even and divisible by dims(1); check setup.h90.'
    if(myid.eq.0.and.(.not.passed)) print*, 'jtot should be even and divisible by dims(2); check setup.h90.'
    if(myid.eq.0.and.(.not.passed)) print*, 'ktot should be even;check setup.h90.'
  end subroutine chk_n
  subroutine chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed)
  implicit none
  character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
  character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
  real(8)         , intent(in), dimension(0:1,3,3) :: bcvel
  real(8)         , intent(in), dimension(0:1,3  ) :: bcpre
  logical         , intent(out) :: passed
  character(len=2) :: bc01v,bc01p
  integer :: ivel,idir,ibound
  passed = .true.
  !
  ! check validity of pressure and velocity BCs
  !
  do ivel = 1,3
    do idir=1,3
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      passed = passed.and.(bc01v.eq.'PP'.or. &
                           bc01v.eq.'ND'.or. &
                           bc01v.eq.'DN'.or. &
                           bc01v.eq.'NN'.or. &
                           bc01v.eq.'DD')
    enddo
  enddo
  if(myid.eq.0.and.(.not.passed)) print*, 'Velocity BCs not valid.'
  do idir=1,3
    bc01p = cbcpre(0,idir)//cbcpre(1,idir)
    passed = passed.and.(bc01p.eq.'PP'.or. &
                         bc01p.eq.'ND'.or. &
                         bc01p.eq.'DN'.or. &
                         bc01p.eq.'NN'.or. &
                         bc01p.eq.'DD')
  enddo
  if(myid.eq.0.and.(.not.passed)) print*, 'Pressure BCs not valid.' 
  do idir=1,3
    ivel = idir
    bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
    bc01p = cbcpre(0,idir)//cbcpre(1,idir)
    passed = passed.and.((bc01v.eq.'PP'.and.bc01p.eq.'PP').or. &
                         (bc01v.eq.'ND'.and.bc01p.eq.'DN').or. &
                         (bc01v.eq.'DN'.and.bc01p.eq.'ND').or. &
                         (bc01v.eq.'NN'.and.bc01p.eq.'DD'))
  enddo
  if(myid.eq.0.and.(.not.passed)) print*, 'Velocity and Pressure BCs not compatible.'
  do idir=1,2
    passed = passed.and.((bcpre(0,idir).eq.0.d0).and.(bcpre(1,idir).eq.0.d0))
  enddo
  if(myid.eq.0.and.(.not.passed)) &
  print*, 'Pressure BCs in directions x and y must be homogeneous (value = 0.d0).'
#ifdef IMPDIFF
  do ivel = 1,3
    do idir=1,2
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      passed = passed.and.(bc01v.ne.'NN'.or. &
    enddo
  enddo
  if(myid.eq.0.and.(.not.passed)) &
  print*, 'Neumann-Neumann velocity BCs with implicit diffusion currently not supported in x and y; only in z.'
  do ivel = 1,3
    do idir=1,2
    passed = passed.and.((bcvel(0,idir,ivel).eq.0.d0).and.(bcvel(1,idir,ivel).eq.0.d0))
    enddo
  enddo
  if(myid.eq.0.and.(.not.passed)) &
  print*, 'Velocity BCs with implicit diffusion in directions x and y must be homogeneous (value = 0.d0).'
#endif
  return 
  end subroutine chk_bc
  subroutine chk_outflow(is_outflow,cbcpre,passed)
  implicit none
  character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
  logical         , intent(in), dimension(0:1,3  ) :: is_outflow
  logical         , intent(out) :: passed
  integer :: ivel,idir,ibound
  passed = .true.
  !
  ! 1) check for compatibility between pressure BCs and outflow BC
  !
  do idir=1,3
    do ibound = 0,1
      passed = passed.and.(cbcpre(ibound,idir).eq.'D'.and.(is_outflow(ibound,idir)))
    enddo
  enddo
  if(myid.eq.0.and.(.not.passed)) &
  print*, 'Dirichlet pressure BC should be an outflow direction; check the BC or is_outflow in bc.h90.'
  return 
  end subroutine chk_outflow
  subroutine chk_forcing(is_force,cbcpre,passed)
  implicit none
  character(len=1), intent(in), dimension(0:1,3) :: cbcpre
  logical         , intent(in), dimension(3) :: is_force
  logical         , intent(out) :: passed
  integer :: ivel,idir,ibound
  passed = .true.
  !
  ! 1) check for compatibility between pressure BCs and forcing BC
  !
  do idir=1,3
    if(is_force(idir)) then
      passed = passed.and.(cbcpre(0,idir)//cbcpre(1,idir).eq.'PP')
    endif
  enddo
  if(myid.eq.0.and.(.not.passed)) &
  print*, 'Flow cannot be forced in a non-periodic direction; check the BCs and is_force in bc.h90.'
  return 
  end subroutine chk_forcing
  subroutine chk_solvers(n,dli,dzci,dzfi,cbcvel,cbcrce,bcvel,bcpre,is_outflow,passed)
  implicit none
  integer, intent(in), dimension(3) :: n
  integer, intent(in), dimension(3) :: n
  real(8), intent(in), dimension(3) :: dli
  real(8), intent(in), dimension(0:) :: dzci,dzfi
  real(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: u,v,w,p,up,vp,wp
  type(C_PTR), dimension(2,2) :: arrplan
  real(8), dimension(n(1),n(2)) :: lambdaxy
  real(8) :: normfftp
  type rhs_bound
    real(8), dimension(n(2),n(3),0:1) :: x
    real(8), dimension(n(1),n(3),0:1) :: y
    real(8), dimension(n(1),n(2),0:1) :: z
  end type rhs_bound 
  type(rhs_bound) :: rhsb
  real(8), dimension(3) :: dl
  real(8) :: dt
  !
  ! initialize velocity below with some random noise
  !
  if(isnoise) then
      call add_noise(n,123,.50d0,up(1:n(1),1:n(2),1:n(3)))
      call add_noise(n,456,.50d0,vp(1:n(1),1:n(2),1:n(3)))
      call add_noise(n,789,.50d0,wp(1:n(1),1:n(2),1:n(3)))
  endif
  !
  ! test pressure correction
  !
  dl  = dli**(-1)
  dzc = dzci**(-1)
  dzf = dzfi**(-1)
  dt  = acos(-1.d0) ! value is irrelevant
  call initsolver(n,dli,dzci,dzfi,cbcpre,bcpre(:,:),lambdaxy,(/'c','c','c'/),a,b,c,arrplan,normfft,rhsb%x,rhsb%y,rhsb%z)
  call fillps(n,dli,dzfi,dtrki,up,vp,wp,p)
  call updt_rhs_b((/'c','c','c'/),cbcpre,n,rhsbp%x,rhsbp%y,rhsbp%z,p(1:n(1),1:n(2),1:n(3)))
  call solver(n,arrplan,normfft,lambdaxy,a,b,c,cbcpre(:,3),(/'c','c','c'/),p(1:n(1),1:n(2),1:n(3)))
  call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
  call correc(n,dli,dzci,dtrk,p,up,vp,wp,u,v,w)
  call bounduvw(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
  call chkdiv(n,dli,dzfi,u,v,w,divtot,divmax)
  passed = passed.and.(divmax.gt.small)
  if(myid.eq.0.and.(.not.passed)) &
  print*, 'Error in pressure correction. Divergence is too large.'
  call fftend(arrplan)
#ifdef IMPDIFF
  alpha = acos(-1.d0) ! irrelevant
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,1),bcvel(:,:,1),lambdaxy,(/'f','c','c'/),a,b,c,arrplan,normfft, &
                  rhsb%x,rhsb%y,rhsb%z)
  up(:,:,:) = up(:,:,:)*alpha
  u( :,:,:) = up(:,:,:)
  b(:) = b(:) + alpha
  call updt_rhs_b((/'f','c','c'/),cbcvel(:,:,1),n,rhsb%x,rhsb%y,rhsb%z,up(1:n(1),1:n(2),1:n(3)))
  call solver(n,arrplan,normfft,lambdaxy,a,b,c,cbcvel(:,3,1),(/'f','c','c'/),up(1:n(1),1:n(2),1:n(3)))
  call fftend(arrplan)
  call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp) ! actually we are only interested in boundary condition in up
  call chkhelmholtz(n,dxi,dyi,dzci,dzfi,alpha,u,up,cbcvel(:,:,1),'c',diffmax)
  passed = passed.and.(diffmax.gt.small)
  if(myid.eq.0.and.(.not.passed)) &
  print*, 'Error in solving the Helmholtz equation in x direction.'
  !
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,2),bcvel(:,:,2),lambdaxy,(/'c','f','c'/),a,b,c,arrplan,normfft, &
                  rhsb%x,rhsb%y,rhsb%z)
  vp(:,:,:) = vp(:,:,:)*alpha
  v( :,:,:) = vp(:,:,:)
  b(:) = b(:) + alpha
  call updt_rhs_b((/'c','f','c'/),cbcvel(:,:,2),n,rhsb%x,rhsb%y,rhsb%z,vp(1:n(1),1:n(2),1:n(3)))
  call solver(n,arrplan,normfft,lambdaxy,a,b,c,cbcvel(:,3,2),(/'f','c','c'/),vp(1:n(1),1:n(2),1:n(3)))
  call fftend(arrplan)
  call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp) ! actually we are only interested in boundary condition in up
  call chkhelmholtz(n,dxi,dyi,dzci,dzfi,alpha,u,up,cbcvel(:,:,2),'c',diffmax)
  passed = passed.and.(diffmax.gt.small)
  if(myid.eq.0.and.(.not.passed)) &
  print*, 'Error in solving the Helmholtz equation in y direction.'
  !
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,3),bcvel(:,:,3),lambdaxy,(/'c','c','f'/),a,b,c,arrplan,normfft, &
                  rhsb%x,rhsb%y,rhsb%z)
  wp(:,:,:) = wp(:,:,:)*alpha
  v( :,:,:) = wp(:,:,:)
  b(:) = b(:) + alpha
  call updt_rhs_b((/'c','c','f'/),cbcvel(:,:,3),n,rhsb%x,rhsb%y,rhsb%z,wp(1:n(1),1:n(2),1:n(3)))
  call solver(n,arrplan,normfft,lambdaxy,a,b,c,cbcvel(:,3,3),(/'f','c','c'/),wp(1:n(1),1:n(2),1:n(3)))
  call fftend(arrplan)
  call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,wp,wp) ! actually we are only interested in boundary condition in up
  call chkhelmholtz(n,dxi,dyi,dzci,dzfi,alpha,u,up,cbcvel(:,:,3),'f',diffmax)
  passed = passed.and.(diffmax.gt.small)
  if(myid.eq.0.and.(.not.passed)) &
  print*, 'Error in solving the Helmholtz equation in y direction.'
#endif
  return
  end subroutine chk_solvers
end module mod_sanity
