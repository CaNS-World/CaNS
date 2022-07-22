!  subroutine initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci,dzfi,cbc,bc,lambdaxy,c_or_f,a,b,c,arrplan,normfft, &
!                        rhsbx,rhsby,rhsbz)
    !
    ! initializes the Poisson/Helmholtz solver
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,n_x_fft,n_y_fft,lo_z,hi_z
    real(rp), intent(in), dimension(3 ) :: dli
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp)        , intent(in), dimension(0:1,3) :: bc
    MYREAL, intent(out), dimension(lo_z(1):,lo_z(2):) :: lambdaxy
    character(len=1), intent(in), dimension(3) :: c_or_f
    MYREAL, intent(out), dimension(:) :: a,b,c
    type(C_PTR), intent(out), dimension(2,2) :: arrplan
    real(rp), intent(out), dimension(:,:,0:) :: rhsbx
    real(rp), intent(out), dimension(:,:,0:) :: rhsby
    real(rp), intent(out), dimension(:,:,0:) :: rhsbz
    real(rp), intent(out) :: normfft
    real(rp), dimension(3)        :: dl
    real(rp), dimension(0:ng(3)+1) :: dzc,dzf
    integer :: i,j
    real(gp), dimension(ng(1))      :: lambdax
    real(gp), dimension(ng(2))      :: lambday
    !
    ! generating eigenvalues consistent with the BCs
    !
    call eigenvalues(ng(1),cbc(:,1),c_or_f(1),lambdax)
    lambdax(:) = lambdax(:)*dli(1)**2
    call eigenvalues(ng(2),cbc(:,2),c_or_f(2),lambday)
    lambday(:) = lambday(:)*dli(2)**2
    !
    ! add eigenvalues
    !
    do j=lo_z(2),hi_z(2)
      do i=lo_z(1),hi_z(1)
        lambdaxy(i,j) = lambdax(i)+lambday(j)
      end do
    end do
    !
    ! compute coefficients for tridiagonal solver
    !
    call tridmatrix(cbc(:,3),ng(3),dli(3),dzci,dzfi,c_or_f(3),a,b,c)
    !
    ! compute values to be added to the right hand side
    !
    dl(:)  = dli( :)**(-1)
    dzc(:) = dzci(:)**(-1)
    dzf(:) = dzfi(:)**(-1)
    call bc_rhs(cbc(:,1),bc(:,1),[dl(1) ,dl(1)      ],[dl(1) ,dl(1)    ],c_or_f(1),rhsbx)
    call bc_rhs(cbc(:,2),bc(:,2),[dl(2) ,dl(2)      ],[dl(2) ,dl(2)    ],c_or_f(2),rhsby)
    if(     c_or_f(3) == 'c') then
      call bc_rhs(cbc(:,3),bc(:,3),[dzc(0),dzc(ng(3)  )],[dzf(1),dzf(ng(3))],c_or_f(3),rhsbz)
    else if(c_or_f(3) == 'f') then
      call bc_rhs(cbc(:,3),bc(:,3),[dzc(1),dzc(ng(3)-1)],[dzf(1),dzf(ng(3))],c_or_f(3),rhsbz)
    end if
    !
    ! prepare ffts
    !
    call fftini(n_x_fft,n_y_fft,cbc(:,1:2),c_or_f(1:2),arrplan,normfft)
!  end subroutine initsolver
