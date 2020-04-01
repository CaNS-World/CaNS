module mod_debug
  use mpi
  use mod_common_mpi, only: myid,ierr!,ijk_min
  use mod_param     , only: dims
  use mod_types
  implicit none
  private
  public chkmean,chk_helmholtz
  contains
  subroutine chkmean(n,dzlzi,p,mean)
    !
    ! compute the mean value of an observable over the entire domain
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: dzlzi
    real(rp), intent(in), dimension(0:,0:,0:) :: p
    real(rp), intent(out) :: mean
    integer :: i,j,k
    mean = 0.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,dzlzi) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:mean)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          mean = mean + p(i,j,k)*dzlzi(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean = mean/(1.*n(1)*dims(1)*n(2)*dims(2))
    return
  end subroutine chkmean
  !
  subroutine chk_helmholtz(n,dli,dzci,dzfi,alpha,fp,fpp,bc,c_or_f,diffmax)
    !
    ! this subroutine checks if the implementation of implicit diffusion is
    ! correct
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(2) :: dli
    real(rp), intent(in) :: alpha
    real(rp), intent(in), dimension(0:) :: dzfi,dzci
    real(rp), intent(in), dimension(0:,0:,0:) :: fp,fpp
    character(len=1), intent(in), dimension(0:1,3) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(out) :: diffmax
    real(rp) :: val
    integer :: i,j,k,im,ip,jm,jp,km,kp
    integer :: idir
    integer, dimension(3) :: q
    !integer :: ii,jj,kk
    q(:) = 0
    do idir = 1,3
      if(bc(1,idir).ne.'P'.and.c_or_f(idir).eq.'f') q(idir) = 1
    enddo
    select case(c_or_f(3))
    !
    ! need to compute the maximum difference!
    !
    case('c')
      diffmax = 0.
      do k=1,n(3)-q(3)
        kp = k + 1
        km = k - 1
        do j=1,n(2)-q(2)
          jp = j + 1
          jm = j - 1
          do i=1,n(1)-q(1)
            ip = i + 1
            im = i - 1
            val =  fpp(i,j,k)+(1./alpha)*( &
                  (fpp(ip,j,k)-2.*fpp(i,j,k)+fpp(im,j,k))*(dli(1)**2) + &
                  (fpp(i,jp,k)-2.*fpp(i,j,k)+fpp(i,jm,k))*(dli(2)**2) + &
                 ((fpp(i,j,kp)-fpp(i,j,k ))*dzci(k ) - &
                  (fpp(i,j,k )-fpp(i,j,km))*dzci(km))*dzfi(k) )
            val = val*alpha
            diffmax = max(diffmax,abs(val-fp(i,j,k)))
            !ii = ijk_min(1)+i
            !jj = ijk_min(2)+j
            !kk = ijk_min(3)+k
            !if(abs(val-fp(i,j,k)).gt.1.e-8) print*, 'Large difference : ', val-fp(i,j,k),ii,jj,kk
          enddo
        enddo
      enddo
    case('f')
      diffmax = 0.
      do k=1,n(3)-q(3)
        kp = k + 1
        km = k - 1
        do j=1,n(2)-q(2)
          jp = j + 1
          jm = j - 1
          do i=1,n(1)-q(1)
            ip = i + 1
            im = i - 1
            val =  fpp(i,j,k)+(1./alpha)*( &
                  (fpp(ip,j,k)-2.*fpp(i,j,k)+fpp(im,j,k))*(dli(1)**2) + &
                  (fpp(i,jp,k)-2.*fpp(i,j,k)+fpp(i,jm,k))*(dli(2)**2) + &
                 ((fpp(i,j,kp)-fpp(i,j,k ))*dzfi(kp) - &
                  (fpp(i,j,k )-fpp(i,j,km))*dzfi(k ))*dzci(k) )
            val = val*alpha
            diffmax = max(diffmax,abs(val-fp(i,j,k)))
            !ii = ijk_min(1)+i
            !jj = ijk_min(2)+j
            !if(abs(val-fp(i,j,k)).gt.1.e-8) print*, 'Large difference : ', val,fp(i,j,k),ii,jj,k
          enddo
        enddo
      enddo
    end select
    call mpi_allreduce(MPI_IN_PLACE,diffmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    return
  end subroutine chk_helmholtz
end module mod_debug
