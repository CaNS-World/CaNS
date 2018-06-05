module mod_debug
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_param     , only: dims
  implicit none
  private
  public chkmean,chkhelmholtz
  contains
  subroutine chkmean(n,dzlzi,p,mean)
    implicit none
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension(0:) :: dzlzi
    real(8), intent(in), dimension(0:,0:,0:) :: p
    real(8), intent(out) :: mean
    integer :: i,j,k
    mean = 0.d0
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
    call mpi_allreduce(MPI_IN_PLACE,mean,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean = mean/(1.d0*n(1)*dims(1)*n(2)*dims(2))
    return
  end subroutine chkmean
  !
  subroutine chkhelmholtz(n,dxi,dyi,dzci,dzfi,alpha,fp,fpp,c_or_f_z)
    implicit none
    integer, intent(in), dimension(3) :: n
    real(8), intent(in) :: dxi,dyi,alpha
    real(8), intent(in), dimension(0:) :: dzfi,dzci
    real(8), intent(in), dimension(0:,0:,0:) :: fp,fpp
    character(len=1), intent(in) :: c_or_f_z
    real(8) :: val
    integer :: i,j,k,im,ip,jm,jp,km,kp
    select case(c_or_f_z)
    case('c')
      do k=1,n(3)-1
        kp = k + 1
        km = k - 1
        do j=1,n(2)-1
          jp = j + 1
          jm = j - 1
          do i=1,n(1)-1
            ip = i + 1
            im = i - 1
            val =  fpp(i,j,k)+(1./alpha)*( &
                  (fpp(ip,j,k)-2.d0*fpp(i,j,k)+fpp(im,j,k))*(dxi**2) + &
                  (fpp(i,jp,k)-2.d0*fpp(i,j,k)+fpp(i,jm,k))*(dyi**2) + &
                 ((fpp(i,j,kp)-fpp(i,j,k))*dzci(k) - &
                  (fpp(i,j,k )-fpp(i,j,km))*dzci(km))*dzfi(k) )
            val = val*alpha
            if(abs(val-fp(i,j,k)).gt.1.e-9) print*, 'Large difference : ', val,fp(i,j,k),i,j,k
          enddo
        enddo
      enddo
    case('f')
      do k=1,n(3)-1
        kp = k + 1
        km = k - 1
        do j=1,n(2)-1
          jp = j + 1
          jm = j - 1
          do i=1,n(1)-1
            ip = i + 1
            im = i - 1
            val =  fpp(i,j,k)+(1./alpha)*( &
                  (fpp(ip,j,k)-2.d0*fpp(i,j,k)+fpp(im,j,k))*(dxi**2) + &
                  (fpp(i,jp,k)-2.d0*fpp(i,j,k)+fpp(i,jm,k))*(dyi**2) + &
                 ((fpp(i,j,kp)-fpp(i,j,k))*dzfi(kp) - &
                  (fpp(i,j,k )-fpp(i,j,km))*dzfi(k))*dzci(k) )
            val = val*alpha
            if(abs(val-fp(i,j,k)).gt.1.e-9) print*, 'Large difference : ', val,fp(i,j,k),i,j,k
          enddo
        enddo
      enddo
    end select
    return
  end subroutine chkhelmholtz
end module mod_debug
