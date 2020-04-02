module mod_output
  use mpi
  use decomp_2d_io
  use mod_param     , only: dims,dx,dy,dz
  use mod_common_mpi, only:ierr,myid,ipencil,ijk_start
  use mod_types
  implicit none
  private
  public out0d,out1d,out1d_2,out2d,out3d
  contains
  subroutine out0d(fname,n,var)
    !
    ! appends the first n entries of an array
    ! var to a file
    ! fname -> name of the file
    ! n     -> number of entries
    ! var   -> input array of real values
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in) :: n
    real(rp), intent(in), dimension(:) :: var
    integer :: iunit
    character(len=30) :: cfmt
    integer :: i
    !
    write(cfmt,'(A,I3,A)') '(',n,'E15.7)'
    iunit = 10
    if (myid .eq. 0) then
      open(iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) (var(i),i=1,n) 
      close(iunit)
    endif
    return
  end subroutine out0d
  !
  subroutine out1d(fname,n,idir,z,dzlzi,p)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! fname -> name of the file
    ! n     -> size of the input array
    ! idir  -> direction of the profile
    ! z     -> z coordinate (grid is non-uniform in z)
    ! dzlzi -> dz/lz weight of a grid cell for averaging over z
    ! p     -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: n
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(0:) :: z,dzlzi
    real(rp), intent(in), dimension(0:,0:,0:) :: p
    real(rp), allocatable, dimension(:) :: p1d
    integer :: i,j,k,ii,jj,kk
    integer :: iunit
    integer, dimension(3) :: ng
    !
    ng(:) = n(:)*dims(:)
    iunit = 10
    select case(idir)
    case(3)
      allocate(p1d(ng(3)))
      p1d(:) = 0.
      do k=1,n(3)
        kk = ijk_start(3) + k
        p1d(kk) = 0.
        do j=1,n(2)
          do i=1,n(1)
            p1d(kk) = p1d(kk) + p(i,j,k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d(:) = p1d(:)/(1.*ng(1)*ng(2))
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do kk=1,ng(3)
          write(iunit,'(2E15.7)') z(kk),p1d(kk)
        enddo
        close(iunit)
      endif
    case(2)
      allocate(p1d(ng(2)))
      p1d(:) = 0.
      do j=1,n(2)
        jj = ijk_start(2)+j
        p1d(jj) = 0.
        do k=1,n(3)
          do i=1,n(1)
            p1d(jj) = p1d(jj) + p(i,j,k)*dzlzi(k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng(2),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d(:) = p1d(:)/(1.*ng(1))
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do jj=1,ng(2)
          write(iunit,'(2E15.7)') (1.*jj-.5)/(1.*ng(2)),p1d(jj)
        enddo
        close(iunit)
      endif
    case(1)
      allocate(p1d(ng(1)))
      p1d(:) = 0.
      do i=1,n(1)
        ii = ijk_start(1)+i
        p1d(ii) = 0.
        do k=1,n(3)
          do j=1,n(2)
            p1d(ii) = p1d(ii) + p(i,j,k)*dzlzi(k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng(1),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d(:) = p1d(:)/(1.*ng(2))
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do i=1,ng(1)
          write(iunit,'(2E15.7)') (1.*i-.5)/(1.*n(1)),p1d(j)
        enddo
        close(iunit)
      endif
    end select
    deallocate(p1d)
  end subroutine out1d
  !
  !
  subroutine out2d(fname,inorm,islice,p)
    !
    ! saves a planar slice of a scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! inorm  -> plane is perpendicular to direction
    !           inorm (1,2,3)
    ! islice -> plane is of constant index islice 
    !           in direction inorm
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in) :: inorm,islice
    real(rp),intent(in), dimension(:,:,:) :: p
    !
    select case(inorm)
    case(1) !normal to x --> yz plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,fname)
    case(2) !normal to y --> zx plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,fname)
    case(3) !normal to z --> xy plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,fname)
    end select
    return
  end subroutine out2d
  !
  subroutine out3d(fname,nskip,p)
    !
    ! saves a 3D scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! nskip  -> array with the step size for which the
    !           field is written; i.e.: (/1,1,1/) 
    !           writes the full field 
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: nskip
    real(rp),intent(in), dimension(:,:,:) :: p
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    !
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fname, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_every(ipencil,p,nskip(1),nskip(2),nskip(3),fname,.true.)
    call MPI_FILE_CLOSE(fh,ierr)
    return
  end subroutine out3d
  !
  subroutine out1d_2(fname,n,idir,z,u,v,w) ! e.g. for a channel with streamwise dir in x
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: n
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(0:) :: z
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), allocatable, dimension(:) :: um,vm,wm,u2,v2,w2,uw
    integer :: i,j,k,kk
    integer :: iunit
    integer, dimension(3) :: ng
    integer :: q
    !
    ng(:) = n(:)*dims(:)
    iunit = 10
    select case(idir)
    case(3)
      q = ng(3)
      allocate(um(0:q+1),vm(0:q+1),wm(0:q+1),u2(0:q+1),v2(0:q+1),w2(0:q+1),uw(0:q+1))
      um(:) = 0.
      vm(:) = 0.
      wm(:) = 0.
      u2(:) = 0.
      v2(:) = 0.
      w2(:) = 0.
      uw(:) = 0.
      do k=1,n(3)
        kk = k+ijk_start(3)
        um(kk) = 0.
        vm(kk) = 0.
        wm(kk) = 0.
        u2(kk) = 0.
        v2(kk) = 0.
        w2(kk) = 0.
        uw(kk) = 0.
        do j=1,n(2)
          do i=1,n(1)
            um(kk) = um(kk) + u(i,j,k)
            vm(kk) = vm(kk) + v(i,j,k)
            wm(kk) = wm(kk) + 0.50*(w(i,j,k-1) + w(i,j,k))
            u2(kk) = u2(kk) + u(i,j,k)**2
            v2(kk) = v2(kk) + v(i,j,k)**2
            w2(kk) = w2(kk) + 0.50*(w(i,j,k)**2+w(i,j,k-1)**2)
            uw(kk) = uw(kk) + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                   (w(i,j,k-1) + w(i,j,k))
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,um(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,wm(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,u2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,v2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,w2(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,uw(1),ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:) = um(:)/(1.*ng(1)*ng(2))
      vm(:) = vm(:)/(1.*ng(1)*ng(2))
      wm(:) = wm(:)/(1.*ng(1)*ng(2))
      u2(:) = sqrt(u2(:)/(1.*ng(1)*ng(2)) - um(:)**2)
      v2(:) = sqrt(v2(:)/(1.*ng(1)*ng(2)) - vm(:)**2)
      w2(:) = sqrt(w2(:)/(1.*ng(1)*ng(2)) - wm(:)**2)
      uw(:) = uw(:)/(1.*ng(1)*ng(2)) - um(:)*wm(:)
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do kk=1,ng(3)
          write(iunit,'(8E15.7)') z(kk),um(kk),vm(kk),wm(kk), &
                                        u2(kk),v2(kk),w2(kk), &
                                        uw(kk)
        enddo
        close(iunit)
      endif
      deallocate(um,vm,wm,u2,v2,w2,uw)
    case(2)
    case(1)
    end select
  end subroutine out1d_2
  !
  subroutine out2d_2(fname,n,idir,z,u,v,w) ! e.g. for a duct
    !
    ! UPDATE CASE WITH STREAMWISE X
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: n
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(0:) :: z
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), allocatable, dimension(:,:) :: um,vm,wm,u2,v2,w2,uv,vw
    integer :: i,j,k,ii,jj,kk
    integer :: iunit
    integer, dimension(3) :: ng
    integer :: p,q
    real(rp) :: x
    !
    ng(:) = n(:)*dims(:)
    iunit = 10
    select case(idir)
    case(3)
      p = ng(1)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),vw(p,q))
      !
      um(:,:) = 0.
      vm(:,:) = 0.
      wm(:,:) = 0.
      u2(:,:) = 0.
      v2(:,:) = 0.
      w2(:,:) = 0.
      uv(:,:) = 0.
      vw(:,:) = 0.
      do k=1,n(3)
        kk = ijk_start(3) + k
        do i=1,n(1)
          ii = ijk_start(1) + i
          um(ii,kk) = 0.
          vm(ii,kk) = 0.
          wm(ii,kk) = 0.
          u2(ii,kk) = 0.
          v2(ii,kk) = 0.
          w2(ii,kk) = 0.
          vw(ii,kk) = 0.
          uv(ii,kk) = 0.
          do j=1,n(2)
            jj = ijk_start(2) + j
            um(ii,kk) = um(ii,kk) + 0.5*(u(i-1,j,k)+u(i,j,k))
            vm(ii,kk) = vm(ii,kk) + v(i,j,k)
            wm(ii,kk) = wm(ii,kk) + 0.5*(w(i,j,k-1)+w(i,j,k))
            u2(ii,kk) = u2(ii,kk) + 0.5*(u(i-1,j,k)**2+u(i,j,k)**2)
            v2(ii,kk) = v2(ii,kk) + v(i,j,k)**2
            w2(ii,kk) = w2(ii,kk) + 0.5*(w(i,j,k-1)**2+w(i,j,k)**2)
            vw(ii,kk) = vw(ii,kk) + 0.25*(v(i,j-1,k) + v(i,j,k))* &
                                         (w(i,j,k-1) + w(i,j,k))
            uv(ii,kk) = uv(ii,kk) + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                         (v(i,j-1,k) + v(i,j,k))
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,um(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,wm(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,u2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,v2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,w2(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vw(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,uv(1,1),ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) =      um(:,:)/(1.*ng(2))
      vm(:,:) =      vm(:,:)/(1.*ng(2))
      wm(:,:) =      wm(:,:)/(1.*ng(2))
      u2(:,:) = sqrt(u2(:,:)/(1.*ng(2)) - um(:,:)**2)
      v2(:,:) = sqrt(v2(:,:)/(1.*ng(2)) - vm(:,:)**2)
      w2(:,:) = sqrt(w2(:,:)/(1.*ng(2)) - wm(:,:)**2)
      vw(:,:) =      vw(:,:)/(1.*ng(2)) - vm(:,:)*wm(:,:)
      uv(:,:) =      uv(:,:)/(1.*ng(2)) - um(:,:)*vm(:,:)
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do kk=1,ng(3)
          do ii=1,ng(1)
            x = (ii-.5)*dx
            write(iunit,'(10E15.7)') x,z(kk),um(ii,kk),vm(ii,kk),wm(ii,kk), &
                                             u2(ii,kk),v2(ii,kk),w2(ii,kk), &
                                             vw(ii,kk),uv(ii,kk)
          enddo
        enddo
        close(iunit)
      endif
      deallocate(um,vm,wm,u2,v2,w2,vw,uv)
    case(2)
    case(1)
    end select
  end subroutine out2d_2
end module mod_output
