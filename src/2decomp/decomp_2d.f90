!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This is the main 2D pencil decomposition module

module decomp_2d

  use MPI

  implicit none

#ifdef GLOBAL_ARRAYS
#include "mafdecls.fh"
#include "global.fh"
#endif

  private        ! Make everything private unless declared public

#ifndef SINGLE
  integer, parameter, public :: mytype = KIND(0.0D0)
  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
  integer, parameter, public :: complex_type = MPI_DOUBLE_COMPLEX
#ifdef GLOBAL_ARRAYS
  integer, parameter, public :: ga_real_type = MT_F_DBL
  integer, parameter, public :: ga_complex_type = MT_F_DCPL
#endif
#else
  integer, parameter, public :: mytype = KIND(0.0)
  integer, parameter, public :: real_type = MPI_REAL
  integer, parameter, public :: complex_type = MPI_COMPLEX
#ifdef GLOBAL_ARRAYS
  integer, parameter, public :: ga_real_type = MT_F_REAL
  integer, parameter, public :: ga_complex_type = MT_F_SCPL
#endif
#endif

  integer, save, public :: mytype_bytes

  ! some key global variables
  integer, save, public :: nx_global, ny_global, nz_global  ! global size

  integer, save, public :: nrank  ! local MPI rank 
  integer, save, public :: nproc  ! total number of processors

  ! parameters for 2D Cartesian topology 
  integer, save, dimension(2) :: dims, coord
  logical, save, dimension(2) :: periodic
  integer, save, public :: DECOMP_2D_COMM_CART_X, &
       DECOMP_2D_COMM_CART_Y, DECOMP_2D_COMM_CART_Z 
  integer, save :: DECOMP_2D_COMM_ROW, DECOMP_2D_COMM_COL

  ! define neighboring blocks (to be used in halo-cell support)
  !  first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
  ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom 
  integer, save, dimension(3,6) :: neighbour 

  ! flags for periodic condition in three dimensions
  logical, save :: periodic_x, periodic_y, periodic_z

#ifdef SHM
  ! derived type to store shared-memory info
  TYPE, public :: SMP_INFO
     integer MPI_COMM          ! SMP associated with this communicator
     integer NODE_ME           ! rank in this communicator
     integer NCPU              ! size of this communicator
     integer SMP_COMM          ! communicator for SMP-node masters
     integer CORE_COMM         ! communicator for cores on SMP-node
     integer SMP_ME            ! SMP-node id starting from 1 ... NSMP
     integer NSMP              ! number of SMP-nodes in this communicator
     integer CORE_ME           ! core id starting from 1 ... NCORE
     integer NCORE             ! number of cores on this SMP-node
     integer MAXCORE           ! maximum no. cores on any SMP-node
     integer N_SND             ! size of SMP shared memory buffer
     integer N_RCV             ! size of SMP shared memory buffer
     integer(8) SND_P          ! SNDBUF address (cray pointer), for real 
     integer(8) RCV_P          ! RCVBUF address (cray pointer), for real
     integer(8) SND_P_c        ! for complex
     integer(8) RCV_P_c        ! for complex
  END TYPE SMP_INFO
#endif

  ! derived type to store decomposition info for a given global data size
  TYPE, public :: DECOMP_INFO
     ! staring/ending index and size of data held by current processor
     integer, dimension(3) :: xst, xen, xsz  ! x-pencil
     integer, dimension(3) :: yst, yen, ysz  ! y-pencil
     integer, dimension(3) :: zst, zen, zsz  ! z-pencil

     ! in addition to local information, processors also need to know 
     ! some global information for global communications to work 

     ! how each dimension is distributed along pencils
     integer, allocatable, dimension(:) :: &
          x1dist, y1dist, y2dist, z2dist, x2dist, z1dist, &
          x1st, y1st, y2st, z2st, x2st, z1st, &
          x1en, y1en, y2en, z2en, x2en, z1en

     ! send/receive buffer counts and displacements for MPI_ALLTOALLV
     integer, allocatable, dimension(:) :: &
          x1cnts, y1cnts, y2cnts, z2cnts, x2cnts, z1cnts
     integer, allocatable, dimension(:) :: &
          x1disp, y1disp, y2disp, z2disp, x2disp, z1disp

     ! buffer counts for MPI_ALLTOALL: either for evenly distributed data
     ! or for padded-alltoall
     integer :: x1count, y1count, y2count, z2count, x2count, z1count
!
     ! buffer counts, displacements and types for MPI_Alltoallw to transform
     ! directly between x- and z-pencils
     ! This is only for the complex datatype
     integer,dimension(:),allocatable::zcnts_xz,ztypes_xz
     integer,dimension(:),allocatable::xcnts_xz,xtypes_xz
#ifdef MPI3
! use MPI_ADDRESS_KIND for MPI_Neighbor_alltoallw call
     integer(kind=MPI_ADDRESS_KIND),dimension(:),allocatable::zdispls_xz,xdispls_xz
     integer :: xtozNeighborComm,ztoxNeighborComm
     integer,dimension(:),allocatable::xranks,zranks
#else
! use default integer for MPI_Alltoallw call
     integer,dimension(:),allocatable::zdispls_xz,xdispls_xz
#endif

     ! evenly distributed data
     logical :: even


#ifdef SHM
     ! For shared-memory implementation

     ! one instance of this derived type for each communicator
     ! shared moemory info, such as which MPI rank belongs to which node
     TYPE(SMP_INFO) :: ROW_INFO, COL_INFO

     ! shared send/recv buffers for ALLTOALLV
     integer, allocatable, dimension(:) :: x1cnts_s, y1cnts_s, &
          y2cnts_s, z2cnts_s
     integer, allocatable, dimension(:) :: x1disp_s, y1disp_s, &
          y2disp_s, z2disp_s
     ! A copy of original buffer displacement (will be overwriten)
     integer, allocatable, dimension(:) :: x1disp_o, y1disp_o, &
          y2disp_o, z2disp_o
#endif
  END TYPE DECOMP_INFO

  ! main (default) decomposition information for global size nx*ny*nz
  TYPE(DECOMP_INFO), save :: decomp_main
  TYPE(DECOMP_INFO), save :: decomp_ph,decomp_sp

  ! staring/ending index and size of data held by current processor
  ! duplicate 'decomp_main', needed by apps to define data structure 
  integer, save, dimension(3), public :: xstart, xend, xsize  ! x-pencil
  integer, save, dimension(3), public :: ystart, yend, ysize  ! y-pencil
  integer, save, dimension(3), public :: zstart, zend, zsize  ! z-pencil

  ! These are the buffers used by MPI_ALLTOALL(V) calls
  integer, save :: decomp_buf_size = 0
  real(mytype),    allocatable, dimension(:) :: work1_r, work2_r
  complex(mytype), allocatable, dimension(:) :: work1_c, work2_c

  ! public user routines
  public :: decomp_2d_init, decomp_2d_finalize, &
       transpose_x_to_y, transpose_y_to_z, &
       transpose_z_to_y, transpose_y_to_x, &
       transpose_z_to_x, transpose_x_to_z, &
#ifdef OCC
       transpose_x_to_y_start, transpose_y_to_z_start, &
       transpose_z_to_y_start, transpose_y_to_x_start, &
       transpose_x_to_y_wait, transpose_y_to_z_wait, &
       transpose_z_to_y_wait, transpose_y_to_x_wait, &
       transpose_test, &
#endif
       decomp_info_init, decomp_info_finalize, partition, &
#ifdef GLOBAL_ARRAYS
       get_global_array, &
#endif
       alloc_x, alloc_y, alloc_z, &
       update_halo, decomp_2d_abort, &
       get_decomp_info


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are routines to perform global data transpositions
  ! 
  !   Four combinations are available, enough to cover all situations
  !    - transpose_x_to_y (X-pencil --> Y-pencil)
  !    - transpose_y_to_z (Y-pencil --> Z-pencil)
  !    - transpose_z_to_y (Z-pencil --> Y-pencil)
  !    - transpose_y_to_x (Y-pencil --> X-pencil)
  !
  !   Generic interface provided here to support multiple data types
  !    - real and complex types supported through generic interface
  !    - single/double precision supported through pre-processing
  !       * see 'mytype' variable at the beginning
  !    - an optional argument can be supplied to transpose data whose 
  !      global size is not the default nx*ny*nz 
  !       * as the case in fft r2c/c2r interface 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface transpose_x_to_y
     module procedure transpose_x_to_y_real
     module procedure transpose_x_to_y_complex
  end interface transpose_x_to_y
  
  interface transpose_y_to_z
     module procedure transpose_y_to_z_real
     module procedure transpose_y_to_z_complex
  end interface transpose_y_to_z
  
  interface transpose_z_to_y
     module procedure transpose_z_to_y_real
     module procedure transpose_z_to_y_complex
  end interface transpose_z_to_y

  interface transpose_y_to_x
     module procedure transpose_y_to_x_real
     module procedure transpose_y_to_x_complex
  end interface transpose_y_to_x

  interface transpose_z_to_x
! not available
     module procedure transpose_z_to_x_real
!     module procedure transpose_z_to_x_complex
  end interface transpose_z_to_x

  interface transpose_x_to_z
! not available
     module procedure transpose_x_to_z_real
!     module procedure transpose_x_to_z_complex
  end interface transpose_x_to_z

#ifdef OCC
  interface transpose_x_to_y_start
     module procedure transpose_x_to_y_real_start
     module procedure transpose_x_to_y_complex_start
  end interface transpose_x_to_y_start

  interface transpose_y_to_z_start
     module procedure transpose_y_to_z_real_start
     module procedure transpose_y_to_z_complex_start
  end interface transpose_y_to_z_start

  interface transpose_z_to_y_start
     module procedure transpose_z_to_y_real_start
     module procedure transpose_z_to_y_complex_start
  end interface transpose_z_to_y_start
     
  interface transpose_y_to_x_start
     module procedure transpose_y_to_x_real_start
     module procedure transpose_y_to_x_complex_start
  end interface transpose_y_to_x_start

  interface transpose_x_to_y_wait
     module procedure transpose_x_to_y_real_wait
     module procedure transpose_x_to_y_complex_wait
  end interface transpose_x_to_y_wait

  interface transpose_y_to_z_wait
     module procedure transpose_y_to_z_real_wait
     module procedure transpose_y_to_z_complex_wait
  end interface transpose_y_to_z_wait

  interface transpose_z_to_y_wait
     module procedure transpose_z_to_y_real_wait
     module procedure transpose_z_to_y_complex_wait
  end interface transpose_z_to_y_wait
     
  interface transpose_y_to_x_wait
     module procedure transpose_y_to_x_real_wait
     module procedure transpose_y_to_x_complex_wait
  end interface transpose_y_to_x_wait
#endif

  interface update_halo
     module procedure update_halo_real
     module procedure update_halo_complex
  end interface update_halo

  interface alloc_x
     module procedure alloc_x_real
     module procedure alloc_x_complex
  end interface alloc_x

  interface alloc_y
     module procedure alloc_y_real
     module procedure alloc_y_complex
  end interface alloc_y

  interface alloc_z
     module procedure alloc_z_real
     module procedure alloc_z_complex
  end interface alloc_z

contains

#ifdef SHM_DEBUG
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! For debugging, print the shared-memory structure
  subroutine print_smp_info(s)
    TYPE(SMP_INFO) :: s
    write(10,*) 'size of current communicator:', s%NCPU
    write(10,*) 'rank in current communicator:', s%NODE_ME
    write(10,*) 'number of SMP-nodes in this communicator:', s%NSMP
    write(10,*) 'SMP-node id (1 ~ NSMP):', s%SMP_ME
    write(10,*) 'NCORE - number of cores on this SMP-node', s%NCORE
    write(10,*) 'core id (1 ~ NCORE):', s%CORE_ME
    write(10,*) 'maximum no. cores on any SMP-node:', s%MAXCORE
    write(10,*) 'size of SMP shared memory SND buffer:', s%N_SND
    write(10,*) 'size of SMP shared memory RCV buffer:', s%N_RCV
  end subroutine print_smp_info
#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to initialise this library
  !   INPUT:
  !     nx, ny, nz   - global data dimension
  !     p_row, p_col - 2D processor grid
  !   OUTPUT:
  !     all internal data structures initialised properly
  !     library ready to use
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_init(nx,ny,nz,p_row,p_col,periodic_bc)

    implicit none

    integer, intent(IN) :: nx,ny,nz,p_row,p_col
    logical, dimension(3), intent(IN), optional :: periodic_bc
    
    integer :: errorcode, ierror, row, col
    
#ifdef SHM_DEBUG
    character(len=80) fname
#endif

    nx_global = nx
    ny_global = ny
    nz_global = nz

    if (present(periodic_bc)) then
       periodic_x = periodic_bc(1)
       periodic_y = periodic_bc(2)
       periodic_z = periodic_bc(3)
    else
       periodic_x = .false.
       periodic_y = .false.
       periodic_z = .false.
    end if

    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

    if (p_row==0 .and. p_col==0) then
       ! determine the best 2D processor grid
       call best_2d_grid(nproc, row, col)
    else
       if (nproc /= p_row*p_col) then
          errorcode = 1
          call decomp_2d_abort(errorcode, &
               'Invalid 2D processor grid - nproc /= p_row*p_col')
       else
          row = p_row
          col = p_col
       end if
    end if
    
    ! Create 2D Catersian topology
    ! Note that in order to support periodic B.C. in the halo-cell code,
    ! need to create multiple topology objects: DECOMP_2D_COMM_CART_?,
    ! corresponding to three pencil orientations. They contain almost
    ! identical topological information but allow different combinations
    ! of periodic conditions.
    dims(1) = row
    dims(2) = col
    periodic(1) = periodic_y
    periodic(2) = periodic_z
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., &  ! do not reorder rank
         DECOMP_2D_COMM_CART_X, ierror)
    periodic(1) = periodic_x
    periodic(2) = periodic_z
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., DECOMP_2D_COMM_CART_Y, ierror)
    periodic(1) = periodic_x
    periodic(2) = periodic_y
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., DECOMP_2D_COMM_CART_Z, ierror)

    call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
    
    ! derive communicators defining sub-groups for ALLTOALL(V)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
         DECOMP_2D_COMM_COL,ierror)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
         DECOMP_2D_COMM_ROW,ierror)

    ! gather information for halo-cell support code
    call init_neighbour
    
    ! actually generate all 2D decomposition information
    call decomp_info_init(nx,ny,nz,decomp_main)
    
    ! make a copy of the decomposition information associated with the
    ! default global size in these global variables so applications can
    ! use them to create data structures 
    xstart = decomp_main%xst
    ystart = decomp_main%yst
    zstart = decomp_main%zst
    xend   = decomp_main%xen
    yend   = decomp_main%yen
    zend   = decomp_main%zen
    xsize  = decomp_main%xsz
    ysize  = decomp_main%ysz
    zsize  = decomp_main%zsz

    decomp_ph=decomp_main

#ifdef SHM_DEBUG
    ! print out shared-memory information
    write(fname,99) nrank
99  format('log',I2.2)
    open(10,file=fname)
    write(10,*)'I am mpi rank ', nrank, 'Total ranks ', nproc
    write(10,*)' '
    write(10,*)'Global data size:'
    write(10,*)'nx*ny*nz', nx,ny,nz
    write(10,*)' '
    write(10,*)'2D processor grid:'
    write(10,*)'p_row*p_col:', dims(1), dims(2)
    write(10,*)' '
    write(10,*)'Portion of global data held locally:'
    write(10,*)'xsize:',xsize
    write(10,*)'ysize:',ysize
    write(10,*)'zsize:',zsize
    write(10,*)' '
    write(10,*)'How pensils are to be divided and sent in alltoallv:'
    write(10,*)'x1dist:',decomp_main%x1dist
    write(10,*)'y1dist:',decomp_main%y1dist
    write(10,*)'y2dist:',decomp_main%y2dist
    write(10,*)'z2dist:',decomp_main%z2dist
    write(10,*)' '
    write(10,*)'######Shared buffer set up after this point######'
    write(10,*)' '
    write(10,*) 'col communicator detais:'
    call print_smp_info(decomp_main%COL_INFO)
    write(10,*)' '
    write(10,*) 'row communicator detais:'
    call print_smp_info(decomp_main%ROW_INFO)
    write(10,*)' '
    write(10,*)'Buffer count and dispalcement of per-core buffers'
    write(10,*)'x1cnts:',decomp_main%x1cnts
    write(10,*)'y1cnts:',decomp_main%y1cnts
    write(10,*)'y2cnts:',decomp_main%y2cnts
    write(10,*)'z2cnts:',decomp_main%z2cnts
    write(10,*)'x1disp:',decomp_main%x1disp
    write(10,*)'y1disp:',decomp_main%y1disp
    write(10,*)'y2disp:',decomp_main%y2disp
    write(10,*)'z2disp:',decomp_main%z2disp
    write(10,*)' '
    write(10,*)'Buffer count and dispalcement of shared buffers'
    write(10,*)'x1cnts:',decomp_main%x1cnts_s
    write(10,*)'y1cnts:',decomp_main%y1cnts_s
    write(10,*)'y2cnts:',decomp_main%y2cnts_s
    write(10,*)'z2cnts:',decomp_main%z2cnts_s
    write(10,*)'x1disp:',decomp_main%x1disp_s
    write(10,*)'y1disp:',decomp_main%y1disp_s
    write(10,*)'y2disp:',decomp_main%y2disp_s
    write(10,*)'z2disp:',decomp_main%z2disp_s
    write(10,*)' '
    close(10)
#endif

    ! determine the number of bytes per float number
    ! do not use 'mytype' which is compiler dependent
    ! also possible to use inquire(iolength=...) 
    call MPI_TYPE_SIZE(real_type,mytype_bytes,ierror)

#ifdef EVEN
    if (nrank==0) write(*,*) 'Padded ALLTOALL optimisation on'
#endif 

    return
  end subroutine decomp_2d_init
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to clean things up
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_finalize

    implicit none
    
    call decomp_info_finalize(decomp_main)

    decomp_buf_size = 0
    deallocate(work1_r, work2_r, work1_c, work2_c)
    
    return
  end subroutine decomp_2d_finalize


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the default decomposition object
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_decomp_info(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(OUT) :: decomp

    decomp = decomp_main

    return
  end subroutine get_decomp_info
    

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Advanced Interface allowing applications to define globle domain of
  ! any size, distribute it, and then transpose data among pencils.
  !  - generate 2D decomposition details as defined in DECOMP_INFO
  !  - the default global data size is nx*ny*nz
  !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
  !  - multiple global sizes can co-exist in one application, each
  !    using its own DECOMP_INFO object
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_init(nx,ny,nz,decomp)

    implicit none
    
    integer, intent(IN) :: nx,ny,nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    integer :: buf_size, status, errorcode

    ! verify the global size can actually be distributed as pencils
    if (nx<dims(1) .or. ny<dims(1) .or. ny<dims(2) .or. nz<dims(2)) then
       errorcode = 6
       call decomp_2d_abort(errorcode, &
            'Invalid 2D processor grid. ' // &
            'Make sure that min(nx,ny) >= p_row and ' // &
            'min(ny,nz) >= p_col')
    end if
    
    if (mod(nx,dims(1))==0 .and. mod(ny,dims(1))==0 .and. &
         mod(ny,dims(2))==0 .and. mod(nz,dims(2))==0) then
       decomp%even = .true.
    else
       decomp%even = .false.
    end if

    ! distribute mesh points
    !allocate(decomp%x1dist(0:dims(1)-1),decomp%y1dist(0:dims(1)-1), &
    !     decomp%y2dist(0:dims(2)-1),decomp%z2dist(0:dims(2)-1))
    allocate(decomp%x1dist(0:dims(1)-1),decomp%y1dist(0:dims(1)-1), &
         decomp%y2dist(0:dims(2)-1),decomp%z2dist(0:dims(2)-1), &
         decomp%x2dist(0:dims(2)-1),decomp%z1dist(0:dims(1)-1))
    allocate(decomp%x1st(0:dims(1)-1),decomp%x1en(0:dims(1)-1), &
         decomp%y1st(0:dims(1)-1),decomp%y1en(0:dims(1)-1), &
         decomp%z1st(0:dims(1)-1),decomp%z1en(0:dims(1)-1))
    allocate(decomp%x2st(0:dims(2)-1),decomp%x2en(0:dims(2)-1), &
         decomp%y2st(0:dims(2)-1),decomp%y2en(0:dims(2)-1), &
         decomp%z2st(0:dims(2)-1),decomp%z2en(0:dims(2)-1))
    call get_dist(nx,ny,nz,decomp)
    
    ! generate partition information - starting/ending index etc.
    call partition(nx, ny, nz, (/ 1,2,3 /), &
         decomp%xst, decomp%xen, decomp%xsz)
    call partition(nx, ny, nz, (/ 2,1,3 /), &
         decomp%yst, decomp%yen, decomp%ysz)
    call partition(nx, ny, nz, (/ 2,3,1 /), &
         decomp%zst, decomp%zen, decomp%zsz)
    
    ! prepare send/receive buffer displacement and count for ALLTOALL(V)
    !allocate(decomp%x1cnts(0:dims(1)-1),decomp%y1cnts(0:dims(1)-1), &
    !     decomp%y2cnts(0:dims(2)-1),decomp%z2cnts(0:dims(2)-1))
    !allocate(decomp%x1disp(0:dims(1)-1),decomp%y1disp(0:dims(1)-1), &
    !     decomp%y2disp(0:dims(2)-1),decomp%z2disp(0:dims(2)-1))
    allocate(decomp%x1cnts(0:dims(1)-1),decomp%y1cnts(0:dims(1)-1), &
         decomp%y2cnts(0:dims(2)-1),decomp%z2cnts(0:dims(2)-1), &
         decomp%z1cnts(0:dims(1)-1),decomp%x2cnts(0:dims(2)-1))
    allocate(decomp%x1disp(0:dims(1)-1),decomp%y1disp(0:dims(1)-1), &
         decomp%y2disp(0:dims(2)-1),decomp%z2disp(0:dims(2)-1), &
         decomp%x2disp(0:dims(2)-1),decomp%z1disp(0:dims(1)-1))
     ! allocate arrays for MPI_ALLTOALLW calls
    allocate(decomp%xcnts_xz(nproc),decomp%zcnts_xz(nproc))
    allocate(decomp%xtypes_xz(nproc),decomp%ztypes_xz(nproc))
    allocate(decomp%xdispls_xz(nproc))
    allocate(decomp%zdispls_xz(nproc))
    call prepare_buffer(decomp)

#ifdef SHM
    ! prepare shared-memory information if required
    call decomp_info_init_shm(decomp)
#endif

    ! allocate memory for the MPI_ALLTOALL(V) buffers
    ! define the buffers globally for performance reason
    
    buf_size = max(decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3), &
         max(decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3), &
         decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)) )
#ifdef EVEN
    ! padded alltoall optimisation may need larger buffer space
    buf_size = max(buf_size, &
         max(decomp%x1count*dims(1),decomp%y2count*dims(2)) ) 
#endif

    ! check if additional memory is required
    ! *** TODO: consider how to share the real/complex buffers 
    if (buf_size > decomp_buf_size) then
       decomp_buf_size = buf_size
       if (allocated(work1_r)) deallocate(work1_r)
       if (allocated(work2_r)) deallocate(work2_r)
       if (allocated(work1_c)) deallocate(work1_c)
       if (allocated(work2_c)) deallocate(work2_c)
       allocate(work1_r(buf_size), STAT=status)
       allocate(work2_r(buf_size), STAT=status)
       allocate(work1_c(buf_size), STAT=status)
       allocate(work2_c(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
    end if

    return
  end subroutine decomp_info_init


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Release memory associated with a DECOMP_INFO object
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_finalize(decomp)

    implicit none

    integer :: i
    integer :: ierror
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    deallocate(decomp%x1dist,decomp%y1dist,decomp%y2dist,decomp%z2dist)
    deallocate(decomp%z1dist,decomp%x2dist)
    deallocate(decomp%x1st,decomp%y1st,decomp%y2st,decomp%z2st)
    deallocate(decomp%z1st,decomp%x2st)
    deallocate(decomp%x1en,decomp%y1en,decomp%y2en,decomp%z2en)
    deallocate(decomp%z1en,decomp%x2en)
    deallocate(decomp%x1cnts,decomp%y1cnts,decomp%y2cnts,decomp%z2cnts)
    deallocate(decomp%z1cnts,decomp%x2cnts)
    deallocate(decomp%x1disp,decomp%y1disp,decomp%y2disp,decomp%z2disp)
    deallocate(decomp%z1disp,decomp%x2disp)
    do i=1,nproc
      if (decomp%ztypes_xz(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%ztypes_xz(i),ierror)
      endif
      if (decomp%xtypes_xz(i).ne.MPI_INTEGER) then
        call MPI_Type_free(decomp%xtypes_xz(i),ierror)
      endif
    enddo
    deallocate(decomp%xcnts_xz,decomp%zcnts_xz)
    deallocate(decomp%xtypes_xz,decomp%ztypes_xz)
    deallocate(decomp%xdispls_xz)
    deallocate(decomp%zdispls_xz)
#ifdef MPI3
    deallocate(decomp%xranks,decomp%zranks)
#endif

!    deallocate(decomp%x1dist,decomp%y1dist,decomp%y2dist,decomp%z2dist)
!    deallocate(decomp%x1cnts,decomp%y1cnts,decomp%y2cnts,decomp%z2cnts)
!    deallocate(decomp%x1disp,decomp%y1disp,decomp%y2disp,decomp%z2disp)

#ifdef SHM
    deallocate(decomp%x1disp_o,decomp%y1disp_o,decomp%y2disp_o, &
         decomp%z2disp_o)
    deallocate(decomp%x1cnts_s,decomp%y1cnts_s,decomp%y2cnts_s, &
         decomp%z2cnts_s)
    deallocate(decomp%x1disp_s,decomp%y1disp_s,decomp%y2disp_s, &
         decomp%z2disp_s)
#endif

    return
  end subroutine decomp_info_finalize


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Find sub-domain information held by current processor
  !   INPUT: 
  !     nx, ny, nz - global data dimension
  !     pdim(3)    - number of processor grid in each dimension, 
  !                  valid values: 1 - distibute locally; 
  !                                2 - distribute across p_row; 
  !                                3 - distribute across p_col
  !   OUTPUT:
  !     lstart(3)  - starting index
  !     lend(3)    - ending index
  !     lsize(3)   - size of the sub-block (redundant) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine partition(nx, ny, nz, pdim, lstart, lend, lsize)

    implicit none

    integer, intent(IN) :: nx, ny, nz
    integer, dimension(3), intent(IN) :: pdim
    integer, dimension(3), intent(OUT) :: lstart, lend, lsize

    integer, allocatable, dimension(:) :: st,en,sz
    integer :: i, gsize

    do i = 1, 3
 
      if (i==1) then
        gsize = nx
      else if (i==2) then
        gsize = ny
      else if (i==3) then
        gsize = nz
      end if

      if (pdim(i) == 1) then        ! all local
        lstart(i) = 1
        lend(i)   = gsize
        lsize(i)  = gsize
      elseif (pdim(i) == 2) then    ! distribute across dims(1)
        allocate(st(0:dims(1)-1))
        allocate(en(0:dims(1)-1))
        allocate(sz(0:dims(1)-1))
        call distribute(gsize,dims(1),st,en,sz)
        lstart(i) = st(coord(1))
        lend(i)   = en(coord(1))
        lsize(i)  = sz(coord(1))
        deallocate(st,en,sz)
      elseif (pdim(i) == 3) then    ! distribute across dims(2)
        allocate(st(0:dims(2)-1))
        allocate(en(0:dims(2)-1))
        allocate(sz(0:dims(2)-1))
        call distribute(gsize,dims(2),st,en,sz)
        lstart(i) = st(coord(2))
        lend(i)   = en(coord(2))
        lsize(i)  = sz(coord(2))
        deallocate(st,en,sz)
      end if    

    end do
    return   

  end subroutine partition

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   - distibutes grid points in one dimension
  !   - handles uneven distribution properly 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine distribute(data1,proc,st,en,sz)
  
    implicit none
    ! data1 -- data size in any dimension to be partitioned
    ! proc  -- number of processors in that dimension
    ! st    -- array of starting index
    ! en    -- array of ending index
    ! sz    -- array of local size  (redundent)
    integer data1,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
    integer i,size1,nl,nu
  
    size1=data1/proc
    nu = data1 - size1 * proc
    nl = proc - nu
    st(0) = 1
    sz(0) = size1
    en(0) = size1
    do i=1,nl-1
      st(i) = st(i-1) + size1
      sz(i) = size1
      en(i) = en(i-1) + size1
    end do
    size1 = size1 + 1
    do i=nl,proc-1
      st(i) = en(i-1) + 1
      sz(i) = size1
      en(i) = en(i-1) + size1
    end do
    en(proc-1)= data1 
    sz(proc-1)= data1-st(proc-1)+1
  
    return
  end subroutine distribute

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Define how each dimension is distributed across processors
  !    e.g. 17 meshes across 4 processor would be distibuted as (4,4,4,5)
  !    such global information is required locally at MPI_ALLTOALLV time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_dist(nx,ny,nz,decomp)

    integer, intent(IN) :: nx, ny, nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp
!    integer, allocatable, dimension(:) :: st,en

!    allocate(st(0:dims(1)-1))
!    allocate(en(0:dims(1)-1))
!    call distribute(nx,dims(1),st,en,decomp%x1dist)
!    call distribute(ny,dims(1),st,en,decomp%y1dist)
!    deallocate(st,en)

!    allocate(st(0:dims(2)-1))
!    allocate(en(0:dims(2)-1))
!    call distribute(ny,dims(2),st,en,decomp%y2dist)
!    call distribute(nz,dims(2),st,en,decomp%z2dist)
!    deallocate(st,en)
    call distribute(nx,dims(1),decomp%x1st,decomp%x1en,decomp%x1dist)
    call distribute(ny,dims(1),decomp%y1st,decomp%y1en,decomp%y1dist)

    call distribute(ny,dims(2),decomp%y2st,decomp%y2en,decomp%y2dist)
    call distribute(nz,dims(2),decomp%z2st,decomp%z2en,decomp%z2dist)

    return
  end subroutine get_dist

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Prepare the send / receive buffers for MPI_ALLTOALLV communications
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine prepare_buffer(decomp)
    
    implicit none
    
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp
!
    integer :: i, k
    integer :: rank_x, rank_z
    integer :: subsize_y, offset_y
    integer :: ierror
#ifdef MPI3
    integer,dimension(nproc) :: xranks,zranks
    integer,dimension(nproc) :: xweights,zweights
    integer :: index_src, index_dest
#endif

!    integer :: i
!
!    ! MPI_ALLTOALLV buffer information
    do i=0, dims(1)-1
       decomp%x1cnts(i) = decomp%x1dist(i)*decomp%xsz(2)*decomp%xsz(3)
       decomp%y1cnts(i) = decomp%ysz(1)*decomp%y1dist(i)*decomp%ysz(3)
       decomp%z1cnts(i) = decomp%zsz(1)*decomp%zsz(2)*decomp%z1dist(i)
       if (i==0) then
          decomp%x1disp(i) = 0  ! displacement is 0-based index
          decomp%y1disp(i) = 0
          decomp%z1disp(i) = 0
       else
          decomp%x1disp(i) = decomp%x1disp(i-1) + decomp%x1cnts(i-1)
          decomp%y1disp(i) = decomp%y1disp(i-1) + decomp%y1cnts(i-1)
          decomp%z1disp(i) = decomp%z1disp(i-1) + decomp%z1cnts(i-1)
       end if
    end do

    do i=0, dims(2)-1
       decomp%x2cnts(i) = decomp%x2dist(i)*decomp%xsz(2)*decomp%xsz(3)
       decomp%y2cnts(i) = decomp%ysz(1)*decomp%y2dist(i)*decomp%ysz(3)
       decomp%z2cnts(i) = decomp%zsz(1)*decomp%zsz(2)*decomp%z2dist(i)
       if (i==0) then
          decomp%x2disp(i) = 0  ! displacement is 0-based index
          decomp%y2disp(i) = 0  ! displacement is 0-based index
          decomp%z2disp(i) = 0
       else
          decomp%x2disp(i) = decomp%x2disp(i-1) + decomp%x2cnts(i-1)
          decomp%y2disp(i) = decomp%y2disp(i-1) + decomp%y2cnts(i-1)
          decomp%z2disp(i) = decomp%z2disp(i-1) + decomp%z2cnts(i-1)
       end if
    end do
!
!    do i=0, dims(1)-1
!       decomp%x1cnts(i) = decomp%x1dist(i)*decomp%xsz(2)*decomp%xsz(3)
!       decomp%y1cnts(i) = decomp%ysz(1)*decomp%y1dist(i)*decomp%ysz(3)
!       if (i==0) then
!          decomp%x1disp(i) = 0  ! displacement is 0-based index
!          decomp%y1disp(i) = 0
!       else
!          decomp%x1disp(i) = decomp%x1disp(i-1) + decomp%x1cnts(i-1)
!          decomp%y1disp(i) = decomp%y1disp(i-1) + decomp%y1cnts(i-1)
!       end if
!    end do
!    
!    do i=0, dims(2)-1
!       decomp%y2cnts(i) = decomp%ysz(1)*decomp%y2dist(i)*decomp%ysz(3)
!       decomp%z2cnts(i) = decomp%zsz(1)*decomp%zsz(2)*decomp%z2dist(i)
!       if (i==0) then
!          decomp%y2disp(i) = 0  ! displacement is 0-based index
!          decomp%z2disp(i) = 0
!       else
!          decomp%y2disp(i) = decomp%y2disp(i-1) + decomp%y2cnts(i-1)
!          decomp%z2disp(i) = decomp%z2disp(i-1) + decomp%z2cnts(i-1)
!       end if
!    end do
    
    ! MPI_ALLTOALL buffer information

    ! For evenly distributed data, following is an easier implementation.
    ! But it should be covered by the more general formulation below.
    !decomp%x1count = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)/dims(1)
    !decomp%y1count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(1) 
    !decomp%y2count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(2)
    !decomp%z2count = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)/dims(2)

    ! For unevenly distributed data, pad smaller messages. Note the 
    ! last blocks along pencils always get assigned more mesh points
    ! for X <=> Y transposes
    decomp%x1count = decomp%x1dist(dims(1)-1) * &
         decomp%y1dist(dims(1)-1) * decomp%xsz(3)
    decomp%y1count = decomp%x1count
    ! for Y <=> Z transposes
    decomp%y2count = decomp%y2dist(dims(2)-1) * &
         decomp%z2dist(dims(2)-1) * decomp%zsz(1)
    decomp%z2count = decomp%y2count
    ! Information for MPI_Alltoallw for complex X <=> Z transposes
    decomp%xdispls_xz(:)=0
    decomp%zdispls_xz(:)=0
    decomp%xcnts_xz(:)=0
    decomp%zcnts_xz(:)=0
    decomp%xtypes_xz(:)=MPI_INTEGER
    decomp%ztypes_xz(:)=MPI_INTEGER
#ifdef MPI3
    index_src=0
    index_dest=0
#endif
    do k=0,dims(1)-1
    do i=0,dims(2)-1
! Actually, rank_x and rank_z are the same..
           call MPI_Cart_rank(DECOMP_2D_COMM_CART_X,(/k,i/),rank_x,ierror)
! this call fails, since DECOMP_2D_COMM_CART_Z is not yet created
!           call MPI_Cart_rank(DECOMP_2D_COMM_CART_Z,(/k,i/),rank_z,ierror)
           rank_z=rank_x
!JD no checks on x- or z-dimension, since we transform from z- into x-pencils, so these always overlap.
           if (decomp%zst(2).le.decomp%y1en(k) .and. &
               decomp%zen(2).ge.decomp%y1st(k)) then
             decomp%zcnts_xz(rank_z+1)=1
             subsize_y=min(decomp%zen(2),decomp%y1en(k))-max(decomp%zst(2),decomp%y1st(k))+1
             offset_y =max(decomp%zst(2),decomp%y1st(k))-decomp%zst(2)
#ifdef MPI3
             index_src=index_src+1
             zranks(index_src)=rank_z
             zweights(index_src)=decomp%zsz(1)*subsize_y*decomp%z2dist(i)
#endif
             call MPI_Type_create_subarray(3,decomp%zsz, &
               (/decomp%zsz(1),subsize_y,decomp%z2dist(i)/), &
               (/0,offset_y,decomp%z2st(i)-decomp%zst(3)/), &
               !MPI_ORDER_FORTRAN,complex_type,decomp%ztypes_xz(rank_z+1),ierror)
               MPI_ORDER_FORTRAN,real_type,decomp%ztypes_xz(rank_z+1),ierror)
             call MPI_Type_commit(decomp%ztypes_xz(rank_z+1),ierror)
!JD send to process with x-pencil defined by (k,i)
!JD x-bounds are taken from the z-pencils
!             send: decomp%zst(1):decomp%zen(1)
!JD y-bounds are the overlapping region of both pencils.
!                   max(decomp%zst(2),decomp%y1st(k)):min(decomp%zen(2),decomp%y1en(k))
!JD z-bounds are taken from the x-pencils.
!                   decomp%z2st(i):decomp%z2en(i)
           endif
!JD no checks on x- or z-dimension, since we transform from z- into x-pencils, so these always overlap.
           if (decomp%xst(2).le.decomp%y2en(i) .and. &
               decomp%xen(2).ge.decomp%y2st(i)) then
             decomp%xcnts_xz(rank_x+1)=1
             subsize_y=min(decomp%xen(2),decomp%y2en(i))-max(decomp%xst(2),decomp%y2st(i))+1
             offset_y =max(decomp%xst(2),decomp%y2st(i))-decomp%xst(2)
#ifdef MPI3
             index_dest=index_dest+1
             xranks(index_dest)=rank_x
             xweights(index_dest)=decomp%x1dist(k)*subsize_y*decomp%xsz(3)
#endif
             call MPI_Type_create_subarray(3,decomp%xsz, &
               (/decomp%x1dist(k),subsize_y,decomp%xsz(3)/), &
               (/decomp%x1st(k)-decomp%xst(1),offset_y,0/), &
               !MPI_ORDER_FORTRAN,complex_type,decomp%xtypes_xz(rank_x+1),ierror)
               MPI_ORDER_FORTRAN,real_type,decomp%xtypes_xz(rank_x+1),ierror)
             call MPI_Type_commit(decomp%xtypes_xz(rank_x+1),ierror)
!JD recv from process with z-pencil defined by (k,i)
!JD x-bounds are taken from the z-pencils
!             send: decomp%x1st(k):decomp%x1en(k)
!JD y-bounds are the overlapping region of both pencils.
!                   max(decomp%xst(2),decomp%y2st(i)):min(decomp%xen(2),decomp%y2en(i))
!JD z-bounds are taken from the x-pencils.
!                   decomp%xst(3):decomp%xen(3)
           endif
         enddo
         enddo

#ifdef MPI3
    allocate(decomp%xranks(index_dest))
    allocate(decomp%zranks(index_src))

    decomp%xranks=xranks(1:index_dest)+1
    decomp%zranks=zranks(1:index_src)+1

    call MPI_Dist_graph_create_adjacent(DECOMP_2D_COMM_CART_X, &
      index_src,zranks(1:index_src),zweights(1:index_src), &
      index_dest,xranks(1:index_dest),xweights(1:index_dest), &
      MPI_INFO_NULL,.true.,decomp%xtozNeighborComm,ierror)
    call MPI_Dist_graph_create_adjacent(DECOMP_2D_COMM_CART_X, &
      index_dest,xranks(1:index_dest),xweights(1:index_dest), &
      index_src,zranks(1:index_src),zweights(1:index_src), &
      MPI_INFO_NULL,.true.,decomp%ztoxNeighborComm,ierror)
#endif
    
    return
  end subroutine prepare_buffer  

#ifdef SHM

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Generate shared-memory information 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_init_shm(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    ! a copy of old displacement array (will be overwritten by shm code)
    allocate(decomp%x1disp_o(0:dims(1)-1),decomp%y1disp_o(0:dims(1)-1), &
         decomp%y2disp_o(0:dims(2)-1),decomp%z2disp_o(0:dims(2)-1))
    decomp%x1disp_o = decomp%x1disp
    decomp%y1disp_o = decomp%y1disp
    decomp%y2disp_o = decomp%y2disp
    decomp%z2disp_o = decomp%z2disp

    call prepare_shared_buffer(decomp%ROW_INFO,DECOMP_2D_COMM_ROW,decomp)
    call prepare_shared_buffer(decomp%COL_INFO,DECOMP_2D_COMM_COL,decomp)

    return
  end subroutine decomp_info_init_shm


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! For shared-memory implementation, prepare send/recv shared buffer
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine prepare_shared_buffer(C,MPI_COMM,decomp)

    implicit none
    
    TYPE(SMP_INFO) :: C
    INTEGER :: MPI_COMM
    TYPE(DECOMP_INFO) :: decomp
    
    INTEGER, ALLOCATABLE :: KTBL(:,:),NARY(:,:),KTBLALL(:,:)
    INTEGER MYSMP, MYCORE, COLOR
    
    integer :: ierror
    
    C%MPI_COMM = MPI_COMM
    CALL MPI_COMM_SIZE(MPI_COMM,C%NCPU,ierror)
    CALL MPI_COMM_RANK(MPI_COMM,C%NODE_ME,ierror)
    C%SMP_COMM  = MPI_COMM_NULL
    C%CORE_COMM = MPI_COMM_NULL
    C%SMP_ME= 0
    C%NCORE = 0
    C%CORE_ME = 0
    C%MAXCORE = 0
    C%NSMP  = 0
    C%N_SND = 0
    C%N_RCV = 0
    C%SND_P = 0
    C%RCV_P = 0
    C%SND_P_c = 0
    C%RCV_P_c = 0
    
    ! get smp-node map for this communicator and set up smp communicators
    CALL GET_SMP_MAP(C%MPI_COMM, C%NSMP, MYSMP, &
         C%NCORE, MYCORE, C%MAXCORE)
    C%SMP_ME = MYSMP + 1
    C%CORE_ME = MYCORE + 1
    ! - set up inter/intra smp-node communicators
    COLOR = MYCORE
    IF (COLOR.GT.0) COLOR = MPI_UNDEFINED
    CALL MPI_Comm_split(C%MPI_COMM, COLOR, MYSMP, C%SMP_COMM, ierror)
    CALL MPI_Comm_split(C%MPI_COMM, MYSMP, MYCORE, C%CORE_COMM, ierror)
    ! - allocate work space
    ALLOCATE(KTBL(C%MAXCORE,C%NSMP),NARY(C%NCPU,C%NCORE))
    ALLOCATE(KTBLALL(C%MAXCORE,C%NSMP))
    ! - set up smp-node/core to node_me lookup table
    KTBL = 0
    KTBL(C%CORE_ME,C%SMP_ME) = C%NODE_ME + 1
    CALL MPI_ALLREDUCE(KTBL,KTBLALL,C%NSMP*C%MAXCORE,MPI_INTEGER, &
         MPI_SUM,MPI_COMM,ierror)
    KTBL=KTBLALL
    !  IF (SUM(KTBL) /= C%NCPU*(C%NCPU+1)/2) &
    !       CALL MPI_ABORT(...
    
    ! compute offsets in shared SNDBUF and RCVBUF
    CALL MAPSET_SMPSHM(C, KTBL, NARY, decomp)
    
    DEALLOCATE(KTBL,NARY)
    
    return
  end subroutine prepare_shared_buffer

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Use Ian Bush's FreeIPC to generate shared-memory information
  !  - system independent solution
  !  - replacing David Tanqueray's implementation in alloc_shm.c
  !    (old C code renamed to get_smp_map2)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_smp_map(comm, nnodes, my_node, ncores, my_core, maxcor)
    
    use FIPC_module
    
    implicit none
    
    integer, intent(IN) :: comm
    integer, intent(OUT) :: nnodes, my_node, ncores, my_core, maxcor
    
    integer :: intra_comm, extra_comm
    integer :: ierror
    
    call FIPC_init(comm, ierror)
    
    ! intra_comm: communicator for processes on this shared memory node
    ! extra_comm: communicator for all rank 0 on each shared memory node
    call FIPC_ctxt_intra_comm(FIPC_ctxt_world, intra_comm, ierror)
    call FIPC_ctxt_extra_comm(FIPC_ctxt_world, extra_comm, ierror)
    
    call MPI_COMM_SIZE(intra_comm,  ncores, ierror)
    call MPI_COMM_RANK(intra_comm, my_core, ierror)
    
    ! only rank 0 on each shared memory node member of extra_comm
    ! for others extra_comm = MPI_COMM_NULL
    if (extra_comm /= MPI_COMM_NULL) then
       call MPI_COMM_SIZE(extra_comm,  nnodes, ierror)
       call MPI_COMM_RANK(extra_comm, my_node, ierror)
    end if
    
    ! other ranks share the same information as their leaders
    call MPI_BCAST( nnodes, 1, MPI_INTEGER, 0, intra_comm, ierror)
    call MPI_BCAST(my_node, 1, MPI_INTEGER, 0, intra_comm, ierror)
    
    ! maxcor
    call MPI_ALLREDUCE(ncores, maxcor, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierror)
    
    call FIPC_finalize(ierror)
    
    return
    
  end subroutine get_smp_map


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set up smp-node based shared memory maps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MAPSET_SMPSHM(C, KTBL, NARY, decomp)
        
    IMPLICIT NONE
    
    TYPE (SMP_INFO) C
    INTEGER KTBL(C%MAXCORE,C%NSMP)
    INTEGER NARY(C%NCPU,C%NCORE)
    TYPE (DECOMP_INFO) :: decomp

    INTEGER i, j, k, l, N, PTR, BSIZ, ierror, status, seed
    character*16 s
 
    BSIZ = C%N_SND
    
    ! a - SNDBUF
    IF (C%MPI_COMM==DECOMP_2D_COMM_COL) THEN
       ALLOCATE(decomp%x1cnts_s(C%NSMP),decomp%x1disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%x1cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       PTR = 0
       DO i=1,C%NSMP
          decomp%x1disp_s(i) = PTR
          N = 0
          DO j=1,C%MAXCORE
             k = KTBL(j,i)
             IF (k > 0) then
                DO l=1,C%NCORE
                   IF (l == C%CORE_ME) decomp%x1disp_o(k-1) = PTR
                   N = N + NARY(k,l)
                   PTR = PTR + NARY(k,l)
                END DO
             END IF
          END DO
          decomp%x1cnts_s(i) = N
       END DO
       decomp%x1disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR
       
    ELSE IF (C%MPI_COMM==DECOMP_2D_COMM_ROW) THEN
       ALLOCATE(decomp%y2cnts_s(C%NSMP),decomp%y2disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%y2cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       PTR = 0
       DO i=1,C%NSMP
          decomp%y2disp_s(i) = PTR
          N = 0
          DO j=1,C%MAXCORE
             k = KTBL(j,i)
             IF (k > 0) then
                DO l=1,C%NCORE
                   IF (l == C%CORE_ME) decomp%y2disp_o(k-1) = PTR
                   N = N + NARY(k,l)
                   PTR = PTR + NARY(k,l)
                END DO
             END IF
          END DO
          decomp%y2cnts_s(i) = N
       END DO
       decomp%y2disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR
    END IF
    
    ! b - RCVBUF
    
    IF (C%MPI_COMM==DECOMP_2D_COMM_COL) THEN
       ALLOCATE(decomp%y1cnts_s(C%NSMP),decomp%y1disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%y1cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       PTR = 0
       DO i=1,C%NSMP
          decomp%y1disp_s(i) = PTR
          N=0
          DO j=1,C%NCORE
             DO l=1,C%MAXCORE
                k = KTBL(l,i)
                IF (k > 0) then
                   IF (j == C%CORE_ME) decomp%y1disp_o(k-1) = PTR
                   N = N + NARY(k,j)
                   PTR = PTR + NARY(k,j)
                END IF
             END DO
          END DO
          decomp%y1cnts_s(i) = N
       END DO
       decomp%y1disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR
       
    ELSE IF (C%MPI_COMM==DECOMP_2D_COMM_ROW) THEN
       ALLOCATE(decomp%z2cnts_s(C%NSMP),decomp%z2disp_s(C%NSMP+1), &
            stat=status)
       CALL MPI_Allgather(decomp%z2cnts, C%NCPU, MPI_INTEGER, &
            NARY, C%NCPU, MPI_INTEGER, C%CORE_COMM, ierror)
       PTR = 0
       DO i=1,C%NSMP
          decomp%z2disp_s(i) = PTR
          N=0
          DO j=1,C%NCORE
             DO l=1,C%MAXCORE
                k = KTBL(l,i)
                IF (k > 0) then
                   IF (j == C%CORE_ME) decomp%z2disp_o(k-1) = PTR
                   N = N + NARY(k,j)
                   PTR = PTR + NARY(k,j)
                END IF
             END DO
          END DO
          decomp%z2cnts_s(i) = N
       END DO
       decomp%z2disp_s(C%NSMP+1) = PTR
       IF (PTR > BSIZ) BSIZ = PTR
       
    END IF
    
    ! check buffer size and (re)-allocate buffer space if necessary
    IF (BSIZ > C%N_SND) then
       IF (C%SND_P /= 0) CALL DEALLOC_SHM(C%SND_P, C%CORE_COMM)
       ! make sure each rank has unique keys to get shared memory
       !IF (C%MPI_COMM==DECOMP_2D_COMM_COL) THEN
       !   seed = nrank+nproc*0+1 ! has to be non-zero
       !ELSE IF (C%MPI_COMM==DECOMP_2D_COMM_ROW) THEN
       !   seed = nrank+nproc*1+1
       !END IF
       status = 1
       !CALL ALLOC_SHM(C%SND_P, BSIZ, real_type, C%CORE_COMM, status, &
       !     seed)
       CALL ALLOC_SHM(C%SND_P, BSIZ, real_type, C%CORE_COMM, status)
       C%N_SND = BSIZ

       IF (C%RCV_P /= 0) CALL DEALLOC_SHM(C%RCV_P, C%CORE_COMM)
       status = 1
       CALL ALLOC_SHM(C%RCV_P, BSIZ, real_type, C%CORE_COMM, status)
       C%N_RCV = BSIZ

       IF (C%SND_P_c /= 0) CALL DEALLOC_SHM(C%SND_P_c, C%CORE_COMM)
       status = 1
       CALL ALLOC_SHM(C%SND_P_c, BSIZ, complex_type, C%CORE_COMM, status)
       C%N_SND = BSIZ

       IF (C%RCV_P_c /= 0) CALL DEALLOC_SHM(C%RCV_P_c, C%CORE_COMM)
       status = 1
       CALL ALLOC_SHM(C%RCV_P_c, BSIZ, complex_type, C%CORE_COMM, status)
       C%N_RCV = BSIZ


    END IF
    
    RETURN
  END SUBROUTINE MAPSET_SMPSHM

#endif


#ifdef GLOBAL_ARRAYS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Create global arrays that mapped to pencil decompisitions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_global_array(ga, ipencil, data_type, opt_decomp)
    
    implicit none

    integer, intent(OUT) :: ga
    integer, intent(IN) :: ipencil ! 1=X-pencil; 2=Y-pencil; 3=Z-pencil
    integer, intent(IN) :: data_type
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: nblock
    integer, allocatable, dimension(:) :: map
    integer :: offset, i, errorcode
    logical :: success

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    ga = ga_create_handle()
    call ga_set_data(ga, 3, &
         (/decomp%xsz(1),decomp%ysz(2),decomp%zsz(3)/), data_type)
    allocate(map(1+dims(1)+dims(2)))

    ! generate the GA irreg distribution parameters using 
    ! 2DECOMP's decomposition information
    if (ipencil==1) then  ! X-pencil
       nblock(1) = 1
       nblock(2) = dims(1)
       nblock(3) = dims(2)
       map(1) = 1
       offset = nblock(1)+1
       do i=0, dims(1)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%y1dist(i-1)
          end if
       end do
       offset = nblock(1) + nblock(2) + 1
       do i=0, dims(2)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%z2dist(i-1)
          end if
       end do
    else if (ipencil==2) then  ! Y-pencil
       nblock(1) = dims(1)
       nblock(2) = 1
       nblock(3) = dims(2)
       offset = 1
       do i=0, dims(1)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%x1dist(i-1)
          end if
       end do
       map(nblock(1)+1) = 1
       offset = nblock(1) + nblock(2) + 1
       do i=0, dims(2)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%z2dist(i-1)
          end if
       end do
    else if (ipencil==3) then  ! Z-pencil
       nblock(1) = dims(1)
       nblock(2) = dims(2)
       nblock(3) = 1
       offset = 1
       do i=0, dims(1)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%x1dist(i-1)
          end if
       end do
       offset = nblock(1)+1
       do i=0, dims(2)-1
          if (i==0) then
             map(offset+i) = 1
          else
             map(offset+i) = map(offset+i-1) + decomp%y2dist(i-1)
          end if
       end do
       map(nblock(1)+nblock(2)+1) = 1
    end if

    call ga_set_irreg_distr(ga, map, nblock)
    success = ga_allocate(ga)
    if (.not.success) then
       errorcode = 7
       call decomp_2d_abort(errorcode, &
            'Failed to create global arrays')
    end if

    deallocate(map)

    return
  end subroutine get_global_array

#endif


#ifdef OCC
  ! For non-blocking communication code, progress the comminication stack
  subroutine transpose_test(handle)

    implicit none

    integer :: handle, ierror

    call NBC_TEST(handle,ierror)

    return
  end subroutine transpose_test
#endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transposition routines 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "transpose_x_to_y.f90"
#include "transpose_y_to_z.f90"
#include "transpose_z_to_y.f90"
#include "transpose_y_to_x.f90"
#include "transpose_x_to_z.f90"
#include "transpose_z_to_x.f90"


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Auto-tuning algorithm to select the best 2D processor grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine best_2d_grid(iproc, best_p_row, best_p_col)

    implicit none

    integer, intent(IN) :: iproc
    integer, intent(OUT) :: best_p_row, best_p_col

    integer, allocatable, dimension(:) :: factors
    double precision :: t1, t2, best_time
    integer :: nfact, i, row, col, ierror, errorcode

    real(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3

    TYPE(DECOMP_INFO) :: decomp

    if (nrank==0) write(*,*) 'In auto-tuning mode......'

    best_time = huge(t1)
    best_p_row = -1
    best_p_col = -1
    
    i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i))
    call findfactor(iproc, factors, nfact)
    if (nrank==0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

!    Make initial communication to un-bias results

    row = factors(1)
    col = iproc / row

       if (min(nx_global,ny_global)>=row .and. &
            min(ny_global,nz_global)>=col) then

          ! 2D Catersian topology
          dims(1) = row
          dims(2) = col
          periodic(1) = .false.
          periodic(2) = .false.
          call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
               .false.,DECOMP_2D_COMM_CART_X, ierror)
          call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)

          ! communicators defining sub-groups for ALLTOALL(V)
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
               DECOMP_2D_COMM_COL,ierror)
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
               DECOMP_2D_COMM_ROW,ierror)

          ! generate 2D decomposition information for this row*col
          call decomp_info_init(nx_global,ny_global,nz_global,decomp)

          ! arrays for X,Y and Z-pencils
          allocate(u1(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
          allocate(u2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))

          ! timing the transposition routines
          t1 = MPI_WTIME()
          call transpose_x_to_y(u1,u2,decomp)
          call transpose_y_to_x(u2,u1,decomp)
          t2 = MPI_WTIME() - t1

          deallocate(u1,u2)
          call decomp_info_finalize(decomp)
          call MPI_Comm_free(DECOMP_2D_COMM_COL,ierror)
          call MPI_Comm_free(DECOMP_2D_COMM_ROW,ierror)
          call MPI_Comm_free(DECOMP_2D_COMM_CART_X, ierror)

       end if


    do i=1, nfact

       row = factors(i)
       col = iproc / row

       ! enforce the limitation of 2D decomposition
       if (min(nx_global,ny_global)>=row .and. &
            min(ny_global,nz_global)>=col) then

          ! 2D Catersian topology
          dims(1) = row
          dims(2) = col
          periodic(1) = .false.
          periodic(2) = .false.
          call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
               .false.,DECOMP_2D_COMM_CART_X, ierror)
          call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
          
          ! communicators defining sub-groups for ALLTOALL(V)
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
               DECOMP_2D_COMM_COL,ierror)
          call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
               DECOMP_2D_COMM_ROW,ierror)
          
          ! generate 2D decomposition information for this row*col
          call decomp_info_init(nx_global,ny_global,nz_global,decomp)

          ! arrays for X,Y and Z-pencils
          allocate(u1(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
          allocate(u2(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
          allocate(u3(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))

          ! timing the transposition routines
          t1 = MPI_WTIME()
          call transpose_x_to_y(u1,u2,decomp)
          call transpose_y_to_z(u2,u3,decomp)
          call transpose_z_to_y(u3,u2,decomp)
          call transpose_y_to_x(u2,u1,decomp)
          t2 = MPI_WTIME() - t1

          deallocate(u1,u2,u3)
          call decomp_info_finalize(decomp)
          call MPI_Comm_free(DECOMP_2D_COMM_COL,ierror)
          call MPI_Comm_free(DECOMP_2D_COMM_ROW,ierror)
          call MPI_Comm_free(DECOMP_2D_COMM_CART_X, ierror)

          call MPI_ALLREDUCE(t2,t1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                   MPI_COMM_WORLD,ierror)
          t1 = t1 / dble(nproc)

          if (nrank==0) then
             write(*,*) 'processor grid', row, ' by ', col, ' time=', t1
          end if

          if (best_time > t1) then
             best_time = t1
             best_p_row = row
             best_p_col = col
          end if

       end if
       
    end do ! loop through processer grid

    deallocate(factors)

    if (best_p_row/=-1) then
       if (nrank==0) then
          write(*,*) 'the best processor grid is probably ', &
               best_p_row, ' by ', best_p_col
       end if
    else
       errorcode = 9
       call decomp_2d_abort(errorcode, &
            'The processor-grid auto-tuning code failed. ' // &
            'The number of processes requested is probably too large.')
    end if

    return
  end subroutine best_2d_grid

#include "factor.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "halo.f90"


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Error handling
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_abort(errorcode, msg)

    implicit none

    integer, intent(IN) :: errorcode
    character(len=*), intent(IN) :: msg

    integer :: ierror
    
    if (nrank==0) then
       write(*,*) '2DECOMP&FFT ERROR - errorcode: ', errorcode
       write(*,*) 'ERROR MESSAGE: ' // msg
    end if
    call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)

    return
  end subroutine decomp_2d_abort


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routines to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "alloc.f90"
    
  
end module decomp_2d

