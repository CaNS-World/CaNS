module diezdecomp_api_cans
  use mpi
  use, intrinsic :: iso_fortran_env, only: i8 => int64, sp => real32, dp => real64
#if defined(_OPENACC)
  use openacc, only: acc_handle_kind
#else
  use, intrinsic :: iso_fortran_env, only: acc_handle_kind => int64
#endif
  use diezdecomp_core
  implicit none
#if defined(_DIEZDECOMP_SINGLE)
  integer, parameter :: rp = sp
  integer, parameter :: MPI_REAL_RP = MPI_REAL
#else
  integer, parameter :: rp = dp
  integer, parameter :: MPI_REAL_RP = MPI_DOUBLE_PRECISION
#endif

  private

  ! --------------- public variables ---------------
  public :: diezdecompPencilInfo                    , &
            diezdecompHandle                        , &
            diezdecompGridDesc                      , &
            diezdecomp_boilerplate_transpose        , &
            diezdecomp_boilerplate_halos            , &
            diezdecomp_get_workspace_size_transposes, &
            diezdecomp_get_workspace_size_halos     , &
            diezdecomp_rank_null                    , &
            diezdecompTransposeXtoY                 , &
            diezdecompTransposeYtoX                 , &
            diezdecompTransposeYtoZ                 , &
            diezdecompTransposeZtoY                 , &
            diezdecompTransposeXtoZ                 , &
            diezdecompTransposeZtoX                 , &
            diezdecompUpdateHalosX                  , &
            diezdecompUpdateHalosY                  , &
            diezdecompUpdateHalosZ                  , &
            diezdecompGetShiftedRank                , &
            diezdecompGridDescCreate                , &
            diezdecomp_halo_mpi_mode                , &
            diezdecompGetPencilInfo                 , &
            diezdecomp_autotune_mode_trials         , &
            diezdecomp_n_warmup                     , &
            diezdecomp_halo_autotuned_pack          , &
            diezdecomp_allow_autotune_reorder       , &
            diezdecomp_summary_halo_autotuning      , &
            diezdecomp_summary_transp_autotuning

  ! these parameters can be changed to control internal behavior (beyond the API)
  integer, save            :: diezdecomp_rank_null              = -1
  logical, save            :: diezdecomp_enable_alltoallv       = .false.
  integer, save            :: diezdecomp_halo_mpi_mode          = 2
  integer, save            :: diezdecomp_halo_autotuned_pack    = 2
  logical, save            :: diezdecomp_allow_autotune_reorder = .false.

  ! --------------- types ---------------
  type diezdecompPencilInfo
    integer              :: lo(3), hi(3), shape(3), internal_order(3), internal_invorder(3), pidx(2), order_halo(3), &
                            shape_mpi_ranks(3), size, offset6(0:2,0:1)
    logical              :: is_bound(0:1,3)
    integer, allocatable :: mpi_ranks(:,:,:), flat_mpi_ranks(:,:)
  end type

  type diezdecompHandle
    logical :: placeholder ! might be used later
  end type

  type diezdecompGridDesc
    type(diezdecompPencilInfo)    :: all_ap(0:2)
    integer                       :: irank, nproc
    type(diezdecomp_props_transp) :: obj_tr(0:2,0:2)
    integer                       :: abs_reorder(0:2,0:2,0:2) = -1
  end type

  contains

  ! ---------------------------------------- transpose halos -------------------------------------------
  function diezdecompTransposeXtoY(ch,gd,pa,pb,work,dtype,&
                                   input_halo_extents, output_halo_extents, input_padding, output_padding, stream) result(ierr)
    implicit none
    type(diezdecompHandle)             :: ch
    integer                            :: dtype, ierr, i_halo(0:2),o_halo(0:2),i_pad(0:2),o_pad(0:2)
    type(diezdecompGridDesc)           :: gd
    real(rp)                           :: pa(:,:,:), pb(:,:,:), work(:)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    integer, optional                  :: input_halo_extents(0:2)  ,&
                                          output_halo_extents(0:2) ,&
                                          input_padding(0:2)       ,&
                                          output_padding(0:2)
    logical                            :: is_device_synchronize
    stream_internal       = diezdecomp_stream_default
    is_device_synchronize = (.not.present(stream))
    if (present(stream)) stream_internal = stream
    if (present(input_halo_extents )) then ; i_halo = input_halo_extents  ; else ; i_halo = -1 ; end if
    if (present(output_halo_extents)) then ; o_halo = output_halo_extents ; else ; o_halo = -1 ; end if
    if (present(input_padding      )) then ; i_pad  = input_padding       ; else ; i_pad  = -1 ; end if
    if (present(output_padding     )) then ; o_pad  = output_padding      ; else ; o_pad  = -1 ; end if
    call diezdecomp_boilerplate_transpose(gd, pa, pb, work, 0, 1, diezdecomp_enable_alltoallv, stream_internal, &
                                          is_device_synchronize, i_halo, o_halo, i_pad, o_pad)
    ierr = 0
  end function

  function diezdecompTransposeYtoX(ch,gd,pa,pb,work,dtype,&
                                   input_halo_extents, output_halo_extents, input_padding, output_padding, stream) result(ierr)
    implicit none
    type(diezdecompHandle)             :: ch
    integer                            :: dtype, ierr, i_halo(0:2),o_halo(0:2),i_pad(0:2),o_pad(0:2)
    type(diezdecompGridDesc)           :: gd
    real(rp)                           :: pa(:,:,:), pb(:,:,:), work(:)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    integer, optional                  :: input_halo_extents(0:2)  ,&
                                          output_halo_extents(0:2) ,&
                                          input_padding(0:2)       ,&
                                          output_padding(0:2)
    logical                            :: is_device_synchronize
    stream_internal       = diezdecomp_stream_default
    is_device_synchronize = (.not.present(stream))
    if (present(stream)) stream_internal = stream
    if (present(input_halo_extents )) then ; i_halo = input_halo_extents  ; else ; i_halo = -1 ; end if
    if (present(output_halo_extents)) then ; o_halo = output_halo_extents ; else ; o_halo = -1 ; end if
    if (present(input_padding      )) then ; i_pad  = input_padding       ; else ; i_pad  = -1 ; end if
    if (present(output_padding     )) then ; o_pad  = output_padding      ; else ; o_pad  = -1 ; end if
    call diezdecomp_boilerplate_transpose(gd, pa, pb, work, 1, 0, diezdecomp_enable_alltoallv, stream_internal, &
                                          is_device_synchronize, i_halo, o_halo, i_pad, o_pad)
    ierr = 0
  end function

  function diezdecompTransposeYtoZ(ch,gd,pa,pb,work,dtype,&
                                   input_halo_extents, output_halo_extents, input_padding, output_padding, stream) result(ierr)
    implicit none
    type(diezdecompHandle)             :: ch
    integer                            :: dtype, ierr, i_halo(0:2),o_halo(0:2),i_pad(0:2),o_pad(0:2)
    type(diezdecompGridDesc)           :: gd
    real(rp)                           :: pa(:,:,:), pb(:,:,:), work(:)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    integer, optional                  :: input_halo_extents(0:2)  ,&
                                          output_halo_extents(0:2) ,&
                                          input_padding(0:2)       ,&
                                          output_padding(0:2)
    logical                            :: is_device_synchronize
    stream_internal       = diezdecomp_stream_default
    is_device_synchronize = (.not.present(stream))
    if (present(stream)) stream_internal = stream
    if (present(input_halo_extents )) then ; i_halo = input_halo_extents  ; else ; i_halo = -1 ; end if
    if (present(output_halo_extents)) then ; o_halo = output_halo_extents ; else ; o_halo = -1 ; end if
    if (present(input_padding      )) then ; i_pad  = input_padding       ; else ; i_pad  = -1 ; end if
    if (present(output_padding     )) then ; o_pad  = output_padding      ; else ; o_pad  = -1 ; end if
    call diezdecomp_boilerplate_transpose(gd, pa, pb, work, 1, 2, diezdecomp_enable_alltoallv, stream_internal, &
                                          is_device_synchronize, i_halo, o_halo, i_pad, o_pad)
    ierr = 0
  end function

  function diezdecompTransposeZtoY(ch,gd,pa,pb,work,dtype,&
                                   input_halo_extents, output_halo_extents, input_padding, output_padding, stream) result(ierr)
    implicit none
    type(diezdecompHandle)             :: ch
    integer                            :: dtype, ierr, i_halo(0:2),o_halo(0:2),i_pad(0:2),o_pad(0:2)
    type(diezdecompGridDesc)           :: gd
    real(rp)                           :: pa(:,:,:), pb(:,:,:), work(:)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    integer, optional                  :: input_halo_extents(0:2)  ,&
                                          output_halo_extents(0:2) ,&
                                          input_padding(0:2)       ,&
                                          output_padding(0:2)
    logical                            :: is_device_synchronize
    stream_internal       = diezdecomp_stream_default
    is_device_synchronize = (.not.present(stream))
    if (present(stream)) stream_internal = stream
    if (present(input_halo_extents )) then ; i_halo = input_halo_extents  ; else ; i_halo = -1 ; end if
    if (present(output_halo_extents)) then ; o_halo = output_halo_extents ; else ; o_halo = -1 ; end if
    if (present(input_padding      )) then ; i_pad  = input_padding       ; else ; i_pad  = -1 ; end if
    if (present(output_padding     )) then ; o_pad  = output_padding      ; else ; o_pad  = -1 ; end if
    call diezdecomp_boilerplate_transpose(gd, pa, pb, work, 2, 1, diezdecomp_enable_alltoallv, stream_internal, &
                                          is_device_synchronize, i_halo, o_halo, i_pad, o_pad)
    ierr = 0
  end function

  function diezdecompTransposeXtoZ(ch,gd,pa,pb,work,dtype,&
                                   input_halo_extents, output_halo_extents, input_padding, output_padding, stream) result(ierr)
    implicit none
    type(diezdecompHandle)             :: ch
    integer                            :: dtype, ierr, i_halo(0:2),o_halo(0:2),i_pad(0:2),o_pad(0:2)
    type(diezdecompGridDesc)           :: gd
    real(rp)                           :: pa(:,:,:), pb(:,:,:), work(:)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    integer, optional                  :: input_halo_extents(0:2)  ,&
                                          output_halo_extents(0:2) ,&
                                          input_padding(0:2)       ,&
                                          output_padding(0:2)
    logical                            :: is_device_synchronize
    stream_internal       = diezdecomp_stream_default
    is_device_synchronize = (.not.present(stream))
    if (present(stream)) stream_internal = stream
    if (present(input_halo_extents )) then ; i_halo = input_halo_extents  ; else ; i_halo = -1 ; end if
    if (present(output_halo_extents)) then ; o_halo = output_halo_extents ; else ; o_halo = -1 ; end if
    if (present(input_padding      )) then ; i_pad  = input_padding       ; else ; i_pad  = -1 ; end if
    if (present(output_padding     )) then ; o_pad  = output_padding      ; else ; o_pad  = -1 ; end if
    call diezdecomp_boilerplate_transpose(gd, pa, pb, work, 0, 2,.false., stream_internal, &
                                          is_device_synchronize, i_halo, o_halo, i_pad, o_pad)
    ierr = 0
  end function

  function diezdecompTransposeZtoX(ch,gd,pa,pb,work,dtype,&
                                   input_halo_extents, output_halo_extents, input_padding, output_padding, stream) result(ierr)
    implicit none
    type(diezdecompHandle)             :: ch
    integer                            :: dtype, ierr, i_halo(0:2),o_halo(0:2),i_pad(0:2),o_pad(0:2)
    type(diezdecompGridDesc)           :: gd
    real(rp)                           :: pa(:,:,:), pb(:,:,:), work(:)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    integer, optional                  :: input_halo_extents(0:2)  ,&
                                          output_halo_extents(0:2) ,&
                                          input_padding(0:2)       ,&
                                          output_padding(0:2)
    logical                            :: is_device_synchronize
    stream_internal       = diezdecomp_stream_default
    is_device_synchronize = (.not.present(stream))
    if (present(stream)) stream_internal = stream
    if (present(input_halo_extents )) then ; i_halo = input_halo_extents  ; else ; i_halo = -1 ; end if
    if (present(output_halo_extents)) then ; o_halo = output_halo_extents ; else ; o_halo = -1 ; end if
    if (present(input_padding      )) then ; i_pad  = input_padding       ; else ; i_pad  = -1 ; end if
    if (present(output_padding     )) then ; o_pad  = output_padding      ; else ; o_pad  = -1 ; end if
    call diezdecomp_boilerplate_transpose(gd, pa, pb, work, 2, 0,.false., stream_internal, &
                                          is_device_synchronize, i_halo, o_halo, i_pad, o_pad)
    ierr = 0
  end function

  ! ---------------------------------------- api halos -------------------------------------------
  function diezdecompUpdateHalosX(ch,gd,p,work,dtype,nh_xyz,periods,ii,offset6_orig, stream) result(ierr)
    implicit none
    type(diezdecompHandle)             :: ch
    type(diezdecompGridDesc)           :: gd
    real(rp), contiguous               :: p(:,:,:), work(:)
    integer                            :: dtype, nh_xyz(3), ii, ierr,offset6(0:2,0:1)
    logical                            :: periods(3)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    logical                            :: is_device_synchronize
    integer, optional                  :: offset6_orig(0:2,0:1)
    if (present(offset6_orig)) then  ; offset6 =  offset6_orig
    else                             ; offset6 = 0                 ; end if
    stream_internal       = diezdecomp_stream_default
    is_device_synchronize = (.not.present(stream))
    if (present(stream)) stream_internal = stream
    call diezdecomp_boilerplate_halos(gd, p, work, nh_xyz, offset6, periods, ii, 1, stream_internal, is_device_synchronize)
    ierr = 0
  end function

  function diezdecompUpdateHalosY(ch,gd,p,work,dtype,nh_xyz,periods,ii,offset6_orig, stream) result(ierr)
    implicit none
    type(diezdecompHandle)             :: ch
    type(diezdecompGridDesc)           :: gd
    real(rp), contiguous               :: p(:,:,:), work(:)
    integer                            :: dtype, nh_xyz(3), ii, ierr,offset6(0:2,0:1)
    logical                            :: periods(3)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    logical                            :: is_device_synchronize
    integer, optional                  :: offset6_orig(0:2,0:1)
    if (present(offset6_orig)) then  ; offset6 =  offset6_orig
    else                             ; offset6 = 0                 ; end if
    stream_internal       = diezdecomp_stream_default
    is_device_synchronize = (.not.present(stream))
    if (present(stream)) stream_internal = stream
    call diezdecomp_boilerplate_halos(gd, p, work, nh_xyz, offset6, periods, ii, 2, stream_internal, is_device_synchronize)
    ierr = 0
  end function

  function diezdecompUpdateHalosZ(ch,gd,p,work,dtype,nh_xyz,periods,ii,offset6_orig, stream) result(ierr)
    implicit none
    type(diezdecompHandle)             :: ch
    type(diezdecompGridDesc)           :: gd
    real(rp), contiguous               :: p(:,:,:), work(:)
    integer                            :: dtype, nh_xyz(3), ii, ierr,offset6(0:2,0:1)
    logical                            :: periods(3)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    logical                            :: is_device_synchronize
    integer, optional                  :: offset6_orig(0:2,0:1)
    if (present(offset6_orig)) then  ; offset6 =  offset6_orig
    else                             ; offset6 = 0                 ; end if
    stream_internal       = diezdecomp_stream_default
    is_device_synchronize = (.not.present(stream))
    if (present(stream)) stream_internal = stream
    call diezdecomp_boilerplate_halos(gd, p, work, nh_xyz, offset6, periods, ii, 3, stream_internal, is_device_synchronize)
    ierr = 0
  end function

  ! ---------------------------------------- diezdecomp size calculation -------------------------------------------
  function diezdecomp_get_workspace_size_halos(ch,gd,ipencil,nh_xyz,wsize) result(ierr)
    implicit none
    type(diezdecompHandle)   :: ch
    type(diezdecompGridDesc) :: gd
    integer(i8)              :: wsize
    integer                  :: ii, jj, kk, side_numel, ipencil, nh_xyz(0:2), internal_order(0:2),ierr
    ierr  = 0
    wsize = 1
    do ii=0,2
      if      (ii == 0) then ; jj = 1; kk = 2
      else if (ii == 1) then ; jj = 0; kk = 2
      else if (ii == 2) then ; jj = 0; kk = 1 ; end if
      internal_order(0:2) =  gd%all_ap(ipencil-1)%internal_order(1:3)
      side_numel          = (gd%all_ap(ipencil-1)%shape(diezdecomp_get_index_int2int(internal_order,jj)+1) + 2*nh_xyz(jj))* &
                            (gd%all_ap(ipencil-1)%shape(diezdecomp_get_index_int2int(internal_order,kk)+1) + 2*nh_xyz(kk))* &
                            nh_xyz(ii)
      wsize               = max(wsize, int(4*side_numel,i8))
    end do
  ierr = 0
  end function

  function diezdecomp_get_workspace_size_transposes(ch, gd, wsize) result(ierr)
    implicit none
    type(diezdecompHandle)   :: ch
    type(diezdecompGridDesc) :: gd
    integer(i8)              :: wsize, sx, sz, sy
    integer                  :: ierr
    ierr  = 0
    wsize = 1
    sx    = int(gd%all_ap(0)%size,i8)
    sy    = int(gd%all_ap(1)%size,i8)
    sz    = int(gd%all_ap(2)%size,i8)
    wsize = max(wsize, sx+sy, sx+sz, sy+sz)
  end function

  ! ---------------------------------------- boilerplate transpose/halos -------------------------------------------
  subroutine diezdecomp_boilerplate_transpose(gd, p_in, p_out, work, ii, jj, &
                                              allow_alltoallv, stream, is_device_synchronize, i_halo, o_halo, i_pad, o_pad)
    implicit none
    type(diezdecompGridDesc) :: gd
    real(rp)                 :: p_in(:,:,:), p_out(:,:,:), work(:)
    integer                  :: ii, jj, kk, ii_jj_kk(0:2), abs_reorder(0:2), &
                                offset6_in(0:2,0:1), offset6_out(0:2,0:1), i_halo(0:2), o_halo(0:2), i_pad(0:2), o_pad(0:2), &
                                zeros6(0:2,0:1)
    logical                  :: allow_alltoallv, is_device_synchronize
    integer(acc_handle_kind) :: stream
    if (.not.(gd%obj_tr(ii,jj)%initialized)) then
      kk       = 3 - ii - jj
      ii_jj_kk = [ii,jj,kk]
      zeros6   = 0
      abs_reorder = gd%abs_reorder(ii,jj,:)
      call diezdecomp_fill_transp_props(gd%obj_tr(ii,jj)                            , &
                                        abs_reorder                                 , &
                                        gd%all_ap(ii)%internal_order                , &
                                        gd%all_ap(jj)%internal_order                , &
                                        gd%all_ap(ii)%lo - 1, gd%all_ap(ii)%hi - 1  , &
                                        gd%all_ap(jj)%lo - 1, gd%all_ap(jj)%hi - 1  , &
                                        ii_jj_kk                                    , &
                                        gd%irank, gd%nproc,allow_alltoallv          , &
                                        gd%all_ap(ii)%shape, zeros6                 , &
                                        gd%all_ap(jj)%shape, zeros6                 , &
                                        diezdecomp_allow_autotune_reorder, stream   )
    end if
    offset6_in  = gd%obj_tr(ii,jj)%offset6_in
    offset6_out = gd%obj_tr(ii,jj)%offset6_out
    if ((i_halo(0)>=0).and.(i_pad(0)>=0)) then
      gd%obj_tr(ii,jj)%offset6_in(:,0)  =  i_halo
      gd%obj_tr(ii,jj)%offset6_in(:,1)  =  i_halo + i_pad
    end if
    if ((o_halo(0)>=0).and.(o_pad(0)>=0)) then
      gd%obj_tr(ii,jj)%offset6_out(:,0)  =  o_halo
      gd%obj_tr(ii,jj)%offset6_out(:,1)  =  o_halo + o_pad
    end if
    call diezdecomp_transp_execute(gd%obj_tr(ii,jj), p_in, p_out, stream, work)
    gd%obj_tr(ii,jj)%offset6_in  = offset6_in
    gd%obj_tr(ii,jj)%offset6_out = offset6_out
    if (is_device_synchronize) then
      !$acc wait(stream)
    end if
  end subroutine

  subroutine diezdecomp_boilerplate_halos(gd, p, work, nh_xyz, offset6, periods, ii, axis, stream, is_device_synchronize, &
                                          print_summary)
    implicit none
    type(diezdecompGridDesc)    :: gd
    real(rp)                    :: p(:,:,:), work(:)
    integer                     :: nh_xyz(3), ii, axis, p_shape(0:2), offset6(0:2,0:1)
    logical                     :: periods(3), is_per, is_device_synchronize
    logical, optional           :: print_summary
    type(diezdecomp_props_halo) :: this
    integer(acc_handle_kind)    :: stream
    ! call check_bounds(ii  , 1, 3)
    ! call check_bounds(axis, 1, 3)

    if (.not.(this%initialized)) then
      is_per  = periods(ii)
      p_shape = [size(p,1), size(p,2), size(p,3)]
      if (axis==0) error stop 'ERROR: axis out-of-scope' ! , axis
      call diezdecomp_fill_halo_props(this, ii-1, is_per, &
                                      gd%all_ap(axis-1)%mpi_ranks, &
                                      gd%all_ap(axis-1)%flat_mpi_ranks, &
                                      gd%all_ap(axis-1)%shape_mpi_ranks, &
                                      nh_xyz, gd%all_ap(axis-1)%order_halo, &
                                      p_shape, offset6, gd%irank, gd%nproc, diezdecomp_halo_mpi_mode, &
                                      diezdecomp_halo_autotuned_pack)
    end if
    call diezdecomp_halos_execute(this, p, work, stream)
    if (is_device_synchronize) then
      !$acc wait(stream)
    end if
    if (present(print_summary)) then
      if (print_summary) call diezdecomp_summary_halo_autotuning(this)
    end if
  end subroutine

  ! ---------------------------------------- diezdecomp wrappers ---------------------------------------------
  subroutine diezdecompGridDescCreate(gd, pdims, gdims, gdims_dist, transpose_axis_contiguous, periods, ipencil, offset6_xyz_orig)
    implicit none
    type(diezdecompGridDesc) :: gd
    integer                  :: pdims(0:1), gdims(0:2), gdims_dist(0:2), ipencil, axis, mpi_ierr, offset6_xyz(0:2,0:1)
    integer, optional        :: offset6_xyz_orig(0:2,0:1)
    logical                  :: transpose_axis_contiguous(0:2)
    logical                  :: periods(3)
    if (present(offset6_xyz_orig)) then  ; offset6_xyz =  offset6_xyz_orig
    else                                 ; offset6_xyz = 0                 ; end if
    ! call check_bounds(ipencil, 1, 3, 103)

    call mpi_comm_rank(mpi_comm_world, gd%irank, mpi_ierr)
    call mpi_comm_size(mpi_comm_world, gd%nproc, mpi_ierr)
    do axis=0,2
      call summary_cudecompGetPencilInfo(gd%all_ap(axis), axis, pdims, gdims, gdims_dist, transpose_axis_contiguous, periods, &
                                        ipencil, offset6_xyz, gd%irank, gd%nproc)
    end do
    block 
      integer :: i,j,k
      do     i=0,2
        do   j=0,2
          do k=0,2
            ! default initialization: gd%abs_reorder matches output "internal_order"
            gd%abs_reorder(i,j,k) = gd%all_ap(j)%internal_order(k+1)
          end do 
        end do 
      end do 
    end block 
  end subroutine

  function diezdecompGetPencilInfo(ch, gd ,ap, axis)  result(ierr)
    implicit none
    type(diezdecompHandle)   :: ch
    type(diezdecompGridDesc) :: gd
    type(diezdecompPencilInfo) :: ap
    integer                  :: axis, ierr
    ! call check_bounds(axis, 1, 3, 121)
    ap = gd%all_ap(axis-1)
    ierr = 0
  end function

  function diezdecompGetShiftedRank(ch, gd, axis, dim, d, is_per, nb_val)  result(ierr)
    implicit none
    type(diezdecompHandle)   :: ch
    type(diezdecompGridDesc) :: gd
    integer                  :: axis, dim, d, nb_val, ierr
    logical                  :: is_per
    ! call check_bounds(axis,  1, 3, 131)
    ! call check_bounds(dim ,  1, 3, 132)
    ! call check_bounds(d   , -1, 1, 133)
    if (gd%all_ap(axis-1)%is_bound(max(0,d),dim)) then
      nb_val = diezdecomp_rank_null
    else
      nb_val = diezdecomp_get_rank_id(gd%all_ap(axis-1)%mpi_ranks, &
                                      gd%all_ap(axis-1)%flat_mpi_ranks, &
                                      gd%all_ap(axis-1)%shape_mpi_ranks, gd%irank, dim-1, d)
    endif
    ierr = 0
  end function

  subroutine diezdecomp_fill_is_bound(this, flat_mpi_ranks, shape_mpi_ranks, periods, ipencil, irank)
    implicit none
    type(diezdecompPencilInfo) :: this
    integer                  :: irank, i0, j0, k0, i, ind, ipencil
    integer                  :: shape_mpi_ranks(0:2), flat_mpi_ranks(0:,0:)
    logical                  :: periods(3)

    ! call check_bounds(ipencil, 1, 3)
    i0 = flat_mpi_ranks(irank,0) + 1
    j0 = flat_mpi_ranks(irank,1) + 1
    k0 = flat_mpi_ranks(irank,2) + 1

    this%is_bound(:,:) = .false.
    do i=1,3
      if (i==1) ind = i0
      if (i==2) ind = j0
      if (i==3) ind = k0
      if (.not.periods(i)) then
        if (ind == 1                   ) this%is_bound(0,i) = .true.
        if (ind == shape_mpi_ranks(i-1)) this%is_bound(1,i) = .true.
      end if
      if (i == ipencil) then
        this%is_bound(0,i) = .true.
        this%is_bound(1,i) = .true.
      end if
    end do
  end subroutine

  ! ---------------------------------------- ported routines -------------------------------------------
  subroutine summary_cudecompGetPencilInfo(this, axis, pdims, gdims, gdims_dist, transpose_axis_contiguous, periods, ipencil, &
                                          offset6_xyz, irank, nproc)
    !  (Needed for compatibility with CaNS)
    !  Summary of external behavior in:
    !    cuDecomp: An Adaptive Pencil Decomposition Library for NVIDIA GPUs
    !    https://github.com/NVIDIA/cuDecomp/blob/main/src/cudecomp.cc
    implicit none
    integer                  :: lo(0:2), hi(0:2), shape(0:2), order(0:2), invorder(0:2), pidx(0:1), &
                                axis, pdims(0:1), gdims(0:2), gdims_dist(0:2), offset6_xyz(0:2,0:1)
    logical                  :: transpose_axis_contiguous(0:2)
    integer                  :: i,j,k,d,odd,irank,nproc
    type(diezdecompPencilInfo) :: this
    ! last 3 arguments (used at the end of this routine for simplicity)
    integer                  :: ipencil
    logical                  :: periods(3)

    ! call check_bounds(axis   , 0,       2)
    ! call check_bounds(ipencil, 1,       3)
    ! call check_bounds(irank  , 0, nproc-1)

    pidx(0) = irank / pdims(1)
    pidx(1) = mod(irank, pdims(1))

    do i=0,2
      if (transpose_axis_contiguous(axis)) then
        order(i) = mod(axis + i, 3)
      else
        order(i) = i
      endif
        invorder(order(i)) = i
      end do
    j = 0
    do i=0,2
      k = invorder(i)
      if (i /= axis) then
        d   = gdims_dist(i) / pdims(j)
        odd = mod(gdims_dist(i), pdims(j))
        shape(k) = d
        if (pidx(j) < odd)                               shape(k) = shape(k) + 1
        if (pidx(j) == min(pdims(j), gdims_dist(i)) - 1) shape(k) = shape(k) + gdims(i) - gdims_dist(i)
        lo(k) = pidx(j) * d + min(pidx(j), odd)
        j       = j + 1
      else
        shape(k) = gdims(i)
        lo(k)    = 0
      end if
      hi(k)      = lo(k) + shape(k) - 1
    end do
    ! copy data (with shift)
    do i=0,2
      this%lo               (i+1) = lo      (i) + 1
      this%hi               (i+1) = hi      (i) + 1
      this%shape            (i+1) = shape   (i)
      this%internal_order   (i+1) = order   (i)
      this%internal_invorder(i+1) = invorder(i)
    end do
    this%pidx(1) = pidx(0)
    this%pidx(2) = pidx(1)
    ! ------ diezDecomp objects ------
    ! fill is_bound and mpi_ranks
    call diezdecomp_fill_mpi_ranks([lo(invorder(0)), lo(invorder(1)), lo(invorder(2))], &
                                   this%mpi_ranks, this%flat_mpi_ranks, this%shape_mpi_ranks, irank, nproc)
    call diezdecomp_fill_is_bound(this, this%flat_mpi_ranks, this%shape_mpi_ranks, periods, ipencil, irank)
    this%order_halo = [0, 1, 2]
    this%size       = product(this%shape)
    do i=0,2
      this%offset6(i,:) = offset6_xyz(order(i),:)
    end do
  end subroutine

end module