module diezdecomp_api_generic
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
#else
  integer, parameter :: rp = dp
#endif
  real(dp), parameter :: LARGE = huge(1._dp)

  private

  ! --------------- public variables ---------------
  public :: diezdecomp_generic_fill_tr_obj         , &
            diezdecomp_generic_fill_hl_obj         , &
            diezdecomp_track_mpi_decomp            , &
            diezdecomp_transp_execute_generic_buf  , &
            diezdecomp_transp_execute_generic_nobuf, &
            diezdecomp_halos_execute_generic       , &
            diezdecomp_props_halo                  , &
            diezdecomp_props_transp                , &
            diezdecomp_parsed_mpi_ranks            , &
            diezdecomp_autotune_mode_trials        , &
            diezdecomp_n_warmup                    , &
            diezdecomp_summary_halo_autotuning     , &
            diezdecomp_summary_transp_autotuning

  ! --------------- types ---------------
  type diezdecomp_parsed_mpi_ranks
    integer              :: shape_mpi_ranks(0:2), irank, nproc
    integer, allocatable :: mpi_ranks(:,:,:), flat_mpi_ranks(:,:)
  end type

  contains

  subroutine diezdecomp_transp_execute_generic_buf(this, p_in, p_out, work, offset_3x2_in, offset_3x2_out, stream)
    implicit none
    type(diezdecomp_props_transp)      :: this
    real(rp), target                   :: p_in(:,:,:), p_out(:,:,:)
    real(rp), target                   :: work(0:)
    integer                            ::    offset6_in(0:2,0:1),    offset6_out(0:2,0:1)
    integer                 , optional :: offset_3x2_in(0:2,0:1), offset_3x2_out(0:2,0:1)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    stream_internal  =  diezdecomp_stream_default
    if (present(stream)) stream_internal = stream
    if (present(offset_3x2_in )) then ; offset6_in  = this%offset6_in  ; this%offset6_in  = offset_3x2_in  ; end if
    if (present(offset_3x2_out)) then ; offset6_out = this%offset6_out ; this%offset6_out = offset_3x2_out ; end if
    call diezdecomp_transp_execute(this, p_in, p_out, stream_internal, work)
    if (present(offset_3x2_in )) then ; this%offset6_in  = offset6_in  ; end if ! restore original offsets
    if (present(offset_3x2_out)) then ; this%offset6_out = offset6_out ; end if
    if (.not.present(stream)) then
      !$acc wait(stream_internal)
    end if
  end subroutine

  subroutine diezdecomp_transp_execute_generic_nobuf(this, p_in, p_out, offset_3x2_in, offset_3x2_out, stream)
    implicit none
    type(diezdecomp_props_transp)      :: this
    real(rp), target                   :: p_in(0:*), p_out(0:*)
    integer                            ::    offset6_in(0:2,0:1),    offset6_out(0:2,0:1)
    integer                 , optional :: offset_3x2_in(0:2,0:1), offset_3x2_out(0:2,0:1)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    stream_internal  =  diezdecomp_stream_default
    if (present(stream)) stream_internal = stream
    if (present(offset_3x2_in )) then ; offset6_in  = this%offset6_in  ; this%offset6_in  = offset_3x2_in  ; end if
    if (present(offset_3x2_out)) then ; offset6_out = this%offset6_out ; this%offset6_out = offset_3x2_out ; end if
    call diezdecomp_transp_execute(this, p_in, p_out, stream_internal)
    if (present(offset_3x2_in )) then ; this%offset6_in  = offset6_in  ; end if ! restore original offsets
    if (present(offset_3x2_out)) then ; this%offset6_out = offset6_out ; end if
    if (.not.present(stream)) then
      !$acc wait(stream_internal)
    end if
  end subroutine

  subroutine diezdecomp_halos_execute_generic(this, p, work, offset3_start, stream)
    implicit none
    type(diezdecomp_props_halo) :: this
    real(rp), contiguous :: p(0:,0:,0:), work(0:)
    integer                            :: offset3(0:2)
    integer                 , optional :: offset3_start(0:2)
    integer(acc_handle_kind), optional :: stream
    integer(acc_handle_kind)           :: stream_internal
    ! note:
    !   1) "offset3_start" are the offsets of p(0:,0:,0:) at the begging of each dimension
    !   2) technically, knowing the offsets at the end (offset3_end) is not necessary,
    !      because "p" is passed as a 3-D array (assumed shape is not used)
    stream_internal  =  diezdecomp_stream_default
    if (present(stream)) stream_internal = stream
    if (present(offset3_start)) then ; offset3  = this%offset3  ; this%offset3  = offset3_start  ; end if
    call diezdecomp_halos_execute(this, p, work, stream_internal)
    if (present(offset3_start)) then ; this%offset3 = offset3  ; end if ! restore original offsets
    if (.not.present(stream)) then
      !$acc wait(stream_internal)
    end if
  end subroutine

  subroutine diezdecomp_generic_fill_tr_obj(tr, obj_in, obj_out, sp_in_full, offset6_in, sp_out_full, offset6_out, &
                                            order_in, order_out, order_intermediate, &
                                            allow_alltoallv, ii, jj, wsize, allow_autotune_reorder, stream)
    implicit none
    type(diezdecomp_props_transp)      ::  tr
    type(diezdecomp_parsed_mpi_ranks)  ::  obj_in, obj_out
    integer(i8)                        ::  wsize
    integer                            ::  sp_in_full(0:2), sp_out_full(0:2), &
                                           order_in(0:2), order_out(0:2), order_intermediate(0:2), &
                                           lo_in(0:2), hi_in(0:2), offset6_in(0:2,0:1), offset6_out(0:2,0:1), &
                                           lo_out(0:2), hi_out(0:2), &
                                           ii, jj, kk, ii_jj_kk(0:2)
    logical                            ::  allow_alltoallv, allow_autotune_reorder
    integer(acc_handle_kind), optional ::  stream
    integer(acc_handle_kind)           ::  stream_internal
    stream_internal = diezdecomp_stream_default
    if (present(stream)) stream_internal = stream
    if (.not.(tr%initialized)) then
      block
        integer :: sp_in(0:2), sp_out(0:2)
        sp_in  = [sp_in_full(0)  - offset6_in(0,0)  - offset6_in(0,1)  , &
                  sp_in_full(1)  - offset6_in(1,0)  - offset6_in(1,1)  , &
                  sp_in_full(2)  - offset6_in(2,0)  - offset6_in(2,1)  ]
        sp_out = [sp_out_full(0) - offset6_out(0,0) - offset6_out(0,1) , &
                  sp_out_full(1) - offset6_out(1,0) - offset6_out(1,1) , &
                  sp_out_full(2) - offset6_out(2,0) - offset6_out(2,1) ]
        wsize = product(sp_in) + product(sp_out)
        if (wsize<1) wsize=1
        call get_lo_hi_3d(lo_in , hi_in , obj_in , sp_in , order_in )
        call get_lo_hi_3d(lo_out, hi_out, obj_out, sp_out, order_out)
      end block
      kk       = 3 - ii - jj
      ii_jj_kk = [ii,jj,kk]
      call diezdecomp_fill_transp_props(tr                                     , &
                                        order_intermediate                     , &
                                        order_in                               , &
                                        order_out                              , &
                                        lo_in, hi_in                           , &
                                        lo_out, hi_out                         , &
                                        ii_jj_kk                               , &
                                        obj_in%irank, obj_in%nproc, allow_alltoallv, &
                                        sp_in_full , offset6_in         , &
                                        sp_out_full, offset6_out        , &
                                        allow_autotune_reorder, stream_internal)
    end if
    ! now tr is ready to call:
    !   call diezdecomp_transp_execute_generic(tr, p_in, p_out, work, stream) ! size(work) --> wsize
  end subroutine

  subroutine diezdecomp_generic_fill_hl_obj(hl, obj_in, p_shape, offset6, ii, nh_xyz, order_halo, periodic, wsize, &
                                            force_halo_sync, autotuned_pack)
    implicit none
    integer(i8)                        ::  wsize
    type(diezdecomp_props_halo)        ::  hl
    type(diezdecomp_parsed_mpi_ranks)  ::  obj_in
    integer                            :: nh_xyz(0:2), ii, p_shape(0:2), order_halo(0:2), force_halo_sync, &
                                          autotuned_pack, offset6(0:2,0:1)
    logical                            :: periodic(0:2), is_per
    ! call check_bounds(ii  , 0, 2)
    if (.not.(hl%initialized)) then
      is_per  = periodic(ii)
      call diezdecomp_fill_halo_props(hl, ii, is_per, &
                                      obj_in%mpi_ranks, &
                                      obj_in%flat_mpi_ranks, &
                                      obj_in%shape_mpi_ranks, &
                                      nh_xyz, order_halo, &
                                      p_shape, offset6, obj_in%irank, obj_in%nproc, force_halo_sync, autotuned_pack)
    end if
    block
      integer :: jj,kk
      if      (ii == 0) then ; jj = 1; kk = 2
      else if (ii == 1) then ; jj = 0; kk = 2
      else if (ii == 2) then ; jj = 0; kk = 1 ; end if
      wsize = p_shape(diezdecomp_get_index_int2int(order_halo,jj))*p_shape(diezdecomp_get_index_int2int(order_halo,kk))*4*nh_xyz(ii)
      if (wsize<1) wsize=1
    end block
    ! now hl is ready to call:
    !   call diezdecomp_halos_execute(hl, p_in, work, stream) ! size(work) --> wsize
  end subroutine

  subroutine diezdecomp_track_mpi_decomp(lo_loc, obj_in, irank, nproc, order)
    type(diezdecomp_parsed_mpi_ranks)  ::  obj_in
    integer                            ::  lo_loc(0:2), lo(0:2), irank, nproc
    integer, optional                  ::  order(0:2)
    lo  =  lo_loc
    if (present(order)) lo  =  [lo_loc(diezdecomp_get_index_int2int(order,0)), &
                                lo_loc(diezdecomp_get_index_int2int(order,1)), &
                                lo_loc(diezdecomp_get_index_int2int(order,2))]
    obj_in%irank  =  irank
    obj_in%nproc  =  nproc
    call diezdecomp_fill_mpi_ranks(lo, obj_in%mpi_ranks, obj_in%flat_mpi_ranks, obj_in%shape_mpi_ranks, irank, nproc)
  end subroutine

  ! ---------------------------------------- diezdecomp utilities -------------------------------------------
  subroutine get_lo_hi_3d(lo, hi, obj, sz3_0, order)
    integer                          :: lo(0:2), hi(0:2), order(0:2), sz3(0:2), sz3_0(0:2), &
                                        ii, jj, kk, ipos, iaxis, mpi_ierr
    type(diezdecomp_parsed_mpi_ranks):: obj
    integer, allocatable             :: buf_flat(:), buf_3d(:,:,:)

    sz3 = [sz3_0(diezdecomp_get_index_int2int(order,0)), &
           sz3_0(diezdecomp_get_index_int2int(order,1)), &
           sz3_0(diezdecomp_get_index_int2int(order,2))]

    allocate(buf_3d(0:obj%shape_mpi_ranks(0)-1, &
                    0:obj%shape_mpi_ranks(1)-1, &
                    0:obj%shape_mpi_ranks(2)-1), &
             buf_flat(0:(3*obj%nproc-1)))
    call MPI_Allgather(sz3, 3, mpi_int, buf_flat, 3, mpi_int, mpi_comm_world, mpi_ierr)
    do   iaxis = 0, 2
      block
        integer :: i,j,k
        do ipos  = 0, obj%nproc-1
          i = obj%flat_mpi_ranks(ipos,0)
          j = obj%flat_mpi_ranks(ipos,1)
          k = obj%flat_mpi_ranks(ipos,2)
          buf_3d(i,j,k) = buf_flat(3*ipos+iaxis)
        end do
      end block
      ii = obj%flat_mpi_ranks(obj%irank,0)
      jj = obj%flat_mpi_ranks(obj%irank,1)
      kk = obj%flat_mpi_ranks(obj%irank,2)
      lo(iaxis) = 0
      hi(iaxis) = 0
      if      (iaxis == 0) then
        block
          integer :: i
          do i=0,ii
            lo(iaxis) = hi(iaxis)
            hi(iaxis) = hi(iaxis) + buf_3d(i,jj,kk)
          end do
        end block
      else if (iaxis == 1) then
        block
          integer :: j
          do j=0,jj
            lo(iaxis) = hi(iaxis)
            hi(iaxis) = hi(iaxis) + buf_3d(ii,j,kk)
          end do
        end block
      else if (iaxis == 2) then
        block
          integer :: k
          do k=0,kk
            lo(iaxis) = hi(iaxis)
            hi(iaxis) = hi(iaxis) + buf_3d(ii,jj,k)
          end do
        end block
      end if
    end do
    hi = hi - 1
    deallocate(buf_3d, buf_flat)
    block
      integer :: A(0:2)
      A(0:2) = lo(0:2) ; lo = [A(order(0)), A(order(1)), A(order(2))]
      A(0:2) = hi(0:2) ; hi = [A(order(0)), A(order(1)), A(order(2))]
    end block
  end subroutine

end module
