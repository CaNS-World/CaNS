module diezdecomp_core
  use mpi
! #if defined(_DIEZDECOMP_NCCL)
!   use nccl
! #endif
  use, intrinsic :: iso_fortran_env, only: i8 => int64, sp => real32, dp => real64
  use, intrinsic :: iso_c_binding  , only: c_loc,c_size_t
#if defined(_OPENACC)
  use openacc, only: acc_handle_kind
#else
  use, intrinsic :: iso_fortran_env, only: acc_handle_kind => int64
#endif
  implicit none
#if defined(_DIEZDECOMP_SINGLE)
  integer, parameter :: rp = sp
  integer, parameter :: MPI_REAL_RP = MPI_REAL
#else
  integer, parameter :: rp = dp
  integer, parameter :: MPI_REAL_RP = MPI_DOUBLE_PRECISION
#endif
  real(dp), parameter            :: LARGE = huge(1._dp)
  integer(acc_handle_kind), save :: diezdecomp_stream_default = 1
  logical, save                  :: diezdecomp_mode_bench_avg = .false.
  private

  ! --------------- public variables ---------------
  public :: diezdecomp_props_halo              , &
            diezdecomp_props_transp            , &
            diezdecomp_fill_halo_props         , &
            diezdecomp_fill_transp_props       , &
            diezdecomp_halos_execute           , &
            diezdecomp_transp_execute          , &
            diezdecomp_fill_mpi_ranks          , &
            diezdecomp_get_index_int2int       , &
            diezdecomp_get_rank_id             , &
            diezdecomp_stream_default          , &
            diezdecomp_autotune_mode_trials    , &
            diezdecomp_n_warmup                , &
            diezdecomp_summary_halo_autotuning , &
            diezdecomp_summary_transp_autotuning

! #if defined(_DIEZDECOMP_NCCL)
!   type(ncclUniqueId)       , save :: nccl_id
!   type(ncclResult)         , save :: nccl_stat
!   type(ncclComm)           , save :: nccl_comm
!   logical                  , save :: nccl_is_init = .false.
!   integer(cuda_stream_kind), save :: nccl_stream
! #endif

  integer, save :: diezdecomp_autotune_mode_trials = 20
  integer, save :: diezdecomp_n_warmup             = 3

  ! --------------- types ---------------

  type diezdecomp_props_halo
    logical :: args_active(0:1)
    integer :: irank, loc_shape(0:2), nh_ijk(0:2), side_numel, nh_dir, dir_pick  , &
               args_other_id(0:1), args_inds_recv_i0(0:1), args_inds_recv_i1(0:1), &
                                   args_inds_send_i0(0:1), args_inds_send_i1(0:1), &
               args_tag_send(0:1), args_tag_recv(0:1), &
               args_pos_send(0:1), args_pos_recv(0:1), args_rel_pos_align, &
               all_request(0:3), all_status(MPI_STATUS_SIZE,0:3), autotuned_pack, &
               offset3(0:2), halo_mpi_mode, nproc
    logical :: initialized = .false.
    real(rp) :: elapsed_sendrecv, elapsed_isendirecv, wtime_2, wtime_3
  end type

  type diezdecomp_props_transp
    integer              :: irank, nproc, mpi_comm_kk, irank_comm_kk, &
                            size_localSend, size_localRecv, pos_start_localSend, pos_start_localRecv, &
                            send_mode_op_batched, send_mode_op_simul, send_autotuned , send_i2b_Mij(0:2,0:2), &
                            recv_mode_op_batched, recv_mode_op_simul, recv_autotuned , recv_i2b_Mij(0:2,0:2), &
                            best_recv_mode_op_simul(0:5), best_recv_autotuned(0:5), best_recv_mode_op_batched(0:5), &
                            best_send_mode_op_simul(0:5), best_send_autotuned(0:5), best_send_mode_op_batched(0:5), &
                            sp_in_0, sp_in_1, sp_in_2, sp_out_0, sp_out_1, sp_out_2, order_a(0:2), order_b(0:2),&
                            offset6_in(0:2,0:1), offset6_out(0:2,0:1), &
                            numel_p_in, numel_p_out, chosen_reorder(0:2), reorder_ijk(0:2,0:5)
    integer, allocatable :: send_ids(:), send_sizes(:), send_tags(:), send_pos_start(:), &
                            recv_ids(:), recv_sizes(:), recv_tags(:), recv_pos_start(:), &
                            send_i0_st(:), send_i1_st(:), send_i2_st(:), send_lshape(:,:), &
                            recv_i0_st(:), recv_i1_st(:), recv_i2_st(:), recv_lshape(:,:), &
                            send_b_i0(:), send_b_i1(:), send_b_i2(:), &
                            recv_b_i0(:), recv_b_i1(:), recv_b_i2(:), &
                            send_all_n0(:,:,:), send_all_n1(:,:,:), send_all_n2(:,:,:), send_all_pos_start(:,:,:),&
                            recv_all_n0(:,:,:), recv_all_n1(:,:,:), recv_all_n2(:,:,:), recv_all_pos_start(:,:,:),&
                            send_i2b_s0(:), send_i2b_i0a(:), &
                            send_i2b_s1(:), send_i2b_i1a(:), &
                            send_i2b_s2(:), send_i2b_i2a(:), &
                            recv_i2b_s0(:), recv_i2b_i0a(:), &
                            recv_i2b_s1(:), recv_i2b_i1a(:), &
                            recv_i2b_s2(:), recv_i2b_i2a(:), &
                            all_request(:), all_status(:,:), all_rank_comm_kk(:)
    logical              :: initialized = .false.
    logical              :: use_alltoallv, use_isendirecv, use_nccl, chosen_backend, single_device, &
                            apply_autotune_reorder
    real(rp)             :: elapsed_time_reorder(0:5), elapsed_alltoallv,elapsed_isendirecv,elapsed_nccl,&
                            recv_elapsed_simul_top, recv_elapsed_batched_top, recv_elapsed_simul_modes(0:5) , &
                            send_elapsed_simul_top, send_elapsed_batched_top, send_elapsed_simul_modes(0:5) , &
                            recv_elapsed_batched_modes(0:5) , &
                            send_elapsed_batched_modes(0:5)
  end type

  contains

  ! ---------------------------------------- diezdecomp halo ------------------------------------------------
  subroutine diezdecomp_fill_halo_props(this, ii, is_per, ranks_ii, flat_ranks_ii, shape_ranks_ii, nh_xyz, abs_order, &
                                 loc_shape_full, offset6, irank, nproc, halo_mpi_mode, autotuned_pack)
    implicit none
    type(diezdecomp_props_halo) :: this
    integer                    :: ii, nh_xyz(0:2), abs_order(0:2), irank, jj, kk, rel_pos, L, p0, p1, &
                                  k, size_dir, loc_shape(0:2), loc_shape_full(0:2), loc_shape_padded(0:2), nproc, &
                                  ranks_ii(0:,0:,0:), flat_ranks_ii(0:,0:), shape_ranks_ii(0:2), halo_mpi_mode, &
                                  autotuned_pack, offset6(0:2,0:1)
    logical                    :: is_per
    ! call check_bounds(ii, 0, 2)
    loc_shape_padded = [loc_shape_full(0) - offset6(0,0) - offset6(0,1),&
                        loc_shape_full(1) - offset6(1,0) - offset6(1,1),&
                        loc_shape_full(2) - offset6(2,0) - offset6(2,1)]
    this%irank = irank
    this%nproc = nproc

    this%elapsed_sendrecv   = -1
    this%elapsed_isendirecv = -1
    this%wtime_2            = -1
    this%wtime_3            = -1

    do k=0,2
      this%nh_ijk(k) = nh_xyz(abs_order(k))
      loc_shape(k)   = loc_shape_padded(k) - 2*this%nh_ijk(k)
    end do

    this%loc_shape(:) = loc_shape(:)

    if      (ii == 0) then ; jj = 1; kk = 2
    else if (ii == 1) then ; jj = 0; kk = 2
    else if (ii == 2) then ; jj = 0; kk = 1
    else                   ; error stop 'ERROR: unrecognized ii' ! , ii
    end if

    this%nh_dir             = nh_xyz(ii)
    this%side_numel         = loc_shape_padded(diezdecomp_get_index_int2int(abs_order,jj))* &
                              loc_shape_padded(diezdecomp_get_index_int2int(abs_order,kk))* this%nh_dir
    rel_pos                 = get_position_in_3darr(flat_ranks_ii, irank, ii)
    this%dir_pick           = diezdecomp_get_index_int2int(abs_order,ii)

    size_dir                = shape_ranks_ii(ii)
    this%args_rel_pos_align = rel_pos
    if (mod(rel_pos,2) == 0) then
      p0 = 0  ;  p1 = 1
    else
      p0 = 1  ;  p1 = 0 ! this is necessary for force_sync_halo, and good for isend/irecv
    end if

    this%args_active      (p0) = is_per.or.(rel_pos>0)
    this%args_other_id    (p0) = diezdecomp_get_rank_id(ranks_ii, flat_ranks_ii, shape_ranks_ii, irank, ii, -1)
    this%args_inds_recv_i0(p0) = 0
    this%args_inds_recv_i1(p0) = this%nh_dir - 1
    this%args_inds_send_i0(p0) = this%nh_dir
    this%args_inds_send_i1(p0) = 2*this%nh_dir - 1
    this%args_tag_send    (p0) = 2
    this%args_tag_recv    (p0) = 1
    this%args_pos_send    (p0) = 0
    this%args_pos_recv    (p0) = 2*this%side_numel
    L                          = this%loc_shape(this%dir_pick)

    this%args_active      (p1) = is_per.or.(rel_pos<(size_dir-1))
    this%args_other_id    (p1) = diezdecomp_get_rank_id(ranks_ii, flat_ranks_ii, shape_ranks_ii, irank, ii,  1)
    this%args_inds_recv_i0(p1) = L + this%args_inds_send_i0(p0)
    this%args_inds_recv_i1(p1) = L + this%args_inds_send_i1(p0)
    this%args_inds_send_i0(p1) = L + this%args_inds_recv_i0(p0)
    this%args_inds_send_i1(p1) = L + this%args_inds_recv_i1(p0)
    this%args_tag_send    (p1) = 1
    this%args_tag_recv    (p1) = 2
    this%args_pos_send    (p1) =   this%side_numel
    this%args_pos_recv    (p1) = 3*this%side_numel

    this%halo_mpi_mode  = halo_mpi_mode
    this%autotuned_pack = autotuned_pack
    this%offset3(0:2)   = offset6(0:2,0)
    this%initialized    = .true.
  end subroutine

  subroutine diezdecomp_halos_execute(this, p, buffer, stream)
    implicit none
    type(diezdecomp_props_halo) :: this
    real(rp) :: p(0:,0:,0:), buffer(0:*)
    real(dp) :: wtime_2, wtime_3, temp_real
    integer  :: i,mpi_ierr
    integer(acc_handle_kind) :: stream
    if (this%nh_dir>0) then
      wtime_2 = 0
      wtime_3 = 0
      do i = 0,1
        if (this%args_active(i)) then
          call pack_slices_halo(p, buffer, this%args_pos_send(i), this%args_inds_send_i0(i), &
                                                                  this%args_inds_send_i1(i), &
                                                                  this%dir_pick, this%loc_shape, this%nh_ijk, .true., &
                                                                  this%offset3, this%autotuned_pack,wtime_2,wtime_3, stream)
        end if
      end do
      call comm_halo_mpi(this, buffer, stream)
      do i = 0,1
        if (this%args_active(i)) then
          call pack_slices_halo(p, buffer, this%args_pos_recv(i), this%args_inds_recv_i0(i), &
                                                                  this%args_inds_recv_i1(i), &
                                                                  this%dir_pick, this%loc_shape, this%nh_ijk, .false.,&
                                                                  this%offset3, this%autotuned_pack,wtime_2,wtime_3, stream)

        end if
      end do
      if ((this%autotuned_pack < 2).or.(3 < this%autotuned_pack)) then
        temp_real = wtime_2/(this%nproc*1._dp)
        call MPI_Allreduce(temp_real, wtime_2, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        temp_real = wtime_3/(this%nproc*1._dp)
        call MPI_Allreduce(temp_real, wtime_3, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        this%autotuned_pack = 3
        if (wtime_2<wtime_3) this%autotuned_pack = 2
        this%wtime_2  =  wtime_2
        this%wtime_3  =  wtime_3
      end if
    end if
  end subroutine

  recursive subroutine comm_halo_mpi(this, buffer, stream)
    implicit none
    type(diezdecomp_props_halo) :: this
    real(rp) :: buffer(0:*)
    integer :: i, mpi_ierr, ireq, mode_best, mode_try,pos
    integer(acc_handle_kind) :: stream
    real(dp) :: elapsed_time,elapsed_best,temp_wtime, temp_real

    if ((this%halo_mpi_mode<1).or.(2<this%halo_mpi_mode)) then
      elapsed_best = LARGE
      do mode_try=1,2
        this%halo_mpi_mode = mode_try
        do pos=1,diezdecomp_n_warmup ; call comm_halo_mpi(this, buffer, stream) ; end do
        if (diezdecomp_mode_bench_avg) then ; elapsed_time = 0.
        else                                ; elapsed_time = LARGE ; end if
        do pos=1,diezdecomp_autotune_mode_trials
          !$acc wait(stream)
          temp_wtime = MPI_Wtime()
          call comm_halo_mpi(this, buffer, stream)
          !$acc wait(stream)
          temp_wtime = MPI_Wtime() - temp_wtime
          if (diezdecomp_mode_bench_avg) then
            elapsed_time = elapsed_time + temp_wtime/(diezdecomp_autotune_mode_trials*1._dp)
          else
            elapsed_time = min(elapsed_time, temp_wtime)
          end if
        end do
        temp_real = elapsed_time/(this%nproc*1._dp)
        call MPI_Allreduce(temp_real, elapsed_time, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        if (mode_try==1) this%elapsed_sendrecv   = elapsed_time
        if (mode_try==2) this%elapsed_isendirecv = elapsed_time
        if (elapsed_time<elapsed_best) then
          elapsed_best = elapsed_time
          mode_best    = mode_try
        end if
      end do
      this%halo_mpi_mode = mode_best
    else
    ireq = 0
    !$acc wait(stream)
    do i = 0,1
      if (this%args_active(i)) then
        if (this%irank /= this%args_other_id(i)) then
          if (this%halo_mpi_mode==1) then
            !$acc host_data use_device(buffer)
call MPI_Sendrecv(buffer(this%args_pos_send(i)), this%side_numel, MPI_REAL_RP, this%args_other_id(i), this%args_tag_send(i),&
            buffer(this%args_pos_recv(i)), this%side_numel, MPI_REAL_RP, this%args_other_id(i), this%args_tag_recv(i),&
            mpi_comm_world, this%all_status(:,0), mpi_ierr)
            !$acc end host_data
          else if (this%halo_mpi_mode==2) then
call isendrecv_1d_buffer(.false., buffer, this%args_pos_recv(i), this%side_numel, &
                         this%args_other_id(i),this%args_tag_recv(i), mpi_comm_world, this%all_request(ireq  ))
call isendrecv_1d_buffer(.true., buffer, this%args_pos_send(i), this%side_numel, &
                         this%args_other_id(i),this%args_tag_send(i), mpi_comm_world, this%all_request(ireq+1))
            ireq = ireq + 2
          else
            error stop 'ERROR: this%halo_mpi_mode out-of-range'
          end if
        end if
      end if
    end do
      do i = 0,1
      if (this%args_active(i)) then
        if (this%irank == this%args_other_id(i)) then
if (this%irank /= this%args_other_id(1-i)) error stop 'ERROR: halos ids:' ! , this%irank,this%args_other_id
          call copy_1d_buffer(buffer, buffer, this%args_pos_send(i), this%args_pos_recv(1-i), this%side_numel, stream)
        end if
      end if
    end do
    if (ireq>0) call mpi_waitall(ireq, this%all_request, this%all_status, mpi_ierr)
    !$acc wait(stream)
    end if
  end subroutine

  subroutine pack_slices_halo(p, buffer, pos_ini, inds_pick_i0, inds_pick_i1, dir_pick, loc_shape, nh_ijk, mode_fwd, &
                              offset3, autotuned, wtime_2, wtime_3, stream)
    real(rp) :: p(0:,0:,0:), buffer(0:*)
    integer  :: loc_shape(0:), &
                pos_ini, inds_pick_i0, inds_pick_i1, dir_pick, nh_ijk(0:2), &
                autotuned,iter1, mode_iter, &
                offset3(0:2)
    logical  :: mode_fwd
    integer(acc_handle_kind) :: stream
    real(dp) :: elapsed_time, temp_wtime, wtime_2, wtime_3
    if ((autotuned < 2).or.(3 < autotuned)) then
      do mode_iter=2,3
        if (diezdecomp_mode_bench_avg) then ; elapsed_time = 0.
        else                                ; elapsed_time = LARGE ; end if
        do iter1 = 1,diezdecomp_autotune_mode_trials
          if      (mode_iter==2) then
            !$acc wait(stream)
            temp_wtime = MPI_Wtime()
            call pack_slices_2d(p, buffer, pos_ini, inds_pick_i0, inds_pick_i1, dir_pick, loc_shape, nh_ijk, mode_fwd, &
                                offset3, stream)
            !$acc wait(stream)
            temp_wtime = MPI_Wtime() - temp_wtime
          else if (mode_iter==3) then
            !$acc wait(stream)
            temp_wtime = MPI_Wtime()
            call pack_slices_3d(p, buffer, pos_ini, inds_pick_i0, inds_pick_i1, dir_pick, loc_shape, nh_ijk, mode_fwd, &
                                offset3, stream)
            !$acc wait(stream)
            temp_wtime = MPI_Wtime() - temp_wtime
          else
            error stop 'ERROR: Unrecognized mode_iter'
          end if
          if (diezdecomp_mode_bench_avg) then
            elapsed_time = elapsed_time + temp_wtime/(diezdecomp_autotune_mode_trials*1._dp)
          else
            elapsed_time = min(elapsed_time, temp_wtime)
          end if
        end do
        if (mode_iter==2) wtime_2 = wtime_2 + elapsed_time
        if (mode_iter==3) wtime_3 = wtime_3 + elapsed_time
      end do
    else
      !$acc wait(stream)
      if      (autotuned == 2) then
        call pack_slices_2d(p, buffer, pos_ini, inds_pick_i0, inds_pick_i1, dir_pick, loc_shape, nh_ijk, mode_fwd, &
                            offset3, stream)
      else if (autotuned == 3) then
        call pack_slices_3d(p, buffer, pos_ini, inds_pick_i0, inds_pick_i1, dir_pick, loc_shape, nh_ijk, mode_fwd, &
                            offset3, stream)
      else
        error stop 'ERROR: autotuned out-of-range (pack_slices_halo)'
      endif
      !$acc wait(stream)
    end if
  end subroutine

  subroutine pack_slices_2d(p, buffer, pos_ini, inds_pick_i0, inds_pick_i1, dir_pick, loc_shape, nh_ijk, mode_fwd, &
                            offset3, stream)
    implicit none
    real(rp) :: p(0:,0:,0:), buffer(0:*)
    integer  :: loc_shape(0:), &
                pos_ini, inds_pick_i0, inds_pick_i1, dir_pick, ni, nj, nk, nh_ijk(0:2), &
                i,j,k,pos,pos_start, offset3(0:2), di,dj,dk
    logical  :: mode_fwd
    integer(acc_handle_kind) :: stream
    ! call check_bounds(dir_pick,0, 2)
    di = offset3(0) ; dj = offset3(1) ; dk = offset3(2)
    ni = loc_shape(0)+ 2*nh_ijk(0)
    nj = loc_shape(1)+ 2*nh_ijk(1)
    nk = loc_shape(2)+ 2*nh_ijk(2)
    pos = pos_ini + 0
    if (dir_pick == 0) then
      do i=inds_pick_i0,inds_pick_i1
        pos_start = pos + 0
        if (mode_fwd) then
          !$acc parallel loop gang vector collapse(2) default(present) async(stream)
          do k = 0, nk-1
          do j = 0, nj-1
            buffer(pos_start + k*nj + j) = p(i+di,j+dj,k+dk)
          end do
          end do
        else
          !$acc parallel loop gang vector collapse(2) default(present) async(stream)
          do k = 0, nk-1
          do j = 0, nj-1
            p(i+di,j+dj,k+dk) = buffer(pos_start + k*nj + j)
          end do
          end do
        end if
        pos = pos + nk*nj
      end do
    else if (dir_pick == 1) then
      do j=inds_pick_i0,inds_pick_i1
        pos_start = pos + 0
        if (mode_fwd) then
          !$acc parallel loop gang vector collapse(2) default(present) async(stream)
          do k = 0, nk-1
          do i = 0, ni-1
            buffer(pos_start + k*ni + i) = p(i+di,j+dj,k+dk)
          end do
          end do
        else
          !$acc parallel loop gang vector collapse(2) default(present) async(stream)
          do k = 0, nk-1
          do i = 0, ni-1
            p(i+di,j+dj,k+dk) = buffer(pos_start + k*ni + i)
          end do
          end do
        end if
        pos = pos + nk*ni
      end do
    else if (dir_pick == 2) then
      do k=inds_pick_i0,inds_pick_i1
        pos_start = pos + 0
        if (mode_fwd) then
          !$acc parallel loop gang vector collapse(2) default(present) async(stream)
          do j = 0, nj-1
          do i = 0, ni-1
            buffer(pos_start + j*ni + i) = p(i+di,j+dj,k+dk)
          end do
          end do
        else
          !$acc parallel loop gang vector collapse(2) default(present) async(stream)
          do j = 0, nj-1
          do i = 0, ni-1
            p(i+di,j+dj,k+dk) = buffer(pos_start + j*ni + i)
          end do
          end do
        end if
        pos = pos + nj*ni
      end do
    else
      error stop 'ERROR: dir_pick not detected: ' ! , dir_pick
    end if
  end subroutine

  subroutine pack_slices_3d(p, buffer, pos_ini, inds_pick_i0, inds_pick_i1, dir_pick, loc_shape, nh_ijk, mode_fwd, &
                            offset3, stream)
    implicit none
    real(rp) :: p(0:,0:,0:), buffer(0:*)
    integer  :: loc_shape(0:), &
                pos_ini, inds_pick_i0, inds_pick_i1, dir_pick, ni, nj, nk, nh_ijk(0:2), &
                i,j,k,pos_start,i0,i1,j0,j1,k0,k1,nji,di,dj,dk, offset3(0:2)
    logical  :: mode_fwd
    integer(acc_handle_kind) :: stream
    ! call check_between(dir_pick,0, 2)
    ni = loc_shape(0)+ 2*nh_ijk(0)
    nj = loc_shape(1)+ 2*nh_ijk(1)
    nk = loc_shape(2)+ 2*nh_ijk(2)
    i0 = 0 ; i1 = ni-1
    j0 = 0 ; j1 = nj-1
    k0 = 0 ; k1 = nk-1
    di = offset3(0) ; dj = offset3(1) ; dk = offset3(2)
    if      (dir_pick ==0) then  ; i0 = inds_pick_i0 ; i1 = inds_pick_i1
    else if (dir_pick ==1) then  ; j0 = inds_pick_i0 ; j1 = inds_pick_i1
    else if (dir_pick ==2) then  ; k0 = inds_pick_i0 ; k1 = inds_pick_i1
    else                         ; error stop 'ERROR: dir_pick not detected: ' ! , dir_pick
    end if
    ni  = i1 - i0 + 1
    nj  = j1 - j0 + 1
    nk  = k1 - k0 + 1
    nji = nj*ni
    if      (dir_pick ==0) then ; pos_start = pos_ini - i0
    else if (dir_pick ==1) then ; pos_start = pos_ini - j0*ni
    else if (dir_pick ==2) then ; pos_start = pos_ini - k0*nji
    end if
    if (mode_fwd) then
      !$acc parallel loop gang vector collapse(3) default(present) async(stream)
      do k = k0,k1
      do j = j0,j1
      do i = i0,i1
        buffer(pos_start + k*nji + j*ni + i) = p(i+di,j+dj,k+dk)
      end do
      end do
      end do
    else
      !$acc parallel loop gang vector collapse(3) default(present) async(stream)
      do k = k0,k1
      do j = j0,j1
      do i = i0,i1
        p(i+di,j+dj,k+dk) = buffer(pos_start + k*nji + j*ni + i)
      end do
      end do
      end do
    end if
  end subroutine

  ! ---------------------------------------- diezdecomp transpose -------------------------------------------
  subroutine diezdecomp_fill_transp_props(this, abs_reorder, loc_order_a, loc_order_b, lo_a, hi_a, lo_b, hi_b, &
                                  ii_jj_kk, irank, nproc, allow_alltoallv, sp_in_full, offset6_in, &
                                  sp_out_full, offset6_out, allow_autotune_reorder, stream)
    implicit none
    type(diezdecomp_props_transp) :: this
    integer                  :: ks, kr, kreq, pos, mpi_ierr, &
                                nproc, irank, prev_pos, jj_pos, kk_pos, &
                                abs_reorder(0:2), loc_order_a(0:2), loc_order_b(0:2), &
                                lo_a(0:2),     hi_a(0:2),     lo_b(0:2),     hi_b(0:2), &
                                sp_in_0, sp_in_1, sp_in_2, sp_out_0, sp_out_1, sp_out_2, ii_jj_kk(0:2), &
                                offset6_in(0:2,0:1), offset6_out(0:2,0:1), sp_in_full(0:2), sp_out_full(0:2)
    integer, allocatable     :: all_lo_a(:,:), all_hi_a(:,:), all_lo_b(:,:), all_hi_b(:,:), &
                                temp_buffer(:,:), temp_buffer_sort(:,:), aux1(:,:), aux2(:,:)
    logical                  :: allow_alltoallv, allow_autotune_reorder
    integer(acc_handle_kind) :: stream
    this%irank = irank
    this%nproc = nproc
    this%order_a(0:2) = loc_order_a(0:2)
    this%order_b(0:2) = loc_order_b(0:2)

    sp_in_0  = sp_in_full(0)  - offset6_in(0,0)  - offset6_in(0,1)
    sp_in_1  = sp_in_full(1)  - offset6_in(1,0)  - offset6_in(1,1)
    sp_in_2  = sp_in_full(2)  - offset6_in(2,0)  - offset6_in(2,1)
    sp_out_0 = sp_out_full(0) - offset6_out(0,0) - offset6_out(0,1)
    sp_out_1 = sp_out_full(1) - offset6_out(1,0) - offset6_out(1,1)
    sp_out_2 = sp_out_full(2) - offset6_out(2,0) - offset6_out(2,1)

    this%offset6_in  = offset6_in
    this%offset6_out = offset6_out

    this%sp_in_0 = sp_in_0
    this%sp_in_1 = sp_in_1
    this%sp_in_2 = sp_in_2

    this%sp_out_0 = sp_out_0
    this%sp_out_1 = sp_out_1
    this%sp_out_2 = sp_out_2

    this%numel_p_in  = sp_in_0  * sp_in_1  * sp_in_2
    this%numel_p_out = sp_out_0 * sp_out_1 * sp_out_2

    block
      integer :: i,j,k,pos
      this%elapsed_time_reorder               = -1
      this%elapsed_alltoallv                  = -1
      this%elapsed_isendirecv                 = -1
      this%elapsed_nccl                       = -1
      this%reorder_ijk                        = -1
      this%chosen_reorder                     = -1
      this%send_elapsed_simul_top             = -1
      this%send_elapsed_batched_top           = -1
      this%send_elapsed_simul_modes           = -1
      this%send_elapsed_batched_modes         = -1
      this%recv_elapsed_simul_top             = -1
      this%recv_elapsed_batched_top           = -1
      this%recv_elapsed_simul_modes           = -1
      this%recv_elapsed_batched_modes         = -1
      pos = -1
      do i=0,2
        do j=0,2
          if (i /= j) then
            pos = pos + 1
            k   = 3-i-j
            this%reorder_ijk(0:2,pos) = [i,j,k]
          end if
        end do
      end do
      if (pos /= 5) error stop 'ERROR: autotune reorder: pos /= 5'
    end block

    ! i) intersect ijk coordinates to find subdomains
    allocate(all_lo_a        (0:2, 0:nproc-1), &
             all_hi_a        (0:2, 0:nproc-1), &
             all_lo_b        (0:2, 0:nproc-1), &
             all_hi_b        (0:2, 0:nproc-1), &
             temp_buffer     (0:nproc-1,0:7) , &
             temp_buffer_sort(0:nproc-1,0:7) )

    call MPI_Allgather(lo_a, 3, mpi_int, all_lo_a, 3, mpi_int, mpi_comm_world, mpi_ierr)
    call MPI_Allgather(hi_a, 3, mpi_int, all_hi_a, 3, mpi_int, mpi_comm_world, mpi_ierr)
    call MPI_Allgather(lo_b, 3, mpi_int, all_lo_b, 3, mpi_int, mpi_comm_world, mpi_ierr)
    call MPI_Allgather(hi_b, 3, mpi_int, all_hi_b, 3, mpi_int, mpi_comm_world, mpi_ierr)

    ! send information
    call transp_match_intervals(lo_a, hi_a, all_lo_b, all_hi_b, &
                             temp_buffer, temp_buffer_sort, ks, nproc, loc_order_a, loc_order_b)
    allocate(this%send_i0_st    (0:ks    ), &
             this%send_i1_st    (0:ks    ), &
             this%send_i2_st    (0:ks    ), &
             this%send_b_i0     (0:ks    ), &
             this%send_b_i1     (0:ks    ), &
             this%send_b_i2     (0:ks    ), &
             this%send_ids      (0:ks    ), &
             this%send_lshape   (0:ks,0:2), &
             this%send_sizes    (0:ks    ), &
             this%send_pos_start(0:ks    ), &
             this%send_tags     (0:ks    ))
    this%send_i0_st (0:ks    ) = temp_buffer(0:ks,  0) - lo_a(0)
    this%send_i1_st (0:ks    ) = temp_buffer(0:ks,  1) - lo_a(1)
    this%send_i2_st (0:ks    ) = temp_buffer(0:ks,  2) - lo_a(2)
    this%send_ids   (0:ks    ) = temp_buffer(0:ks,  3)
    this%send_lshape(0:ks,0:2) = temp_buffer(0:ks,4:6)

    ! recv information
    call transp_match_intervals(lo_b, hi_b, all_lo_a, all_hi_a, &
                                temp_buffer, temp_buffer_sort, kr, nproc, loc_order_b, loc_order_a)
    allocate(this%recv_i0_st    (0:kr    ), &
             this%recv_i1_st    (0:kr    ), &
             this%recv_i2_st    (0:kr    ), &
             this%recv_b_i0     (0:kr    ), &
             this%recv_b_i1     (0:kr    ), &
             this%recv_b_i2     (0:kr    ), &
             this%recv_ids      (0:kr    ), &
             this%recv_lshape   (0:kr,0:2), &
             this%recv_sizes    (0:kr    ), &
             this%recv_pos_start(0:kr    ), &
             this%recv_tags     (0:kr    ))
    this%recv_i0_st (0:kr    ) = temp_buffer(0:kr,  0) - lo_b(0)
    this%recv_i1_st (0:kr    ) = temp_buffer(0:kr,  1) - lo_b(1)
    this%recv_i2_st (0:kr    ) = temp_buffer(0:kr,  2) - lo_b(2)
    this%recv_ids   (0:kr    ) = temp_buffer(0:kr,  3)
    this%recv_lshape(0:kr,0:2) = temp_buffer(0:kr,4:6)

    ! ii) complete recv/send variables
    kreq = ks + kr + 1
    allocate(this%all_request(0:kreq), this%all_status(MPI_STATUS_SIZE, 0:kreq), this%all_rank_comm_kk(0:nproc-1))

    this%best_send_mode_op_batched = -1
    this%best_recv_mode_op_batched = -1
    this%best_send_mode_op_simul   = -1
    this%best_recv_mode_op_simul   = -1
    this%best_send_autotuned       = -1
    this%best_recv_autotuned       = -1

    prev_pos = 0
    do pos=0,ks
      this%send_sizes(pos)     = product(this%send_lshape(pos,0:2))
      this%send_tags(pos)      = 1 + this%irank * nproc + this%send_ids(pos) ! 1 + sender * nproc + recv
      this%send_pos_start(pos) = prev_pos
      prev_pos                 = prev_pos + this%send_sizes(pos)
    end do
    prev_pos = 0
    do pos=0,kr
      this%recv_sizes(pos)     = product(this%recv_lshape(pos,0:2))
      this%recv_tags(pos)      = 1 + this%recv_ids(pos) * nproc + this%irank ! 1 + sender * nproc + recv
      this%recv_pos_start(pos) = prev_pos
      prev_pos                 = prev_pos + this%recv_sizes(pos)
    end do

    ! iii) single_device calls (if applicable)
    this%single_device = ((size(this%send_ids,1) == 1).and.&
                          (size(this%recv_ids,1) == 1).and.&
                          (this%send_ids(0) == this%irank).and.&
                          (this%recv_ids(0) == this%irank))

    !this%recv_tags(:) = 0
    !this%send_tags(:) = 0

    this%use_isendirecv = .true. ! always possible
    this%chosen_backend = .false.
!#if defined(_DIEZDECOMP_NCCL)
!    this%use_nccl       = .true.
!    if (.not.nccl_is_init) then
!      nccl_is_init = .true.
!      if (this%irank == 0) nccl_stat = ncclGetUniqueId(nccl_id)
!      call MPI_Bcast(nccl_id%internal, int(sizeof(nccl_id%internal)), MPI_BYTE, 0, mpi_comm_world, mpi_ierr)
!      nccl_stat = ncclCommInitRank(nccl_comm, nproc, nccl_id, this%irank)
!      mpi_ierr = cudaStreamCreateWithFlags(nccl_stream, cudaStreamNonBlocking)
!    end if
!#else
    this%use_nccl       = .false.
!#endif

    ! iv) setup alltoallv communication (if applicable)
    this%use_alltoallv = allow_alltoallv.and.(.not.this%single_device)
    if (ks /= kr) then
      this%use_alltoallv = .false.
    else
      do pos=0,ks
        if (this%send_ids(pos) /= this%recv_ids(pos)) this%use_alltoallv = .false.
      end do
    endif
    block
      integer :: vmin,vmax
      call MPI_Allreduce(ks, vmax, 1, mpi_int, mpi_max, mpi_comm_world, mpi_ierr)
      call MPI_Allreduce(ks, vmin, 1, mpi_int, mpi_min, mpi_comm_world, mpi_ierr)
      if (vmin /= vmax) this%use_alltoallv = .false.
    end block
    if (this%use_alltoallv) then
      block
        integer :: k0,k1,s(0:2), temp_int, ii,jj,kk
        ii = ii_jj_kk(0)
        jj = ii_jj_kk(1)
        kk = ii_jj_kk(2)
        k0 = -1
        k1 = -1
        s  = [sp_in_0, sp_in_1, sp_in_2]
        call check_a2av_alignment(this%send_ids, ks, ii, jj, kk, k0, k1, s, this%order_a, all_lo_a, all_hi_a, &
                                  this%use_alltoallv, this%irank, jj_pos)
        s  = [sp_out_0, sp_out_1, sp_out_2]
        call check_a2av_alignment(this%recv_ids, kr, jj, ii, kk, k0, k1, s, this%order_b, all_lo_b, all_hi_b, &
                                  this%use_alltoallv, this%irank, temp_int)
        kk_pos = k0
      end block
    end if
    ! deallocate buffers
    deallocate(all_lo_a        , &
               all_hi_a        , &
               all_lo_b        , &
               all_hi_b        , &
               temp_buffer     , &
               temp_buffer_sort)
    block
      logical :: temp_bool
      temp_bool = this%use_alltoallv
      call MPI_Allreduce(temp_bool, this%use_alltoallv, 1, mpi_logical, mpi_land, mpi_comm_world, mpi_ierr)
    end block

    if (this%use_alltoallv) then
      call MPI_Comm_split(mpi_comm_world, kk_pos+1, jj_pos+1 , this%mpi_comm_kk, mpi_ierr)
      write(*,*) 'INFO: using alltoallv', this%irank
      ! note: mpi implementations rarely simplify unnecessary (local -> local) calls
    else
      this%mpi_comm_kk = mpi_comm_world
    !  write(*,*) 'INFO: using isend/irecv', this%irank
    end if

    call MPI_Comm_Rank(this%mpi_comm_kk, this%irank_comm_kk, mpi_ierr)
    call MPI_Allgather(this%irank_comm_kk   , 1, mpi_int,&
                       this%all_rank_comm_kk, 1, mpi_int, mpi_comm_world, mpi_ierr)
    ! v) configure local -> local transfers (if applicable)
    this%pos_start_localSend = -1
    this%pos_start_localRecv = -1
    do pos=0,ks
      if (this%send_ids(pos) == this%irank) then
        if (this%pos_start_localSend>=0) error stop 'ERROR: this%pos_start_localSend >= 0:' ! , this%pos_start_localSend
        this%pos_start_localSend = this%send_pos_start(pos)
        this%size_localSend      = this%send_sizes(pos)
      end if
    end do
    do pos=0,kr
      if (this%recv_ids(pos) == this%irank) then
        if (this%pos_start_localRecv>=0) error stop 'ERROR: this%pos_start_localRecv >= 0:' ! , this%pos_start_localRecv
        this%pos_start_localRecv = this%recv_pos_start(pos)
        this%size_localRecv      = this%recv_sizes(pos)
      end if
    end do

    ! vi) block variables
    allocate(this%send_i2b_s0(minval(this%send_i0_st):maxval(this%send_i0_st+this%send_lshape(:,0)-1)),&
             this%send_i2b_s1(minval(this%send_i1_st):maxval(this%send_i1_st+this%send_lshape(:,1)-1)),&
             this%send_i2b_s2(minval(this%send_i2_st):maxval(this%send_i2_st+this%send_lshape(:,2)-1)),&
             this%recv_i2b_s0(minval(this%recv_i0_st):maxval(this%recv_i0_st+this%recv_lshape(:,0)-1)),&
             this%recv_i2b_s1(minval(this%recv_i1_st):maxval(this%recv_i1_st+this%recv_lshape(:,1)-1)),&
             this%recv_i2b_s2(minval(this%recv_i2_st):maxval(this%recv_i2_st+this%recv_lshape(:,2)-1)),&

             this%send_i2b_i0a(minval(this%send_i0_st):maxval(this%send_i0_st+this%send_lshape(:,0)-1)),&
             this%send_i2b_i1a(minval(this%send_i1_st):maxval(this%send_i1_st+this%send_lshape(:,1)-1)),&
             this%send_i2b_i2a(minval(this%send_i2_st):maxval(this%send_i2_st+this%send_lshape(:,2)-1)),&
             this%recv_i2b_i0a(minval(this%recv_i0_st):maxval(this%recv_i0_st+this%recv_lshape(:,0)-1)),&
             this%recv_i2b_i1a(minval(this%recv_i1_st):maxval(this%recv_i1_st+this%recv_lshape(:,1)-1)),&
             this%recv_i2b_i2a(minval(this%recv_i2_st):maxval(this%recv_i2_st+this%recv_lshape(:,2)-1)),&

             aux1(max(ks,kr)+1,2),&
             aux2(max(ks,kr)+1,2))

    call mk_blocks_1d(this%send_i0_st, this%send_lshape(:,0), this%send_b_i0, this%send_i2b_s0, this%send_i2b_i0a, aux1, aux2)
    call mk_blocks_1d(this%send_i1_st, this%send_lshape(:,1), this%send_b_i1, this%send_i2b_s1, this%send_i2b_i1a, aux1, aux2)
    call mk_blocks_1d(this%send_i2_st, this%send_lshape(:,2), this%send_b_i2, this%send_i2b_s2, this%send_i2b_i2a, aux1, aux2)
    call mk_blocks_1d(this%recv_i0_st, this%recv_lshape(:,0), this%recv_b_i0, this%recv_i2b_s0, this%recv_i2b_i0a, aux1, aux2)
    call mk_blocks_1d(this%recv_i1_st, this%recv_lshape(:,1), this%recv_b_i1, this%recv_i2b_s1, this%recv_i2b_i1a, aux1, aux2)
    call mk_blocks_1d(this%recv_i2_st, this%recv_lshape(:,2), this%recv_b_i2, this%recv_i2b_s2, this%recv_i2b_i2a, aux1, aux2)

    deallocate(aux1,aux2)

    !$acc enter data copyin(this)
    !$acc enter data copyin(this%send_i2b_s0 , this%send_i2b_s1 , this%send_i2b_s2 ) async(stream)
    !$acc enter data copyin(this%send_i2b_i0a, this%send_i2b_i1a, this%send_i2b_i2a) async(stream)
    !$acc enter data copyin(this%recv_i2b_s0 , this%recv_i2b_s1 , this%recv_i2b_s2 ) async(stream)
    !$acc enter data copyin(this%recv_i2b_i0a, this%recv_i2b_i1a, this%recv_i2b_i2a) async(stream)
    !$acc wait(stream)

    allocate(this%send_all_n0(minval(this%send_b_i0):maxval(this%send_b_i0),&
                              minval(this%send_b_i1):maxval(this%send_b_i1),&
                              minval(this%send_b_i2):maxval(this%send_b_i2)))
    allocate(this%recv_all_n0(minval(this%recv_b_i0):maxval(this%recv_b_i0),&
                              minval(this%recv_b_i1):maxval(this%recv_b_i1),&
                              minval(this%recv_b_i2):maxval(this%recv_b_i2)))
    allocate(this%send_all_n1       , mold=this%send_all_n0)
    allocate(this%send_all_n2       , mold=this%send_all_n0)
    allocate(this%send_all_pos_start, mold=this%send_all_n0)
    allocate(this%recv_all_n1       , mold=this%recv_all_n0)
    allocate(this%recv_all_n2       , mold=this%recv_all_n0)
    allocate(this%recv_all_pos_start, mold=this%recv_all_n0)

    block
      integer :: abs0_reorder(0:2)
      if (((minval(abs0_reorder)==0).and.(maxval(abs0_reorder) == 2).and.(sum(abs0_reorder) == 3)).and.&
           (.not.this%single_device)) then
        abs0_reorder = abs_reorder
      else
        abs0_reorder = this%order_b(0:2)
      end if
      !$acc wait(stream)
      call diezdecomp_update_order_intermediate(this, abs0_reorder)
    end block

    block
      logical :: temp_bool
      temp_bool = (allow_autotune_reorder.and.(.not.this%single_device))
      call MPI_Allreduce(temp_bool, this%apply_autotune_reorder, 1, mpi_logical, mpi_land, mpi_comm_world, mpi_ierr)
    end block
    this%initialized            = .true.

  end subroutine

  subroutine diezdecomp_update_order_intermediate(this, abs_reorder, mode_send)
    integer                        ::  k, reorder_send(0:2), reorder_recv(0:2), abs_reorder(0:2)
    type(diezdecomp_props_transp)  ::  this
    logical, optional              ::  mode_send
    logical                        ::  to_send, to_recv

    to_send = .true.
    to_recv = .true.
    if (present(mode_send)) then
      to_send = (     mode_send)
      to_recv = (.not.mode_send)
    end if

    if (to_send.and.to_recv) this%chosen_reorder(0:2) = abs_reorder(0:2)

    if (.not.((minval(abs_reorder)==0).and.(maxval(abs_reorder) == 2).and.(sum(abs_reorder) == 3))) then
      error stop '((minval(abs_reorder) < 0).or.(minval(abs_reorder) == maxval(abs_reorder)))'
    end if
    if (to_send) then
      do k=0,2
        reorder_send(k) = diezdecomp_get_index_int2int(this%order_a, abs_reorder(k))
      end do
      call get_i2b_Mij(this%send_i2b_Mij, this%order_a, abs_reorder)
      call detailed_block_info(reorder_send     , this%send_pos_start, &
                               this%send_i0_st  , this%send_i1_st    , this%send_i2_st  , this%send_lshape,&
                               this%send_b_i0   , this%send_b_i1     , this%send_b_i2   , &
                               this%send_all_n0 , this%send_all_n1   , this%send_all_n2 , this%send_all_pos_start)
      this%send_mode_op_batched = -1 ! mode_op: 1 <-> 6
      this%send_mode_op_simul   = -1 ! mode_op: 1 <-> 6
      this%send_autotuned       = -1 ! [mode_op chosen -1/1, simultaneous = 1/batched = 2]
    end if
    if (to_recv) then
      do k=0,2
        reorder_recv(k) = diezdecomp_get_index_int2int(this%order_b, abs_reorder(k))
      end do
      call get_i2b_Mij(this%recv_i2b_Mij, this%order_b, abs_reorder)
      call detailed_block_info(reorder_recv     , this%recv_pos_start, &
                               this%recv_i0_st  , this%recv_i1_st    , this%recv_i2_st   , this%recv_lshape,&
                               this%recv_b_i0   , this%recv_b_i1     , this%recv_b_i2    , &
                               this%recv_all_n0 , this%recv_all_n1   , this%recv_all_n2  , this%recv_all_pos_start)
      this%recv_mode_op_batched = -1 ! mode_op: 1 <-> 6
      this%recv_mode_op_simul   = -1 ! mode_op: 1 <-> 6
      this%recv_autotuned       = -1 ! [mode_op chosen -1/1, simultaneous = 1/batched = 2]
    end if
  end subroutine

  subroutine diezdecomp_transp_execute(this, p_in, p_out, stream, buffer)
    implicit none
    type(diezdecomp_props_transp)   :: this
    real(rp), target                :: p_in(0:*), p_out(0:*)
    real(rp), target, optional      :: buffer(0:)
    real(rp), pointer, dimension(:) :: buf_in, buf_out
    integer(acc_handle_kind)        :: stream
    associate(s0 => this%numel_p_in  , &
              s1 => this%numel_p_out )
      if (present(buffer)) then
          buf_in (0:s0-1) => buffer(0 : s0   -1 )
          buf_out(0:s1-1) => buffer(s0:(s0+s1-1))
      else
          buf_in (0:s0-1) => p_out(0:s0-1)
          buf_out(0:s1-1) => p_in (0:s1-1)
      end if
    end associate

    call wrapper_batched_transpose(this, p_in, buf_in, .true., this%offset6_in, this%elapsed_time_reorder, stream)

    if (.not.this%chosen_backend) call comm_choose_backend(this, buf_out, buf_in, stream)
    if (this%use_alltoallv) then
      call comm_transpose_mpi_alltoallv(this, buf_out, buf_in, stream)
    else
      call comm_transpose_imode(this, buf_out, buf_in, this%use_nccl, stream)
    endif

    call wrapper_batched_transpose(this, p_out, buf_out, .false., this%offset6_out, this%elapsed_time_reorder, stream)

    if (this%apply_autotune_reorder) then
      block
        integer  :: i,j,k, abs_order(0:2),pos,mpi_ierr
        real(rp) :: elapsed_best, temp_real
        elapsed_best = LARGE
        do pos=0,5
          i = this%reorder_ijk(0,pos)
          j = this%reorder_ijk(1,pos)
          k = this%reorder_ijk(2,pos)
          temp_real = this%elapsed_time_reorder(pos)/(this%nproc*1._dp)
          call MPI_Allreduce(temp_real, this%elapsed_time_reorder(pos), 1, MPI_REAL_RP, mpi_sum, mpi_comm_world, mpi_ierr)
          if (this%elapsed_time_reorder(pos)<elapsed_best) then
            elapsed_best = this%elapsed_time_reorder(pos)
            abs_order    = [i,j,k]
            call diezdecomp_update_order_intermediate(this, abs_order)
            this%send_mode_op_batched   = this%best_send_mode_op_batched(pos)
            this%send_mode_op_simul     = this%best_send_mode_op_simul(pos)
            this%send_autotuned         = this%best_send_autotuned(pos)
            this%recv_mode_op_batched   = this%best_recv_mode_op_batched(pos)
            this%recv_mode_op_simul     = this%best_recv_mode_op_simul(pos)
            this%recv_autotuned         = this%best_recv_autotuned(pos)
          end if
        end do
        !$acc wait(stream)
        this%apply_autotune_reorder = .false.
      end block
    end if
  end subroutine

  ! ---------------------------------------- batched transpose implementations ------------------------------------------------

  recursive subroutine wrapper_batched_transpose(this, p_in, p_out, mode_fwd, offset6, elapsed_time_reorder, stream)
    implicit none
    real(rp)                      ::  p_in(0:*), p_out(0:*), elapsed_time_reorder(0:5), elapsed_time, temp_wtime,temp_x6(0:5)
    type(diezdecomp_props_transp) ::  this
    logical                       ::  mode_fwd
    integer                       ::  pos,iter1,abs_reorder(0:2), offset6(0:2,0:1)
    integer(acc_handle_kind)      ::  stream

    if (this%apply_autotune_reorder) then
      this%apply_autotune_reorder = .false.
      if (mode_fwd) elapsed_time_reorder = 0
      do pos=0,5
        abs_reorder(0:2) = this%reorder_ijk(0:2,pos)
        !$acc wait(stream)
        call diezdecomp_update_order_intermediate(this, abs_reorder, mode_fwd)
        call wrapper_batched_transpose(this, p_in, p_out, mode_fwd, offset6, temp_x6, stream)
        !$acc wait(stream)
        if (diezdecomp_mode_bench_avg) then ; elapsed_time = 0.
        else                                ; elapsed_time = LARGE ; end if
        do iter1 = 1,diezdecomp_autotune_mode_trials
          !$acc wait(stream)
          temp_wtime = MPI_Wtime()
          call wrapper_batched_transpose(this, p_in, p_out, mode_fwd, offset6, temp_x6, stream)
          !$acc wait(stream)
          temp_wtime = MPI_Wtime() - temp_wtime
          if (diezdecomp_mode_bench_avg) then
            elapsed_time = elapsed_time + temp_wtime/(diezdecomp_autotune_mode_trials*1._dp)
          else
            elapsed_time = min(elapsed_time, temp_wtime)
          end if
        end do
        elapsed_time_reorder(pos) = elapsed_time_reorder(pos) + elapsed_time
        if (mode_fwd) then
          this%best_send_mode_op_batched(pos) = this%send_mode_op_batched
          this%best_send_mode_op_simul(pos)   = this%send_mode_op_simul
          this%best_send_autotuned(pos)       = this%send_autotuned
        else
          this%best_recv_mode_op_batched(pos) = this%recv_mode_op_batched
          this%best_recv_mode_op_simul(pos)   = this%recv_mode_op_simul
          this%best_recv_autotuned(pos)       = this%recv_autotuned
        end if
      end do
      this%apply_autotune_reorder = .true.
    else
      if (mode_fwd) then
        call wrapper_batched_transpose_inner(this%send_sizes, p_in, p_out, .true., offset6, &
                                             this%send_mode_op_batched, this%send_lshape, &
                                             this%send_i0_st, this%send_i1_st, this%send_i2_st, &
                                             this%send_b_i0, this%send_b_i1, this%send_b_i2, &
                                             this%send_i2b_s0   , this%send_i2b_s1   , this%send_i2b_s2, &
                                             this%send_i2b_i0a, this%send_i2b_i1a, this%send_i2b_i2a, this%send_i2b_Mij, &
                                             this%send_all_n0, this%send_all_n1, this%send_all_n2, this%send_all_pos_start, &
                                             this%send_mode_op_simul, this%send_autotuned, &
                                             this%sp_in_0, this%sp_in_1, this%sp_in_2, &
                                             this%send_elapsed_batched_top, this%send_elapsed_simul_top, &
                                             this%send_elapsed_batched_modes, &
                                             this%send_elapsed_simul_modes, &
                                             stream)
      else
        call wrapper_batched_transpose_inner(this%recv_sizes, p_in, p_out, .false., offset6, &
                                             this%recv_mode_op_batched, this%recv_lshape, &
                                             this%recv_i0_st, this%recv_i1_st, this%recv_i2_st, &
                                             this%recv_b_i0, this%recv_b_i1, this%recv_b_i2, &
                                             this%recv_i2b_s0   , this%recv_i2b_s1   , this%recv_i2b_s2   , &
                                             this%recv_i2b_i0a, this%recv_i2b_i1a, this%recv_i2b_i2a, this%recv_i2b_Mij, &
                                             this%recv_all_n0, this%recv_all_n1, this%recv_all_n2, this%recv_all_pos_start, &
                                             this%recv_mode_op_simul, this%recv_autotuned, &
                                             this%sp_out_0, this%sp_out_1, this%sp_out_2, &
                                             this%recv_elapsed_batched_top, this%recv_elapsed_simul_top, &
                                             this%recv_elapsed_batched_modes, &
                                             this%recv_elapsed_simul_modes, &
                                             stream)
      end if
    end if
    !$acc wait(stream)
  endsubroutine

  subroutine wrapper_batched_transpose_inner(sizes, p_in, p_out, mode_fwd, offset6, mode_op_batched, loc_shape, &
                                                 i0_start, i1_start, i2_start, block_i0, block_i1, block_i2, &
                                                 i2b_s0, i2b_s1, i2b_s2, &
                                                 i2b_i0a, i2b_i1a, i2b_i2a, i2b_Mij, &
                                                 all_n0, all_n1, all_n2, all_pos_start, mode_op_simul, mode_transp, &
                                                 sp_in_0, sp_in_1, sp_in_2, &
                                                 elapsed_batched_top, elapsed_simul_top, elapsed_batched_modes, &
                                                 elapsed_simul_modes, &
                                                 stream)
    implicit none
    real(rp) :: p_in(0:*), p_out(0:*)
    logical  :: mode_fwd
    integer  :: iter0, iter1, pos, &
                i0_start(0:), i1_start(0:), i2_start(0:), block_i0(0:), block_i1(0:), block_i2(0:), &
                i2b_s0(0:) , i2b_s1(0:) , i2b_s2(0:) , offset6(0:2,0:1), &
                i2b_i0a(0:), i2b_i1a(0:), i2b_i2a(0:), i2b_Mij(0:2,0:2), lshape(0:2), &
                all_n0(0:,0:,0:), all_n1(0:,0:,0:), all_n2(0:,0:,0:), all_pos_start(0:,0:,0:), mode_transp, &
                mode_op_batched, sizes(0:), loc_shape(0:,0:), mode_op_simul, sp_in_0, sp_in_1, sp_in_2
    real(rp) :: elapsed_time, temp_wtime, elapsed_batched_top, elapsed_simul_top, elapsed_transp, &
                elapsed_batched_modes(0:) , elapsed_simul_modes(0:)
    integer(acc_handle_kind) :: stream
    logical  :: choose_mode_transp, choose_mode_batched, choose_mode_simul

    elapsed_transp      = LARGE
    choose_mode_transp  = ((mode_transp    <1).or.(mode_transp    >2))
    choose_mode_simul   = ((mode_op_simul  <1).or.(mode_op_simul  >6))
    choose_mode_batched = ((mode_op_batched<1).or.(mode_op_batched>6))

    if (choose_mode_simul.or.choose_mode_transp) then
      elapsed_batched_top = LARGE
      do iter0  = 3,6
        if (diezdecomp_mode_bench_avg) then ; elapsed_time = 0.
        else                                ; elapsed_time = LARGE ; end if
        do iter1 = 1,diezdecomp_autotune_mode_trials
          !$acc wait(stream)
          temp_wtime = MPI_Wtime()
          call flat_3d_transpose_simultaneous(p_in, p_out, mode_fwd, offset6, iter0, &
                                              i2b_s0, i2b_s1, i2b_s2, &
                                              i2b_i0a, i2b_i1a, i2b_i2a, i2b_Mij, &
                                              sp_in_0, sp_in_1, sp_in_2, stream)
          !$acc wait(stream)
          temp_wtime = MPI_Wtime() - temp_wtime
          if (diezdecomp_mode_bench_avg) then
            elapsed_time = elapsed_time + temp_wtime/(diezdecomp_autotune_mode_trials*1._dp)
          else
            elapsed_time = min(elapsed_time, temp_wtime)
          end if
        end do
        elapsed_simul_modes(iter0-1) = elapsed_time
        if (elapsed_time < elapsed_batched_top) then
          if (choose_mode_simul) mode_op_simul = iter0
          elapsed_batched_top                         = elapsed_time
        end if
      end do
      if (elapsed_batched_top < elapsed_transp) then
        if (choose_mode_transp) mode_transp = 1
        elapsed_transp                      = elapsed_batched_top
      end if
    end if

    if (choose_mode_batched.or.choose_mode_transp) then
      elapsed_simul_top = LARGE
      do iter0  = 3,6
        if (diezdecomp_mode_bench_avg) then ; elapsed_time = 0.
        else                                ; elapsed_time = LARGE ; end if
        do iter1 = 1,diezdecomp_autotune_mode_trials
          !$acc wait(stream)
          temp_wtime = MPI_Wtime()
          do pos=0,size(sizes,1)-1
            lshape = loc_shape(pos,:)
            call flat_3d_transpose_batched(p_in, p_out, mode_fwd, offset6, iter0, lshape, &
                                           i0_start(pos), i1_start(pos), &
                                           i2_start(pos), block_i0(pos), &
                                           block_i1(pos), block_i2(pos), &
                                           all_n0, all_n1, all_n2, all_pos_start, sp_in_0, sp_in_1, stream)
          end do
          !$acc wait(stream)
          temp_wtime = MPI_Wtime() - temp_wtime
          if (diezdecomp_mode_bench_avg) then
            elapsed_time = elapsed_time + temp_wtime/(diezdecomp_autotune_mode_trials*1._dp)
          else
            elapsed_time = min(elapsed_time, temp_wtime)
          end if
        end do
        elapsed_batched_modes(iter0-1) =  elapsed_time
        if (elapsed_time < elapsed_simul_top) then
          if (choose_mode_batched) mode_op_batched = iter0
          elapsed_simul_top                             = elapsed_time
        end if
      end do
      if (elapsed_simul_top < elapsed_transp) then
        if (choose_mode_transp) mode_transp = 2
        elapsed_transp                      = elapsed_simul_top
      end if
    end if

    if (mode_transp == 1) then
      call flat_3d_transpose_simultaneous(p_in, p_out, mode_fwd, offset6, mode_op_simul, &
                                          i2b_s0, i2b_s1, i2b_s2, &
                                          i2b_i0a, i2b_i1a, i2b_i2a, i2b_Mij, &
                                          sp_in_0, sp_in_1, sp_in_2, stream)
    else if (mode_transp == 2) then
      do pos=0,size(sizes,1)-1
        lshape = loc_shape(pos,:)
        call flat_3d_transpose_batched(p_in, p_out, mode_fwd, offset6, mode_op_batched, lshape, &
                                       i0_start(pos), i1_start(pos), &
                                       i2_start(pos), block_i0(pos), &
                                       block_i1(pos), block_i2(pos), &
                                       all_n0, all_n1, all_n2, all_pos_start, sp_in_0, sp_in_1, stream)
      end do
    else
      error stop 'ERROR: autotuning mode not recognized!'
    end if
  end subroutine

  subroutine flat_3d_transpose_batched(p_in, p_out, mode_fwd, offset6, mode_op, loc_shape, i0_start, i1_start, i2_start, &
                                       ib, jb, kb, all_n0, all_n1, all_n2, all_pos_start, sp_in_0, sp_in_1, stream)
    implicit none
    real(rp) :: p_out(0:*), p_in(0:*)
    logical  :: mode_fwd
    integer  :: i0_start, i1_start, i2_start, ib, jb, kb, pos_start, mode_op, na, nb, sp_in_0, sp_in_1, &
                i0_last, i1_last, i2_last, di, dj, dk, &
                all_n0(0:,0:,0:), all_n1(0:,0:,0:), all_n2(0:,0:,0:), all_pos_start(0:,0:,0:), &
                loc_shape(0:2), n0, n1, n2, i0,i1,i2, offset6(0:2,0:1)
    integer(acc_handle_kind) :: stream
    di        = offset6(0,0) ; dj = offset6(1,0) ; dk = offset6(2,0)
    na        =     sp_in_0 + offset6(0,0) + offset6(0,1)
    nb        = na*(sp_in_1 + offset6(1,0) + offset6(1,1))
    n0        = all_n0(ib,jb,kb)
    n1        = all_n1(ib,jb,kb)
    n2        = all_n2(ib,jb,kb)
    pos_start = all_pos_start(ib,jb,kb)

    i2_last   = i2_start + loc_shape(2)-1
    i1_last   = i1_start + loc_shape(1)-1
    i0_last   = i0_start + loc_shape(0)-1

    if (mode_fwd) then
      select case (mode_op)
      case (1)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i0 = i0_start, i0_last
        do i1 = i1_start, i1_last
        do i2 = i2_start, i2_last
          p_out(pos_start + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      case (2)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i0 = i0_start, i0_last
        do i2 = i2_start, i2_last
        do i1 = i1_start, i1_last
          p_out(pos_start + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      case (3)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i1 = i1_start, i1_last
        do i0 = i0_start, i0_last
        do i2 = i2_start, i2_last
          p_out(pos_start + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      case (4)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i1 = i1_start, i1_last
        do i2 = i2_start, i2_last
        do i0 = i0_start, i0_last
          p_out(pos_start + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      case (5)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i2 = i2_start, i2_last
        do i0 = i0_start, i0_last
        do i1 = i1_start, i1_last
          p_out(pos_start + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      case (6)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i2 = i2_start, i2_last
        do i1 = i1_start, i1_last
        do i0 = i0_start, i0_last
          p_out(pos_start + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      end select
    else
      select case (mode_op)
      case (1)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i0 = i0_start, i0_last
        do i1 = i1_start, i1_last
        do i2 = i2_start, i2_last
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(pos_start + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      case (2)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i0 = i0_start, i0_last
        do i2 = i2_start, i2_last
        do i1 = i1_start, i1_last
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(pos_start + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      case (3)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i1 = i1_start, i1_last
        do i0 = i0_start, i0_last
        do i2 = i2_start, i2_last
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(pos_start + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      case (4)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i1 = i1_start, i1_last
        do i2 = i2_start, i2_last
        do i0 = i0_start, i0_last
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(pos_start + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      case (5)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i2 = i2_start, i2_last
        do i0 = i0_start, i0_last
        do i1 = i1_start, i1_last
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(pos_start + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      case (6)
        !$acc parallel loop  collapse(3) default(present) async(stream)
        do i2 = i2_start, i2_last
        do i1 = i1_start, i1_last
        do i0 = i0_start, i0_last
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(pos_start + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      end select
    end if
  end subroutine


  subroutine flat_3d_transpose_simultaneous(p_in, p_out, mode_fwd, offset6, mode_op, &
                                            i2b_s0, i2b_s1, i2b_s2, &
                                            i2b_i0a, i2b_i1a, i2b_i2a, i2b_Mij, &
                                            sp_in_0, sp_in_1, sp_in_2, stream)
    implicit none
    real(rp) :: p_out(0:*), p_in(0:*)
    logical  :: mode_fwd
    integer  :: i0_start, i1_start, i2_start, i0_last, i1_last, i2_last, di, dj, dk, &
                              i0, i1, i2, mode_op, sp_in_0, sp_in_1, sp_in_2, na, nb, &
                              i2b_s0(0:)   , i2b_s1(0:)   , i2b_s2(0:), &
                              i2b_i0a(0:), i2b_i1a(0:), i2b_i2a(0:), i2b_Mij(0:2,0:2), &
                              n0, n1, n2, p, m00, m01, m02, &
                                             m10, m11, m12, &
                                             m20, m21, m22, s0,s1,s2, offset6(0:2,0:1), na_in, nb_in
    integer, parameter :: int1 = 1
    integer(acc_handle_kind) :: stream
    i0_start = 0
    i1_start = 0
    i2_start = 0
    i0_last  = sp_in_0-1
    i1_last  = sp_in_1-1
    i2_last  = sp_in_2-1
    di       = offset6(0,0) ; dj = offset6(1,0) ; dk = offset6(2,0)
    na       =     sp_in_0 + offset6(0,0) + offset6(0,1)
    nb       = na*(sp_in_1 + offset6(1,0) + offset6(1,1))
    na_in    = sp_in_0
    nb_in    = sp_in_0 * sp_in_1
    m00      = i2b_Mij(0,0)
    m01      = i2b_Mij(0,1)
    m02      = i2b_Mij(0,2)
    m10      = i2b_Mij(1,0)
    m11      = i2b_Mij(1,1)
    m12      = i2b_Mij(1,2)
    m20      = i2b_Mij(2,0)
    m21      = i2b_Mij(2,1)
    m22      = i2b_Mij(2,2)
    if (mode_fwd) then
      select case (mode_op)
      case (1)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i0 = i0_start, i0_last
        do i1 = i1_start, i1_last
        do i2 = i2_start, i2_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_out(p + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      case (2)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i0 = i0_start, i0_last
        do i2 = i2_start, i2_last
        do i1 = i1_start, i1_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_out(p + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      case (3)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i1 = i1_start, i1_last
        do i0 = i0_start, i0_last
        do i2 = i2_start, i2_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_out(p + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      case (4)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i1 = i1_start, i1_last
        do i2 = i2_start, i2_last
        do i0 = i0_start, i0_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_out(p + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      case (5)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i2 = i2_start, i2_last
        do i0 = i0_start, i0_last
        do i1 = i1_start, i1_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_out(p + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      case (6)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i2 = i2_start, i2_last
        do i1 = i1_start, i1_last
        do i0 = i0_start, i0_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_out(p + i0*n0 + i1*n1 + i2*n2) = p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb)
        end do
        end do
        end do
      end select
    else
      select case (mode_op)
      case (1)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i0 = i0_start, i0_last
        do i1 = i1_start, i1_last
        do i2 = i2_start, i2_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(p + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      case (2)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i0 = i0_start, i0_last
        do i2 = i2_start, i2_last
        do i1 = i1_start, i1_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(p + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      case (3)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i1 = i1_start, i1_last
        do i0 = i0_start, i0_last
        do i2 = i2_start, i2_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(p + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      case (4)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i1 = i1_start, i1_last
        do i2 = i2_start, i2_last
        do i0 = i0_start, i0_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(p + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      case (5)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i2 = i2_start, i2_last
        do i0 = i0_start, i0_last
        do i1 = i1_start, i1_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(p + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      case (6)
        !$acc parallel loop  collapse(3) default(present) private(i0,i1,i2,s0,s1,s2,n0,n1,n2,p) async(stream)
        do i2 = i2_start, i2_last
        do i1 = i1_start, i1_last
        do i0 = i0_start, i0_last
          s0  = i2b_s0(i0) ; s1  = i2b_s1(i1) ; s2  = i2b_s2(i2)
          n0 = (int1+m00*(s0-int1))*(int1+m01*(s1-int1))*(int1+m02*(s2-int1))
          n1 = (int1+m10*(s0-int1))*(int1+m11*(s1-int1))*(int1+m12*(s2-int1))
          n2 = (int1+m20*(s0-int1))*(int1+m21*(s1-int1))*(int1+m22*(s2-int1))
          p = (nb_in-n2)*i2b_i2a(i2) + (na_in*s2-n1)*i2b_i1a(i1) + (s1*s2-n0)*i2b_i0a(i0)
          p_in((i0+di)+(i1+dj)*na+(i2+dk)*nb) = p_out(p + i0*n0 + i1*n1 + i2*n2)
        end do
        end do
        end do
      end select
    end if
  end subroutine


  ! ---------------------------------------- transpose auxiliary functions ------------------------------------------------

  subroutine  transp_match_intervals(lo_a, hi_a, all_lo_b, all_hi_b, &
                                     temp_buffer, temp_buffer_sort, &
                                     i, nproc, loc_order_a, loc_order_b)
    implicit none
    integer :: lo_a(0:), hi_a(0:), all_lo_b(0:,0:), all_hi_b(0:,0:), &
               temp_buffer(0:,0:), temp_buffer_sort(0:,0:), temp_lo(0:2), temp_hi(0:2), &
               i,j,k, nproc, &
               loc_order_a(0:2), loc_order_b(0:2), reorder(0:2)
    do j=0,2
      reorder(j) = diezdecomp_get_index_int2int(loc_order_b, loc_order_a(j))
    end do
    i = 0
    do k=0,nproc-1
      do j=0,2
        temp_lo(j) = max(lo_a(j), all_lo_b(reorder(j),k))
        temp_hi(j) = min(hi_a(j), all_hi_b(reorder(j),k))
      end do
      if ((temp_lo(0) <= temp_hi(0)).and.&
          (temp_lo(1) <= temp_hi(1)).and.&
          (temp_lo(2) <= temp_hi(2))) then
        temp_buffer(i,0:2) = temp_lo
        temp_buffer(i,3)   = k
        temp_buffer(i,4:6) = temp_hi - temp_lo + 1
        i                  = i + 1
      end if
    end do
    i = i - 1
    block
      integer :: ip,p
      do ip=0,i
        p                 = temp_buffer(ip,0)
        temp_buffer(ip,0) = temp_buffer(ip,2)
        temp_buffer(ip,2) = p
      end do
    call mergesort_rows(temp_buffer, 0, i+0, temp_buffer_sort)
      do ip=0,i
        p                 = temp_buffer(ip,0)
        temp_buffer(ip,0) = temp_buffer(ip,2)
        temp_buffer(ip,2) = p
      end do
    end block
  end subroutine

  subroutine detailed_block_info(reorder, pos_start0, i0_start, i1_start, i2_start, loc_shape, block_i0, block_i1, block_i2, &
                                 all_n0, all_n1, all_n2, all_pos_start)
    implicit none
    integer :: pos, ib, jb, kb, reorder_inv(0:2), shape_aux(0:3), reorder(0:2), &
               pos_start0(0:), i0_start(0:), i1_start(0:), i2_start(0:), loc_shape(0:,0:), &
               block_i0(0:), block_i1(0:), block_i2(0:), &
               all_n0(0:,0:,0:), all_n1(0:,0:,0:), all_n2(0:,0:,0:), all_pos_start(0:,0:,0:)
      do pos=0,size(i0_start,1)-1
        reorder_inv = (/ reorder(2), reorder(1), reorder(0) /)
        shape_aux = (/ loc_shape(pos,reorder_inv(0)), loc_shape(pos,reorder_inv(1)), loc_shape(pos,reorder_inv(2)), 1 /)
        ib = block_i0(pos)
        jb = block_i1(pos)
        kb = block_i2(pos)
        all_n0(ib,jb,kb) = product(shape_aux(diezdecomp_get_index_int2int(reorder_inv,0)+1:))
        all_n1(ib,jb,kb) = product(shape_aux(diezdecomp_get_index_int2int(reorder_inv,1)+1:))
        all_n2(ib,jb,kb) = product(shape_aux(diezdecomp_get_index_int2int(reorder_inv,2)+1:))
        all_pos_start(ib,jb,kb) = pos_start0(pos) - (i0_start(pos)*all_n0(ib,jb,kb) + &
                                                     i1_start(pos)*all_n1(ib,jb,kb) + &
                                                     i2_start(pos)*all_n2(ib,jb,kb))
      end do
  end subroutine

  subroutine mk_blocks_1d(i0_start, shape_0, block_i0, i2b_s0, i2b_i0a, buffer, buffer_sort)
    implicit none
    integer :: i0_start(0:), shape_0(0:), block_i0(0:), i2b_s0(0:), i2b_i0a(0:), buffer(0:,0:), buffer_sort(0:,0:),&
    last,i,k,prev,myid
    last = size(i0_start,1)-1
    do k=0,last
      buffer(k,:) = [i0_start(k), k]
    end do
    call mergesort_rows(buffer, 0, last, buffer_sort)
    myid = -1
    prev = -1
    do k=0,last
      if (buffer(k,0) > prev) myid = myid + 1
      block_i0(buffer(k,1)) = myid
      prev = buffer(k,0)
    end do
    do k=0,last
      do i=i0_start(k),i0_start(k)+shape_0(k)-1
        i2b_s0(i)  = shape_0(k)
        i2b_i0a(i) = i0_start(k)
      end do
    end do
  end subroutine

  subroutine get_i2b_Mij(m, s, t)
    integer :: m(0:2,0:2), s(0:2), t(0:2),i,j,a,b
    m = 0
    do i=0,2
      a = diezdecomp_get_index_int2int(t,s(i))
      if ((a-1)>=0) then
        do j=0,a-1
          b = diezdecomp_get_index_int2int(s,t(j))
          m(i,b) = 1
        end do
      end if
    end do
  end subroutine

  ! ---------------------------------------- communication auxiliary routines ------------------------------------------------

  subroutine comm_choose_backend(this, buffer, p_out, stream)
    implicit none
    type(diezdecomp_props_transp) :: this
    real(rp)                 :: buffer(0:*)
    real(rp)                 :: p_out(0:*)
    integer                  :: pos, n_active, mode_best, mpi_ierr
    real(dp)                 :: elapsed_best, elapsed_alltoallv, elapsed_isendirecv, elapsed_nccl, temp_real, temp_wtime
    integer(acc_handle_kind) :: stream
    this%chosen_backend = .true.
    n_active = 0
    if (this%use_alltoallv)  n_active = n_active + 1
    if (this%use_isendirecv) n_active = n_active + 1
    if (this%use_nccl)       n_active = n_active + 1
    elapsed_best       = LARGE
    elapsed_alltoallv  = LARGE
    elapsed_isendirecv = LARGE
    elapsed_nccl       = LARGE
    mode_best          = 2
    if (n_active>1) then
      if (this%use_alltoallv) then
        do pos=1,diezdecomp_n_warmup ; call comm_transpose_mpi_alltoallv(this, buffer, p_out, stream) ; end do
        if (diezdecomp_mode_bench_avg) then ; elapsed_alltoallv = 0.
        else                                ; elapsed_alltoallv = LARGE ; end if
        do pos=1,diezdecomp_autotune_mode_trials
          !$acc wait(stream)
          temp_wtime = MPI_Wtime()
          call comm_transpose_mpi_alltoallv(this, buffer, p_out, stream)
          !$acc wait(stream)
          temp_wtime = MPI_Wtime() - temp_wtime
          if (diezdecomp_mode_bench_avg) then
            elapsed_alltoallv = elapsed_alltoallv + temp_wtime/(diezdecomp_autotune_mode_trials*1._dp)
          else
            elapsed_alltoallv = min(elapsed_alltoallv, temp_wtime)
          end if
        end do
        temp_real = elapsed_alltoallv/(this%nproc*1._dp)
        call MPI_Allreduce(temp_real, elapsed_alltoallv, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        if (elapsed_alltoallv<elapsed_best) mode_best=1
        elapsed_best            = min(elapsed_best, elapsed_alltoallv)
        this%elapsed_alltoallv  = elapsed_alltoallv
      end if

      if (this%use_isendirecv) then
        do pos=1,diezdecomp_n_warmup ; call comm_transpose_imode(this, buffer, p_out, .false., stream) ; end do
        if (diezdecomp_mode_bench_avg) then ; elapsed_isendirecv = 0.
        else                                ; elapsed_isendirecv = LARGE ; end if
        do pos=1,diezdecomp_autotune_mode_trials
          !$acc wait(stream)
          temp_wtime = MPI_Wtime()
          call comm_transpose_imode(this, buffer, p_out, .false., stream)
          !$acc wait(stream)
          temp_wtime = MPI_Wtime() - temp_wtime
          if (diezdecomp_mode_bench_avg) then
            elapsed_isendirecv = elapsed_isendirecv + temp_wtime/(diezdecomp_autotune_mode_trials*1._dp)
          else
            elapsed_isendirecv = min(elapsed_isendirecv, temp_wtime)
          end if
        end do
        temp_real = elapsed_isendirecv/(this%nproc*1._dp)
        call MPI_Allreduce(temp_real, elapsed_isendirecv, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        if (elapsed_isendirecv<elapsed_best) mode_best=2
        elapsed_best            = min(elapsed_best, elapsed_isendirecv)
        this%elapsed_isendirecv = elapsed_isendirecv
      end if
#if defined(_DIEZDECOMP_NCCL)
      if (this%use_nccl) then
        do pos=1,diezdecomp_n_warmup ; call comm_transpose_imode(this, buffer, p_out, .true., stream) ; end do
        if (diezdecomp_mode_bench_avg) then ; elapsed_nccl = 0.
        else                                ; elapsed_nccl = LARGE ; end if
        do pos=1,diezdecomp_autotune_mode_trials
          !$acc wait(stream)
          temp_wtime = MPI_Wtime()
          call comm_transpose_imode(this, buffer, p_out, .true., stream)
          !$acc wait(stream)
          temp_wtime = MPI_Wtime() - temp_wtime
          if (diezdecomp_mode_bench_avg) then
            elapsed_nccl = elapsed_nccl + temp_wtime/(diezdecomp_autotune_mode_trials*1._dp)
          else
            elapsed_nccl = min(elapsed_nccl, temp_wtime)
          end if
        end do
        temp_real = elapsed_nccl/(this%nproc*1._dp)
        call MPI_Allreduce(temp_real, elapsed_nccl, 1, MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm_world, mpi_ierr)
        if (elapsed_nccl<elapsed_best) mode_best=3
        elapsed_best      = min(elapsed_best, elapsed_nccl)
        this%elapsed_nccl = elapsed_nccl
      end if
#endif
      if (mode_best == 1) then ; this%use_alltoallv = .true. ; this%use_isendirecv = .false.; this%use_nccl = .false.;end if
      if (mode_best == 2) then ; this%use_alltoallv = .false.; this%use_isendirecv = .true. ; this%use_nccl = .false.;end if
      if (mode_best == 3) then ; this%use_alltoallv = .false.; this%use_isendirecv = .false.; this%use_nccl = .true. ;end if
      !write(*,*) this%irank, 'chosen backend', mode_best, elapsed_alltoallv, elapsed_isendirecv, elapsed_nccl
    end if
  end subroutine

  subroutine comm_transpose_mpi_alltoallv(this, buffer, p_out, stream)
    implicit none
    type(diezdecomp_props_transp) :: this
    real(rp)                 :: buffer(0:*)
    real(rp)                 :: p_out(0:*)
    integer                  :: mpi_ierr
    integer(acc_handle_kind) :: stream

    !$acc wait(stream)
    !$acc host_data use_device(buffer, p_out)
    call MPI_Alltoallv(p_out , this%send_sizes, this%send_pos_start, MPI_REAL_RP, &
                       buffer, this%recv_sizes, this%recv_pos_start, MPI_REAL_RP, &
                       this%mpi_comm_kk, mpi_ierr)
    !$acc end host_data
  end subroutine

  subroutine comm_transpose_imode(this, buffer, p_out, is_nccl, stream)
    implicit none
    type(diezdecomp_props_transp) :: this
    real(rp)                 :: buffer(0:*)
    real(rp)                 :: p_out(0:*)
    integer                  :: mpi_ierr, pos, ireq
    logical                  :: is_nccl
    integer(acc_handle_kind) :: stream
    !$acc wait(stream)
    ireq = 0
#if defined(_DIEZDECOMP_NCCL)
    if (is_nccl)  nccl_stat = ncclGroupStart()
#endif
    do pos=0,size(this%recv_ids,1)-1
      if (this%recv_ids(pos) /= this%irank) then
        call isendrecv_1d_buffer(.false., buffer, this%recv_pos_start(pos), this%recv_sizes(pos), &
                                 this%recv_ids(pos) , &
                                 this%recv_tags(pos), mpi_comm_world, this%all_request(ireq), is_nccl)
        ireq = ireq + 1
      end if
    end do
    do pos=0,size(this%send_ids,1)-1
      if (this%send_ids(pos) /= this%irank) then
        call isendrecv_1d_buffer(.true., p_out, this%send_pos_start(pos), this%send_sizes(pos), &
                                 this%send_ids(pos) , &
                                 this%send_tags(pos), mpi_comm_world, this%all_request(ireq), is_nccl)
        ireq = ireq + 1
      end if
    end do
    if ((this%pos_start_localSend>=0).and.(this%pos_start_localRecv>=0)) then
      call copy_1d_buffer(p_out, buffer, this%pos_start_localSend, this%pos_start_localRecv, this%size_localSend, stream)
    end if
#if defined(_DIEZDECOMP_NCCL)
    if (is_nccl)  then
      nccl_stat = ncclGroupEnd()
      mpi_ierr = cudaStreamSynchronize(nccl_stream)
      ireq = 0 ! nccl did not use ireq
    end if
#endif
    if (ireq>0) call mpi_waitall(ireq, this%all_request, this%all_status, mpi_ierr)
    !$acc wait(stream)
  end subroutine

  subroutine isendrecv_1d_buffer(mode_fwd,A,i,s,other,tag,comm,req,is_nccl)
    implicit none
    real(rp)          :: A(0:*)
    integer           :: i,s,other,tag,comm,req, i0,mpi_ierr
    logical           :: mode_fwd, is_mpi
    logical, optional :: is_nccl
    is_mpi = .true.
#if defined(_DIEZDECOMP_NCCL)
    if (present(is_nccl))  is_mpi = .not.is_nccl
#endif
    i0 = i
    !$acc host_data use_device(A)
    if (is_mpi) then
      if (mode_fwd) then ; call MPI_ISend(A(i0),s, MPI_REAL_RP, other, tag, comm, req, mpi_ierr)
      else               ; call MPI_IRecv(A(i0),s, MPI_REAL_RP, other, tag, comm, req, mpi_ierr)
      end if
#if defined(_DIEZDECOMP_NCCL)
    else
      if (mode_fwd) then ; nccl_stat=ncclSend(A(i0),s, ncclDouble, other, nccl_comm, nccl_stream)
      else               ; nccl_stat=ncclRecv(A(i0),s, ncclDouble, other, nccl_comm, nccl_stream)
      end if
#endif
    end if
    !$acc end host_data
  end subroutine

  ! ---------------------------------------- diezdecomp get_ranks -------------------------------------------

  subroutine diezdecomp_fill_mpi_ranks(lo, mpi_ranks, flat_mpi_ranks, shape_mpi_ranks, irank, nproc)
    implicit none
    integer              :: lo(0:2), mpi_ierr, irank, nproc, i, j, k, pos, last, shape_mpi_ranks(0:2)
    integer, allocatable :: buffer(:), full_arr(:,:), buffer_sort(:,:), mpi_ranks(:,:,:), flat_mpi_ranks(:,:)

    ! Notes: 1) lo(0:2) must be the global (i,j,k) position of "irank" in the 3-D pencil decomposition (xyz).
    !        2) If lo(0:2) has an order different from (ijk or xyz) (e.g., yxz), then it must be converted before
    !           calling this subroutine (by an API ifle for instance).

    allocate(buffer(0:(nproc*3-1)), full_arr(0:(nproc-1),0:3), buffer_sort(0:(nproc-1),0:3))

    call MPI_Allgather(lo    , 3, mpi_int,&
                       buffer, 3, mpi_int,&
                       mpi_comm_world, mpi_ierr)
    do   i = 0, nproc-1
      do j = 0, 2
        full_arr(i,j) = buffer(3*i + j)
      end do
      full_arr(i,3) = i ! "i" is irank for this buffer position
    end do

    ! check if irank matches full_arr data
    if (maxval(abs(full_arr(irank,0:2)-lo))>0) error stop  'ERROR! full_arr(irank,0:2)' ! , full_arr(irank+1,1:3), lo

    ! sort full_arr (by position)
    call mergesort_rows(full_arr, 0, nproc-1, buffer_sort)

    ! transform "lo" into ijk
    do j=0,2
      last = -1
      pos  = -1
      do i=0,nproc-1
        if (full_arr(i,j) < last) then
          pos = 0
        else if (full_arr(i,j) > last) then
          pos = pos + 1
        end if
        last          = full_arr(i,j)
        full_arr(i,j) = pos
      end do
    end do

    if (allocated(mpi_ranks)) deallocate(mpi_ranks)
    allocate(mpi_ranks(0:maxval(full_arr(:,0)), 0:maxval(full_arr(:,1)), 0:maxval(full_arr(:,2))))

    if (nproc /= (size(mpi_ranks,1)*size(mpi_ranks,2)*size(mpi_ranks,3))) error stop  'ERROR! nproc /= prod(shape(mpi_ranks))'

    do i = 0, nproc-1
      mpi_ranks(full_arr(i,0),full_arr(i,1),full_arr(i,2)) = full_arr(i,3)
    end do

    deallocate(buffer, full_arr, buffer_sort)

    ! complete shape_mpi_ranks, flat_mpi_ranks
      do i=0,2
        shape_mpi_ranks(i) = size(mpi_ranks,i+1)
      end do

      if (allocated(flat_mpi_ranks)) deallocate(flat_mpi_ranks)
      allocate(flat_mpi_ranks(0:nproc-1,0:2))

      do     i=0,size(mpi_ranks,1)-1
        do   j=0,size(mpi_ranks,2)-1
          do k=0,size(mpi_ranks,3)-1
            pos                   = mpi_ranks(i,j,k)
            flat_mpi_ranks(pos,0) = i
            flat_mpi_ranks(pos,1) = j
            flat_mpi_ranks(pos,2) = k
          end do
        end do
      end do
  end subroutine

  function diezdecomp_get_rank_id(ranks_ii, flat_ranks_ii, shape_ranks_ii, irank, ii, d) result(pos)
    implicit none
    integer, intent(in) :: ranks_ii(0:,0:,0:), flat_ranks_ii(0:,0:),shape_ranks_ii(0:2), irank, ii, d
    integer             :: i,j,k, pos

    ! call check_bounds(ii,  0, 2)
    ! call check_bounds( d, -1, 1)

    i = flat_ranks_ii(irank,0)
    j = flat_ranks_ii(irank,1)
    k = flat_ranks_ii(irank,2)
    if (ii == 0) i = modulo(i + d,shape_ranks_ii(0))
    if (ii == 1) j = modulo(j + d,shape_ranks_ii(1))
    if (ii == 2) k = modulo(k + d,shape_ranks_ii(2))
    pos = ranks_ii(i,j,k)
  end function

  function get_position_in_3darr(flat_ranks_ii, irank, ii) result(pos)
    implicit none
    integer, intent(in) :: flat_ranks_ii(0:,0:), irank, ii
    integer             :: pos
    pos = flat_ranks_ii(irank, ii)
  end function

  ! ---------------------------------------- general purpose routines ------------------------------------------------

  subroutine copy_1d_buffer(A, B, ia, ib, sa, stream)
    implicit none
    real(rp) :: A(0:*), B(0:*)
    integer  :: i, ia, ib, sa, di
    integer(acc_handle_kind) :: stream

    if (sa>0) then
      di = ia - ib
      !$acc parallel loop  collapse(1) default(present) async(stream)
      do i = ib, ib+sa-1
        B(i) = A(i+di)
      end do
    end if
  end subroutine

  ! subroutine check_bounds(x, a, b, key0)
  !   implicit none
  !   integer           :: x, a, b, key
  !   integer, optional :: key0 ! key0 is a hint/identifier about the source call
  !   key = -1
  !   if (present(key0)) key = key0
  !   if (.not.((a <= x).and.(x <= b))) then
  !     error stop 'ERROR: out-of-bounds' ! , x, a, b, key
  !   end if
  ! end subroutine

  function diezdecomp_get_index_int2int(A, ref) result(pos)
    implicit none
    ! implementation of python function: pos = list_x.index(value_y)
    integer :: A(0:), ref, pos
    pos = findloc(A,ref,dim=1) - 1
    if (pos == -1) error stop 'ERROR: diezdecomp_get_index_int2int'
  end function

  function is_greater(a, b) result(res)
    implicit none
    logical :: res
    integer :: a(:), b(:), i
    res = .false.
    do i=1,max(size(a,1),size(b,1))
      if (a(i) > b(i)) then
        res = .true.
        return
      else if (a(i) < b(i)) then
        res = .false.
        return
      end if
    end do
    ! reached this line ---> all equal ---> is_greater=false --> ok
  end function

  recursive subroutine mergesort_rows(A,i,j,buffer)
    implicit none
    integer :: i,j,L,mid,i0,j0,k,k2,pos,buffer(0:,:),A(0:,:)
    logical :: aux_bool
    L = j-i+1
    if (L>1) then
      if (L==2) then
        if (is_greater(A(i,:),A(j,:))) then
          buffer(0,:) = A(i,:)
          A(i,:)      = A(j,:)
          A(j,:)      = buffer(0,:)
        end if
      else
        mid = (i+j)/2
        call mergesort_rows(A,     i, mid, buffer)
        call mergesort_rows(A, mid+1,   j, buffer)
        i0 = i
        j0 = mid+1
        k  = 0
        do while (((i0<=mid).or.(j0<=j)).and.(k<=L))
          aux_bool = (j0>j)
          if (.not.aux_bool) then
            aux_bool = (i0<=mid)
            if (aux_bool) aux_bool=aux_bool.and.(.not.is_greater(A(i0,:),A(j0,:)))
          end if
          if (aux_bool) then !((j0>j).or.((i0<=mid).and.(.not.is_greater(A(i0,:),A(j0,:))))) then
            buffer(k,:) = A(i0,:)
            i0          = i0 + 1
          else
            buffer(k,:) = A(j0,:)
            j0          = j0 + 1
          end if
          k = k + 1
        end do
        pos = i
        do k2=0,k-1
          A(pos,:) = buffer(k2,:)
          pos      = pos + 1
        end do
      end if
    end if
  end subroutine

  ! ---------------------- autotuning summaries ----------------------
  subroutine diezdecomp_summary_transp_autotuning(this)
    implicit none
    type(diezdecomp_props_transp) :: this
    integer                       :: pos, use_alltoallv, use_isendirecv, use_nccl, temp_count
    character(len=65)             :: txt
    character(len=84)             :: stats84
    associate(irank => this%irank, nproc => this%nproc)
      if (irank==0) write(*,'(A57)') '------------------- Summary Transpose -------------------'

      call check_unique_int(this%chosen_reorder(0))
      call check_unique_int(this%chosen_reorder(1))
      call check_unique_int(this%chosen_reorder(2))
      if (irank==0) write(*,'(A26,3I2)')            '    Loop order chosen  :  ', this%chosen_reorder
      do pos=0,5
        call get_stats_value(txt, this%elapsed_time_reorder(pos))
        if (irank==0) write(*,'(A26,3I2,A2,A65)') '    Elapsed time, order: (', this%reorder_ijk(:,pos),') ', txt
      end do

      use_alltoallv   =  0
      use_isendirecv  =  0
      use_nccl        =  0
      if (this%use_alltoallv  ) use_alltoallv   =  1
      if (this%use_isendirecv ) use_isendirecv  =  1
      if (this%use_nccl       ) use_nccl        =  1

      call check_unique_int(use_alltoallv)
      call check_unique_int(use_isendirecv)
      call check_unique_int(use_nccl)

      if (this%single_device) then
        if (irank==0)                       write(*,'(A47)') '    MPI mode chosen      : None (single_device)'
      else
        if ((use_alltoallv+use_isendirecv+use_nccl) /= 1) error stop 'ERROR: (use_alltoallv+use_isendirecv+use_nccl) /= 1)'
        if (this%use_alltoallv .and.(irank==0))  write(*,'(A38)') '    MPI mode chosen      : AlltoAllV  '
        if (this%use_isendirecv.and.(irank==0))  write(*,'(A38)') '    MPI mode chosen      : ISend/IRecv'
        if (this%use_nccl      .and.(irank==0))  write(*,'(A38)') '    MPI mode chosen      : NCCL       '
      end if

      call get_stats_value(txt, this%elapsed_alltoallv)
      if (irank==0) write(*,'(A45,A65)') '    Elapsed time       : AlltoAllV (sync)    ', txt

      call get_stats_value(txt, this%elapsed_isendirecv)
      if (irank==0) write(*,'(A45,A65)') '    Elapsed time       : ISend/IRecv (async) ', txt

      call get_stats_value(txt, this%elapsed_nccl)
      if (irank==0) write(*,'(A45,A65)') '    Elapsed time       : NCCL        (async) ', txt

      if (irank==0) write(*,*) ''

      call get_count_int(temp_count, this%send_autotuned, 1)
      call get_stats_modes(stats84, this%send_mode_op_simul)
      if (irank==0) write(*,'(A48,I6,A84)')         '    Transpose mode Simultaneous chosen (send) : ',temp_count,stats84
      call get_stats_value(txt, this%send_elapsed_simul_top)
      if (irank==0) write(*,'(A45,A65)')            '    Elapsed time: Simultaneous (send) (best) ', txt
      if (irank==0) write(*,*) ''
      do pos=0,5
        call get_stats_value(txt, this%send_elapsed_simul_modes(pos))
        if (irank==0) write(*,'(A45,I2,A2,A65)')    '    Elapsed time, Simultaneous (send) (mode: ', pos+1,') ', txt
      end do
      if (irank==0) write(*,*) ''

      call get_count_int(temp_count, this%send_autotuned, 2)
      call get_stats_modes(stats84, this%send_mode_op_batched)
      if (irank==0) write(*,'(A48,I6,A84)')         '    Transpose mode Batched      chosen (send) : ',temp_count,stats84
      call get_stats_value(txt, this%send_elapsed_batched_top)
      if (irank==0) write(*,'(A45,A65)')            '    Elapsed time: Batched      (send) (best) ', txt
      if (irank==0) write(*,*) ''
      do pos=0,5
        call get_stats_value(txt, this%send_elapsed_batched_modes(pos))
        if (irank==0) write(*,'(A45,I2,A2,A65)')    '    Elapsed time, Batched      (send) (mode: ', pos+1,') ', txt
      end do
      if (irank==0) write(*,*) ''

      call get_count_int(temp_count, this%recv_autotuned, 1)
      call get_stats_modes(stats84, this%recv_mode_op_simul)
      if (irank==0) write(*,'(A48,I6,A84)')         '    Transpose mode Simultaneous chosen (recv) : ',temp_count,stats84
      call get_stats_value(txt, this%recv_elapsed_simul_top)
      if (irank==0) write(*,'(A45,A65)')            '    Elapsed time: Simultaneous (recv) (best) ', txt
      if (irank==0) write(*,*) ''
      do pos=0,5
        call get_stats_value(txt, this%recv_elapsed_simul_modes(pos))
        if (irank==0) write(*,'(A45,I2,A2,A65)')    '    Elapsed time, Simultaneous (recv) (mode: ', pos+1,') ', txt
      end do
      if (irank==0) write(*,*) ''

      call get_count_int(temp_count, this%recv_autotuned, 2)
      call get_stats_modes(stats84, this%recv_mode_op_batched)
      if (irank==0) write(*,'(A48,I6,A84)')         '    Transpose mode Batched      chosen (recv) : ',temp_count,stats84
      call get_stats_value(txt, this%recv_elapsed_batched_top)
      if (irank==0) write(*,'(A45,A65)')            '    Elapsed time: Batched      (recv) (best) ', txt
      if (irank==0) write(*,*) ''
      do pos=0,5
        call get_stats_value(txt, this%recv_elapsed_batched_modes(pos))
        if (irank==0) write(*,'(A45,I2,A2,A65)')    '    Elapsed time, Batched      (recv) (mode: ', pos+1,') ', txt
      end do
      if (irank==0) write(*,*) ''
    end associate
  end subroutine

  subroutine get_stats_modes(txt,val)
    character(len=84) :: txt
    integer           :: val, i, p(1:6), other, temp_int, mpi_ierr
    p     = 0
    other = 0
    if ((val<1).or.(6<val)) then
      other = other + 1
    else
        p(val) = p(val) + 1
    end if
    do i=1,6
      temp_int = p(i)
      call MPI_Allreduce(temp_int , p(i), 1, mpi_int, mpi_sum, mpi_comm_world, mpi_ierr)
    end do
    write(txt,'(A5,I6,A6,I6,A6,I6,A6,I6,A6,I6,A6,I6,A6,I6,A1)') ' (1: ',p(1),' , 2: ',p(2),' , 3: ',p(3),&
                                                               ' , 4: ',p(4),' , 5: ',p(5),' , 6: ',p(6),' , ?: ',other,')'
  end subroutine

  subroutine diezdecomp_summary_halo_autotuning(this)
    implicit none
    type(diezdecomp_props_halo) :: this
    character(len=65)           :: txt
    associate(irank => this%irank, nproc => this%nproc, nh_dir => this%nh_dir)
      if (irank==0) write(*,*) ''
      if (irank==0) write(*,'(A57)') '-------------------- Halo Autotuning --------------------'

      call check_unique_int(this%halo_mpi_mode)

      if (this%halo_mpi_mode==1) then
        if (irank==0) write(*,'(A44)')   '    MPI mode chosen    : MPI_SendRecv (sync)'
      else if (this%halo_mpi_mode==2) then
        if (irank==0) write(*,'(A44)')   '    MPI mode chosen    : ISend/IRecv (async)'
      else
        if (.not.((nproc== 1).or.(nh_dir==0))) error stop 'ERROR: this%halo_mpi_mode out-of-range'
        if (irank==0)   write(*,'(A37,I6,A10,I6,A1)') '    MPI mode chosen    : None (nproc=',nproc,') (nh_dir=',nh_dir,')'
      end if
      if (irank==0) write(*,*) ''
      call get_stats_value(txt, this%elapsed_sendrecv)
      if (irank==0) write(*,'(A45,A65)') '    Elapsed time       : MPI_SendRecv (sync) ', txt
      call get_stats_value(txt, this%elapsed_isendirecv)
      if (irank==0) write(*,'(A45,A65)') '    Elapsed time       : ISend/IRecv (async) ', txt
      if (irank==0) write(*,*) ''
      call check_unique_int(this%autotuned_pack)

      if (this%autotuned_pack==2) then
        if (irank==0) write(*,'(A40)')    '    Packing mode chosen: pack_slices_2d '
      else if (this%autotuned_pack==3) then
        if (irank==0) write(*,'(A40)')    '    Packing mode chosen: pack_slices_3d '
      else
        if (nh_dir == 0) then
          if (irank==0) write(*,'(A40)')  '    Packing mode chosen: None (nh_dir=0)'
        else
          error stop 'ERROR: this%autotuned_pack out-of-range'
        end if
      end if
      if (irank==0) write(*,*) ''
      call get_stats_value(txt, this%wtime_2)
      if (irank==0) write(*,'(A45,A65)') '    Elapsed time       : pack_slices_2d      ', txt
      call get_stats_value(txt, this%wtime_3)
      if (irank==0) write(*,'(A45,A65)') '    Elapsed time       : pack_slices_3d      ', txt
      if (irank==0) write(*,*) ''
    end associate
  end subroutine

  subroutine check_unique_int(val)
    integer :: val, mpi_ierr, val_min, val_max
    call MPI_Allreduce(val, val_min, 1, mpi_int, mpi_min, mpi_comm_world, mpi_ierr)
    call MPI_Allreduce(val, val_max, 1, mpi_int, mpi_max, mpi_comm_world, mpi_ierr)
    if (val_min /= val_max) error stop 'ERROR: check_unique_int (val_min /= val_max)'
  end subroutine

  subroutine get_stats_value(txt, val)
    real(rp)          ::  val, val_min, val_max, val_avg
    integer           ::  mpi_ierr, total, int1
    character(len=65) :: txt
    int1 = 1
    call MPI_Allreduce(val , val_min, 1, MPI_REAL_RP, mpi_min, mpi_comm_world, mpi_ierr)
    call MPI_Allreduce(val , val_max, 1, MPI_REAL_RP, mpi_max, mpi_comm_world, mpi_ierr)
    call MPI_Allreduce(val , val_avg, 1, MPI_REAL_RP, mpi_sum, mpi_comm_world, mpi_ierr)
    call MPI_Allreduce(int1, total  , 1, mpi_int    , mpi_sum, mpi_comm_world, mpi_ierr)
    val_avg  =  val_avg / (total+0.d0)
    write(txt,'(A6,E14.6,A8,E14.6,A8,E14.6,A1)') '(min: ',val_min,' , avg: ',val_avg,' , max: ',val_max,')'
  end subroutine

  subroutine get_count_int(n, val, target_val)
    integer :: n,val,target_val,temp_int, mpi_ierr
    temp_int = 0
    n        = 0
    if (val==target_val) temp_int = 1
    call MPI_Allreduce(temp_int, n, 1, mpi_int, mpi_sum, mpi_comm_world, mpi_ierr)
  end subroutine

  subroutine check_a2av_alignment(my_ids, n, ii, jj, kk, k0, k1, s, order, all_lo, all_hi, use_alltoallv, irank, jj_pos)
    integer :: i_loc, j_loc, k_loc, j_prev, n, ii, jj, kk, k0, k1, other, &
               s(0:2), order(0:2), all_lo(0:,0:), all_hi(0:,0:), my_ids(0:), pos, irank, jj_pos
    logical :: use_alltoallv
    i_loc   =  diezdecomp_get_index_int2int(order, ii)
    j_loc   =  diezdecomp_get_index_int2int(order, jj)
    k_loc   =  diezdecomp_get_index_int2int(order, kk)
    j_prev  =  -1
    if (use_alltoallv) then
      do pos=0,n
        other = my_ids(pos)
        if      (all_lo(i_loc,other) /= 0             ) use_alltoallv = .false.
        if      (all_hi(i_loc,other) /= (s(i_loc)-1)  ) use_alltoallv = .false.
        if      (all_lo(j_loc,other)<=j_prev          ) use_alltoallv = .false.
        j_prev = all_hi(j_loc,other)
        if (k0<0) k0 = all_lo(k_loc,other)
        if (k1<0) k1 = all_hi(k_loc,other)
        if (all_lo(k_loc,other) /= k0            ) use_alltoallv = .false.
        if (all_hi(k_loc,other) /= k1            ) use_alltoallv = .false.
        if (other == irank) jj_pos = j_prev
      end do
    end if
  end subroutine
end module
