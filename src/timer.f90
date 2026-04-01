! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
!
! a simple timer, see https://github.com/p-costa/first-timer
!
module mod_timer
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use mpi
#if defined(_OPENACC)
  use openacc, only: acc_handle_kind
  !@cuf use cudafor
#endif
#if defined(_USE_NVTX)
  use mod_nvtx
#endif
  implicit none
  private
  public :: timer_tic,timer_toc,timer_print,timer_cleanup,timer_drain
  !
  logical, parameter :: GPU_DEFAULT_SYNC = .true.
  integer, parameter :: max_name_len = 50
  character(max_name_len), allocatable :: timer_names(:)
  integer , allocatable :: timer_counts(:)
  real(dp), allocatable :: timer_tictoc(:),timer_elapsed_acc(:), &
                                           timer_elapsed_min(:), &
                                           timer_elapsed_max(:)
  integer :: ntimers = 0
#if defined(_OPENACC)
  integer, parameter :: TIMER_MAX_PENDING = 2048
  integer :: ce_slot_counter = 0
  integer :: ce_drain_n = 0
  integer :: ce_drain_slot(    TIMER_MAX_PENDING)
  integer :: ce_drain_timer_idx(TIMER_MAX_PENDING)
  integer, allocatable :: timer_tic_slot(:)
  !@cuf type(cudaEvent), allocatable :: ce_pool_start(:), ce_pool_end(:)
  !@cuf logical :: ce_pool_initialized = .false.
#endif
contains
  subroutine timer_print(myid_arg)
    use, intrinsic :: iso_fortran_env, only: stdo => output_unit
    integer , parameter :: MYID_PRINT = 0
    logical , parameter :: is_verbose_level_1 = .false.
    logical , parameter :: is_verbose_level_2 = .false.
    integer , intent(in), optional :: myid_arg
    real(dp), allocatable :: timing_results_acc(:,:), &
                             timing_results_min(:,:), &
                             timing_results_max(:,:)
    integer  :: i,myid,nproc,ierr
    !
    if(present(myid_arg)) then
      myid = myid_arg
    else
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    end if
#if defined(_OPENACC) && defined(_TIMER)
    call timer_drain()
#endif
    allocate(timing_results_acc(ntimers,3), &
             timing_results_min(ntimers,3), &
             timing_results_max(ntimers,3))
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call MPI_ALLREDUCE(timer_elapsed_acc(:),timing_results_acc(:,1),ntimers,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(timer_elapsed_acc(:),timing_results_acc(:,2),ntimers,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(timer_elapsed_acc(:),timing_results_acc(:,3),ntimers,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    timing_results_acc(:,3) = timing_results_acc(:,3)/nproc
    call MPI_ALLREDUCE(timer_elapsed_min(:),timing_results_min(:,1),ntimers,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(timer_elapsed_min(:),timing_results_min(:,2),ntimers,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(timer_elapsed_min(:),timing_results_min(:,3),ntimers,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    timing_results_min(:,3) = timing_results_min(:,3)/nproc
    call MPI_ALLREDUCE(timer_elapsed_max(:),timing_results_max(:,1),ntimers,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(timer_elapsed_max(:),timing_results_max(:,2),ntimers,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(timer_elapsed_max(:),timing_results_max(:,3),ntimers,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    timing_results_max(:,3) = timing_results_max(:,3)/nproc
    !
    if(myid == MYID_PRINT) then
      write(stdo,*) ''
      write(stdo,*) '*** timing results [s] ***'
      write(stdo,*) ''
      if(nproc == 1.or..not.is_verbose_level_1) then
        do i = 1,ntimers
          write(stdo,'(3A)'      ) 'Label: "',trim(timer_names(i)), '"'
          write(stdo,'(A,3E15.7)') 'Elapsed time:', timing_results_acc(i,3:3)
          write(stdo,'(A,I7)'    ) 'Number of calls:', timer_counts(i)
          write(stdo,'(A,1E15.7)') 'Average elapsed time per task (per call average):',timing_results_acc(i,3:3)/timer_counts(i)
          if(is_verbose_level_2) then
            write(stdo,'(A,1E15.7)') 'Average elapsed time per task (per call minimum):',timing_results_min(i,3:3)
            write(stdo,'(A,1E15.7)') 'Average elapsed time per task (per call maximum):',timing_results_max(i,3:3)
          end if
          write(stdo,*) ''
        end do
      else
        do i = 1,ntimers
          write(stdo,'(3A)'      ) 'Label: "',trim(timer_names(i)), '"'
          write(stdo,'(A,3E15.7)') 'Maximum, minimum, average elapsed time per task:', timing_results_acc(i,1:3)
          write(stdo,'(A,I7)'    ) 'Number of calls:', timer_counts(i)
          write(stdo,'(A,3E15.7)') 'Maximum, minimum, average elapsed time per task (per call average):', &
                                    timing_results_acc(i,1:3)/timer_counts(i)
          if(is_verbose_level_2) then
            write(stdo,'(A,3E15.7)') 'Maximum, minimum, average elapsed time per task (per call minimum):', &
                                      timing_results_min(i,1:3)
            write(stdo,'(A,3E15.7)') 'Maximum, minimum, average elapsed time per task (per call maximum):', &
                                      timing_results_max(i,1:3)
          end if
          write(stdo,*) ''
        end do
      end if
    end if
  end subroutine timer_print
  subroutine timer_tic(timer_name,nvtx_id_fix,nvtx_color,nvtx_id_inc,nvtx_gpu_stream,cuda_stream)
    !@cuf use cudafor
    character(*), intent(in) :: timer_name
    integer         , intent(in   ), optional :: nvtx_id_fix     ! if <= 0, only label and no color
    character(len=1), intent(in   ), optional :: nvtx_color      ! g/b/y/m/c/r/w following matplotlib's convention
    integer         , intent(inout), optional :: nvtx_id_inc     ! to increment the id, e.g.: call timer_tic(name,nvtx_id_inc=i_nvtx)
    integer         , intent(in   ), optional :: nvtx_gpu_stream ! to optionally sync host/device over a stream/queue (asynchronous if < 0)
#if defined(_OPENACC)
    integer(acc_handle_kind), intent(in), optional :: cuda_stream
    integer(acc_handle_kind) :: stream_to_use
#else
    integer, intent(in), optional :: cuda_stream ! unused argument to avoid CPP conditionals in the subroutine definition
#endif
    integer :: idx,nvtx_id
    logical :: is_nvtx,is_gpu_sync
    !@cuf integer :: istat
    is_nvtx = .false.
    if(     present(nvtx_id_inc)) then
      nvtx_id = nvtx_id_inc
      if(nvtx_id == huge(1)) nvtx_id_inc = 0 ! avoid overflow
      nvtx_id_inc = nvtx_id_inc + 1
      is_nvtx = .true.
    else if(present(nvtx_id_fix)) then
      nvtx_id = nvtx_id_fix
      is_nvtx = .true.
    else if(present(nvtx_color )) then
      is_nvtx = .true.
    end if
    is_gpu_sync = GPU_DEFAULT_SYNC
    if(present(nvtx_gpu_stream)) then
      if(nvtx_gpu_stream < 0) then
        is_gpu_sync = .false.
      end if
    end if
    if(is_gpu_sync) then
      if(.not.present(nvtx_gpu_stream)) then
#if   defined(_OPENACC)
        !$acc wait
#elif defined(_CUDA)
        !@cuf istat=cudaDeviceSynchronize()
#endif
      else
#if   defined(_OPENACC)
        !$acc wait(nvtx_gpu_stream)
#elif defined(_CUDA)
        !@cuf istat=cudaStreamSynchronize(nvtx_gpu_stream)
#endif
      end if
    end if
#if defined(_USE_NVTX)
    if(is_nvtx) then
      if(     present(nvtx_color)) then
        call nvtxStartRange(trim(timer_name),color=nvtx_color)
      else if(nvtx_id > 0        ) then
          call nvtxStartRange(trim(timer_name),id=nvtx_id)
      else
        call nvtxStartRange(trim(timer_name))
      end if
    end if
#endif
#if !defined(_TIMER)
    return
#endif
    if(.not.allocated(timer_names)) then
      allocate(timer_names(      0), &
               timer_counts(     0), &
               timer_tictoc(     0), &
               timer_elapsed_acc(0), &
               timer_elapsed_min(0), &
               timer_elapsed_max(0))
#if defined(_OPENACC)
      allocate(timer_tic_slot(0))
#endif
    end if
    idx = timer_search(timer_name)
    if(idx <= 0) then
      ntimers = ntimers + 1
      call concatenate_c(timer_names,timer_name)
      timer_counts      = [timer_counts     ,0          ]
      timer_tictoc      = [timer_tictoc     ,0._dp      ]
      timer_elapsed_acc = [timer_elapsed_acc,0._dp      ]
      timer_elapsed_min = [timer_elapsed_min,huge(0._dp)]
      timer_elapsed_max = [timer_elapsed_max,tiny(0._dp)]
#if defined(_OPENACC)
      timer_tic_slot    = [timer_tic_slot   ,0          ]
#endif
      idx = ntimers
    end if
#if defined(_OPENACC)
    if(present(cuda_stream)) then
      stream_to_use = cuda_stream
    else
      stream_to_use = int(0, acc_handle_kind)
    end if
    !@cuf if(.not.ce_pool_initialized) call timer_init_cuda_events()
    if(ce_slot_counter >= TIMER_MAX_PENDING) then
      error stop 'timer_tic: pending CUDA event pool exhausted; call timer_drain() more often'
    end if
    ce_slot_counter = ce_slot_counter + 1
    timer_tic_slot(idx) = ce_slot_counter
    !@cuf istat = cudaEventRecord(ce_pool_start(ce_slot_counter),stream_to_use)
#else
    timer_tictoc(idx) = MPI_WTIME()
#endif
  end subroutine timer_tic
  subroutine timer_toc(timer_name,nvtx_gpu_stream,ierror,cuda_stream)
    !@cuf use cudafor
    character(*), intent(in) :: timer_name
    integer, intent(in), optional :: nvtx_gpu_stream
    integer, intent(out), optional :: ierror
#if defined(_OPENACC)
    integer(acc_handle_kind), intent(in), optional :: cuda_stream
    integer(acc_handle_kind) :: stream_to_use
#else
    integer, intent(in), optional :: cuda_stream ! unused argument to avoid CPP conditionals in the subroutine definition
#endif
    integer :: idx
    logical :: is_gpu_sync
    !@cuf integer :: istat
    !@cuf integer :: slot_ce
    is_gpu_sync = GPU_DEFAULT_SYNC
    if(present(nvtx_gpu_stream)) then
      if(nvtx_gpu_stream < 0) then
        is_gpu_sync = .false.
      end if
    end if
    if(is_gpu_sync) then
      if(.not.present(nvtx_gpu_stream)) then
#if   defined(_OPENACC)
        !$acc wait
#elif defined(_CUDA)
        !@cuf istat=cudaDeviceSynchronize()
#endif
      else
#if   defined(_OPENACC)
        !$acc wait(nvtx_gpu_stream)
#elif defined(_CUDA)
        !@cuf istat=cudaStreamSynchronize(nvtx_gpu_stream)
#endif
      end if
    end if
#if defined(_USE_NVTX)
    call nvtxEndRange
#endif
#if !defined(_TIMER)
    return
#endif
    if(present(ierror)) ierror = 0
    idx = timer_search(timer_name)
    if(idx > 0) then
#if defined(_OPENACC)
      if(present(cuda_stream)) then
        stream_to_use = cuda_stream
      else
        stream_to_use = int(0, acc_handle_kind)
      end if
      if(ce_drain_n >= TIMER_MAX_PENDING) then
        error stop 'timer_toc: pending CUDA event queue exhausted; call timer_drain() more often'
      end if
      !@cuf slot_ce = timer_tic_slot(idx)
      !@cuf istat = cudaEventRecord(ce_pool_end(slot_ce),stream_to_use)
      !@cuf ce_drain_n = ce_drain_n + 1
      !@cuf ce_drain_slot(ce_drain_n) = slot_ce
      !@cuf ce_drain_timer_idx(ce_drain_n) = idx
#else
      timer_tictoc(idx)      = MPI_WTIME() - timer_tictoc(idx)
      timer_elapsed_acc(idx) =    (timer_elapsed_acc(idx)+timer_tictoc(idx))
      timer_elapsed_min(idx) = min(timer_elapsed_min(idx),timer_tictoc(idx))
      timer_elapsed_max(idx) = max(timer_elapsed_max(idx),timer_tictoc(idx))
      timer_counts(idx)      = timer_counts(idx) + 1
#endif
    else
      if(present(ierror)) ierror = 1
    end if
  end subroutine timer_toc
  subroutine timer_drain()
    !
    ! Resolve all pending CUDA event pairs recorded since the last drain.
    ! Must be called from the host after asynchronous work has been queued.
    !
#if defined(_OPENACC) && defined(_TIMER)
    !@cuf use cudafor
    !@cuf real :: elapsed_ms
    real(dp) :: elapsed
    integer :: i,slot,tidx
    !@cuf integer :: istat
    if(ce_drain_n <= 0) return
    !$acc wait
    do i = 1,ce_drain_n
      slot = ce_drain_slot(i)
      tidx = ce_drain_timer_idx(i)
      elapsed = 0._dp
      !@cuf istat = cudaEventElapsedTime(elapsed_ms,ce_pool_start(slot),ce_pool_end(slot))
      !@cuf elapsed = real(elapsed_ms,dp)*1.e-3_dp
      timer_elapsed_acc(tidx) = timer_elapsed_acc(tidx) + elapsed
      timer_elapsed_min(tidx) = min(timer_elapsed_min(tidx),elapsed)
      timer_elapsed_max(tidx) = max(timer_elapsed_max(tidx),elapsed)
      timer_counts(tidx) = timer_counts(tidx) + 1
    end do
    ce_drain_n = 0
    ce_slot_counter = 0
#endif
  end subroutine timer_drain
  subroutine timer_cleanup
#if !defined(_TIMER)
    return
#endif
#if defined(_OPENACC)
    call timer_drain()
#endif
    if(allocated(timer_names)) then
      deallocate(timer_names,timer_counts,timer_tictoc,timer_elapsed_acc,timer_elapsed_min,timer_elapsed_max)
    end if
#if defined(_OPENACC)
    if(allocated(timer_tic_slot)) deallocate(timer_tic_slot)
    !@cuf call timer_destroy_cuda_events()
#endif
  end subroutine timer_cleanup
  integer function timer_search(timer_name)
    character(*), intent(in) :: timer_name
    integer :: i
    timer_search = -1
    do i = 1,ntimers
      if(timer_names(i) == timer_name) then
        timer_search = i
      end if
    end do
  end function timer_search
  real(dp) function timer_time(timer_name,ierror)
    character(*), intent(in) :: timer_name
    integer, intent(out), optional :: ierror
    integer :: idx
    if(present(ierror)) ierror = 0
    timer_time = -1._dp
    idx = timer_search(timer_name)
    if(idx > 0) then
      timer_time = timer_elapsed_acc(idx)
    else
      if(present(ierror)) ierror = 1
    end if
  end function timer_time
  subroutine concatenate_c(arr,val)
    character(*), intent(inout), allocatable, dimension(:) :: arr
    character(*), intent(in   ) :: val
    character(:), allocatable, dimension(:) :: arr_tmp
    integer :: n
    n = size(arr)
    allocate(arr_tmp,source=arr)
    deallocate(arr); allocate(arr(n+1))
    arr(1:n) = arr_tmp(:); arr(n+1) = val
  end subroutine concatenate_c
  subroutine timer_init_cuda_events()
    !@cuf use cudafor
    !@cuf integer :: i,istat
    !@cuf allocate(ce_pool_start(TIMER_MAX_PENDING),ce_pool_end(TIMER_MAX_PENDING))
    !@cuf do i = 1,TIMER_MAX_PENDING
    !@cuf   istat = cudaEventCreate(ce_pool_start(i))
    !@cuf   istat = cudaEventCreate(ce_pool_end(i))
    !@cuf end do
    !@cuf ce_pool_initialized = .true.
  end subroutine timer_init_cuda_events
  subroutine timer_destroy_cuda_events()
    !@cuf use cudafor
    !@cuf integer :: i,istat
    !@cuf if(.not.ce_pool_initialized) return
    !@cuf do i = 1,TIMER_MAX_PENDING
    !@cuf   istat = cudaEventDestroy(ce_pool_start(i))
    !@cuf   istat = cudaEventDestroy(ce_pool_end(i))
    !@cuf end do
    !@cuf deallocate(ce_pool_start,ce_pool_end)
    !@cuf ce_pool_initialized = .false.
  end subroutine timer_destroy_cuda_events
end module mod_timer
