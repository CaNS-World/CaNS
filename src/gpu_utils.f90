! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_gpu_utils
#if 0
  use mod_types
  !@cuf use cudafor, only:cudaMemcpyAsync,cudaMemcpyHostToDevice,cudaMemcpyDeviceToHost
  implicit none
  private
  public to_gpu_cuda,to_cpu_cuda
  interface to_gpu_cuda
    module procedure to_gpu_cuda_sp,to_gpu_cuda_dp
  end interface
  interface to_cpu_cuda
    module procedure to_cpu_cuda_sp,to_cpu_cuda_dp
  end interface
  contains
   subroutine to_gpu_cuda_dp(a_h,a_d,istream)
     implicit none
     real(dp), intent(inout), allocatable, dimension(..) :: a_h
     real(dp), intent(inout), allocatable, dimension(..) :: a_d
     integer , intent(inout), optional :: istream
     !@cuf integer :: istat
     !@cuf attributes(device) :: a_d
#if defined(_CUDA)
     if(.not.allocated(a_d)) allocate(a_d,mold=a_h)
     if(present(istream)) then
       istat = cudaMemcpyAsync(a_d,a_h,size(a_d),cudaMemcpyHostToDevice,istream)
     else
       a_d = a_h
     end if
#else
     call move_alloc(a_h,a_d)
#endif
   end subroutine to_gpu_cuda_dp
   subroutine to_gpu_cuda_sp(a_h,a_d,istream)
     implicit none
     real(sp), intent(inout), allocatable, dimension(..) :: a_h
     real(sp), intent(inout), allocatable, dimension(..) :: a_d
     integer , intent(inout), optional :: istream
     !@cuf integer :: istat
     !@cuf attributes(device) :: a_d
#if defined(_CUDA)
     if(.not.allocated(a_d)) allocate(a_d,mold=a_h)
     if(present(istream)) then
       istat = cudaMemcpyAsync(a_d,a_h,size(a_d),cudaMemcpyHostToDevice,istream)
     else
       a_d = a_h
     end if
#else
     call move_alloc(a_h,a_d)
#endif
   end subroutine to_gpu_cuda_sp
   subroutine to_cpu_cuda_dp(a_d,a_h,istream)
     implicit none
     real(dp), intent(inout), allocatable, dimension(..) :: a_d
     real(dp), intent(inout), allocatable, dimension(..) :: a_h
     integer , intent(inout), optional :: istream
     !@cuf integer :: istat
     !@cuf attributes(device) :: a_d
#if defined(_CUDA)
     if(.not.allocated(a_d)) allocate(a_d,mold=a_h)
     if(present(istream)) then
       istat = cudaMemcpyAsync(a_h,a_d,size(a_h),cudaMemcpyDevicetoHost,istream)
     else
       a_h = a_d
     end if
#else
     call move_alloc(a_d,a_h)
#endif
   end subroutine to_cpu_cuda_dp
   subroutine to_cpu_cuda_sp(a_d,a_h,istream)
     implicit none
     real(sp), intent(inout), allocatable, dimension(..) :: a_d
     real(sp), intent(inout), allocatable, dimension(..) :: a_h
     integer , intent(inout), optional :: istream
     !@cuf integer :: istat
     !@cuf attributes(device) :: a_d
#if defined(_CUDA)
     if(.not.allocated(a_h)) allocate(a_h,mold=a_d)
     if(present(istream)) then
       istat = cudaMemcpyAsync(a_h,a_d,size(a_h),cudaMemcpyDevicetoHost,istream)
     else
       a_h = a_d
     end if
#else
     call move_alloc(a_d,a_h)
#endif
   end subroutine to_cpu_cuda_sp
#endif
end module mod_gpu_utils
