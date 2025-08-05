! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module diezdecomp
#if defined(_OPENACC) && defined(_USE_DIEZDECOMP)
  !
  ! named constants
  !
  use diezdecomp_api_cans, only: CUDECOMP_RANK_NULL => diezdecomp_rank_null
  !
  ! types
  !
  use diezdecomp_api_cans, only: cudecompHandle     => diezdecompHandle, &
                                 cudecompGridDesc   => diezdecompGridDesc, &
                                 cudecompPencilInfo => diezdecompPencilInfo
  !
  ! initialization functions
  !
  use diezdecomp_api_cans, only: diezdecompGridDescCreate
  !
  ! query functions - pencils and decomposition
  !
  use diezdecomp_api_cans, only: cudecompGetPencilInfo  => diezdecompGetPencilInfo, &
                                 cudecompGetShiftedRank => diezdecompGetShiftedRank
  !
  ! query functions - workspace sizes
  !
  use diezdecomp_api_cans, only: cudecompGetTransposeWorkspaceSize => diezdecomp_get_workspace_size_transposes, &
                                 cudecompGetHaloWorkspaceSize      => diezdecomp_get_workspace_size_halos
  !
  ! collective transposes functions
  !
  use diezdecomp_api_cans, only: cudecompTransposeXtoY => diezdecompTransposeXtoY, &
                                 cudecompTransposeYtoX => diezdecompTransposeYtoX, &
                                 cudecompTransposeYtoZ => diezdecompTransposeYtoZ, &
                                 cudecompTransposeZtoY => diezdecompTransposeZtoY, &
                                                          diezdecompTransposeXtoZ, &
                                                          diezdecompTransposeZtoX
  !
  ! halo update functions
  !
  use diezdecomp_api_cans, only: cudecompUpdateHalosX => diezdecompUpdateHalosX, &
                                 cudecompUpdateHalosY => diezdecompUpdateHalosY, &
                                 cudecompUpdateHalosZ => diezdecompUpdateHalosZ
  implicit none
  public
#endif
end module diezdecomp
