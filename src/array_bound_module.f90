MODULE array_bound_module

  use parallel_module, only: id_grid,ir_grid,id_band,ir_band &
                            ,id_bzsm,ir_bzsm,id_spin,ir_spin &
                            ,myrank_g,myrank_b,myrank_k,myrank_s
  use basic_type_factory

  implicit none

  PRIVATE
  PUBLIC :: set_array_bound
  PUBLIC :: get_grid_range_local, get_grid_range_globl
  PUBLIC :: get_spin_range_local, get_spin_range_globl

  integer,PUBLIC :: ML ,ML_0 ,ML_1
  integer,PUBLIC :: MB ,MB_0 ,MB_1
  integer,PUBLIC :: MBZ,MBZ_0,MBZ_1
  integer,PUBLIC :: MSP,MSP_0,MSP_1

CONTAINS


  SUBROUTINE set_array_bound
    call write_border( 80, " set_array_bound(start)" )
    ML_0  = id_grid(myrank_g)+1
    ML_1  = id_grid(myrank_g)+ir_grid(myrank_g)
    MB_0  = id_band(myrank_b)+1
    MB_1  = id_band(myrank_b)+ir_band(myrank_b)
    MBZ_0 = id_bzsm(myrank_k)+1
    MBZ_1 = id_bzsm(myrank_k)+ir_bzsm(myrank_k)
    MSP_0 = id_spin(myrank_s)+1
    MSP_1 = id_spin(myrank_s)+ir_spin(myrank_s)
    ML  = sum(ir_grid)
    MB  = sum(ir_band)
    MBZ = sum(ir_bzsm)
    MSP = sum(ir_spin)
    call write_border( 80, " set_array_bound(end)" )
  END SUBROUTINE set_array_bound


  SUBROUTINE get_grid_range_local( g )
    implicit none
    type( ArrayRange ) :: g
    g%head = ML_0
    g%tail = ML_1
    call getArraySize( g )
  END SUBROUTINE get_grid_range_local

  SUBROUTINE get_grid_range_globl( g )
    implicit none
    type( ArrayRange ) :: g
    g%head = 1
    g%tail = ML
    call getArraySize( g )
  END SUBROUTINE get_grid_range_globl


  SUBROUTINE get_spin_range_local( g )
    implicit none
    type( ArrayRange ) :: g
    g%head = MSP_0
    g%tail = MSP_1
    call getArraySize( g )
  END SUBROUTINE get_spin_range_local

  SUBROUTINE get_spin_range_globl( g )
    implicit none
    type( ArrayRange ) :: g
    g%head = 1
    g%tail = MSP
    call getArraySize( g )
  END SUBROUTINE get_spin_range_globl


  SUBROUTINE getArraySize( range )
    implicit none
    type( ArrayRange ) :: range
    range%size = range%tail - range%head + 1
  END SUBROUTINE getArraySize


END MODULE array_bound_module
