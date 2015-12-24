MODULE physical_type_methods

  use parallel_module
  use rgrid_variables, only: dV

  implicit none

  PRIVATE
  PUBLIC :: dSpatialIntegral

CONTAINS


  SUBROUTINE dSpatialIntegral( d )
    implicit none
    real(8) :: d(:)
    integer :: i
    call MPI_ALLREDUCE(MPI_IN_PLACE,d,size(d),MPI_REAL8,MPI_SUM,comm_grid,i)
    d=d*dV
  END SUBROUTINE dSpatialIntegral


END MODULE physical_type_methods
