MODULE mpi_wrapper_module

  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: dallgatherv_mpi_wrapper, zallgatherv_mpi_wrapper

CONTAINS


  SUBROUTINE dallgatherv_mpi_wrapper( f, ftot )
    implicit none
    real(8),intent(IN) :: f(:)
    real(8),allocatable,intent(INOUT) :: ftot(:)
    integer :: n,ierr
    if ( .not.allocated(ftot) ) then
       n=sum(ir_grid)
       allocate( ftot(n) ) ; ftot=0.0d0
    end if
    call MPI_ALLGATHERV( f, size(f), MPI_REAL8, ftot, ir_grid, id_grid &
         ,MPI_REAL8, comm_grid, ierr)
  END SUBROUTINE dallgatherv_mpi_wrapper

  SUBROUTINE zallgatherv_mpi_wrapper( f, ftot )
    implicit none
    complex(8),intent(IN) :: f(:)
    complex(8),allocatable,intent(INOUT) :: ftot(:)
    integer :: n,ierr
    if ( .not.allocated(ftot) ) then
       n=sum(ir_grid)
       allocate( ftot(n) ) ; ftot=(0.0d0,0.0d0)
    end if
    call MPI_ALLGATHERV( f, size(f), MPI_COMPLEX16, ftot, ir_grid, id_grid &
         ,MPI_COMPLEX16, comm_grid, ierr)
  END SUBROUTINE zallgatherv_mpi_wrapper


END MODULE mpi_wrapper_module
