MODULE grid_module

  use rgrid_variables
  use parallel_module
  use basic_type_factory

  implicit none

  PRIVATE
  PUBLIC :: grid, get_range_rgrid
  PUBLIC :: get_map_3d_to_1d
  PUBLIC :: mpi_allgatherv_grid, zmpi_allgatherv_grid
  PUBLIC :: inner_product_grid

  type grid
     type( ArrayRange1D ) :: g1
     type( ArrayRange3D ) :: g3
     real(8) :: spacing(3)
     real(8) :: VolumeElement
  end type grid

CONTAINS

  SUBROUTINE get_range_rgrid( rgrid )
    implicit none
    type(grid) :: rgrid
    rgrid%g1%head = Igrid(1,0)
    rgrid%g1%tail = Igrid(2,0)
    rgrid%g1%size = Igrid(2,0)-Igrid(1,0)+1
    rgrid%g1%head_global = 1
    rgrid%g1%tail_global = Ngrid(0)
    rgrid%g1%size_global = Ngrid(0)
    rgrid%g3%x%head = Igrid(1,1)
    rgrid%g3%x%tail = Igrid(2,1)
    rgrid%g3%x%size = Igrid(2,1)-Igrid(1,1)+1
    rgrid%g3%y%head = Igrid(1,2)
    rgrid%g3%y%tail = Igrid(2,2)
    rgrid%g3%y%size = Igrid(2,2)-Igrid(1,2)+1
    rgrid%g3%z%head = Igrid(1,3)
    rgrid%g3%z%tail = Igrid(2,3)
    rgrid%g3%z%size = Igrid(2,3)-Igrid(1,3)+1
    rgrid%g3%x%head_global = 0
    rgrid%g3%x%tail_global = Ngrid(1)-1
    rgrid%g3%x%size_global = Ngrid(1)
    rgrid%g3%y%head_global = 0
    rgrid%g3%y%tail_global = Ngrid(2)-1
    rgrid%g3%y%size_global = Ngrid(2)
    rgrid%g3%z%head_global = 0
    rgrid%g3%z%tail_global = Ngrid(3)-1
    rgrid%g3%z%size_global = Ngrid(3)
    rgrid%spacing(1:3) = Hgrid(1:3)
    rgrid%VolumeElement = dV
  END SUBROUTINE get_range_rgrid

  SUBROUTINE get_map_3d_to_1d( LLL )
     implicit none
     integer,allocatable,intent(OUT) :: LLL(:,:,:)
     integer :: i,i1,i2,i3,ierr
     allocate( LLL(0:Ngrid(1)-1,0:Ngrid(2)-1,0:Ngrid(3)-1) ) ; LLL=0
     i=Igrid(1,0)-1
     do i3=Igrid(1,3),Igrid(2,3)
     do i2=Igrid(1,2),Igrid(2,2)
     do i1=Igrid(1,1),Igrid(2,1)
        i=i+1
        LLL(i1,i2,i3)=i
     end do
     end do
     end do
     call MPI_ALLREDUCE( MPI_IN_PLACE, LLL, size(LLL), MPI_INTEGER, MPI_SUM, comm_grid, ierr )
  END SUBROUTINE get_map_3d_to_1d

  SUBROUTINE mpi_allgatherv_grid( f, g )
    implicit none
    real(8),intent(IN)  :: f(:)
    real(8),intent(OUT) :: g(:)
    integer :: ierr
    call MPI_ALLGATHERV( f, size(f), MPI_REAL8, g, ir_grid, id_grid, &
                         MPI_REAL8, comm_grid, ierr )
  END SUBROUTINE mpi_allgatherv_grid

  SUBROUTINE zmpi_allgatherv_grid( zf, zg )
    implicit none
    complex(8),intent(IN)  :: zf(:)
    complex(8),intent(OUT) :: zg(:)
    integer :: ierr
    call MPI_ALLGATHERV( zf, size(zf), MPI_COMPLEX16, zg, ir_grid, &
                         id_grid, MPI_COMPLEX16, comm_grid, ierr )
  END SUBROUTINE zmpi_allgatherv_grid


  SUBROUTINE inner_product_grid( f, g, c, fgc )
    implicit none
    real(8),intent(IN) :: f(:), g(:), c
    real(8),intent(OUT) :: fgc
    real(8) :: fgc_tmp
    integer :: i
    fgc_tmp=sum(f*g)*c
    call MPI_ALLREDUCE( fgc_tmp, fgc, 1, MPI_REAL8, MPI_SUM, comm_grid, i )
  END SUBROUTINE inner_product_grid


END MODULE grid_module
