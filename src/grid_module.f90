MODULE grid_module

  use rgrid_variables
  use parallel_module
  use basic_type_factory

  implicit none

  PRIVATE
  PUBLIC :: grid, get_range_rgrid
  PUBLIC :: construct_map_3d_to_1d_grid
  PUBLIC :: construct_map_1d_to_3d_grid
  PUBLIC :: get_map_3d_to_1d_grid
  PUBLIC :: mpi_allgatherv_grid, zmpi_allgatherv_grid
  PUBLIC :: inner_product_grid

  type grid
     type( ArrayRange1D ) :: g1
     type( ArrayRange3D ) :: g3
     real(8) :: spacing(3)
     real(8) :: VolumeElement
     integer :: comm
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
    rgrid%comm = comm_grid
  END SUBROUTINE get_range_rgrid


  SUBROUTINE construct_map_3d_to_1d_grid( Ngrid, Igrid, comm, LLL )
    implicit none
    integer,intent(IN) :: Ngrid(0:3), Igrid(2,0:3), comm
    integer,allocatable,intent(INOUT) :: LLL(:,:,:)
    integer :: i,i1,i2,i3
    if ( .not.allocated(LLL) ) then
       allocate( LLL(0:Ngrid(1)-1,0:Ngrid(2)-1,0:Ngrid(3)-1) )
    end if
    LLL=0
    i=Igrid(1,0)-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       LLL(i1,i2,i3)=i
    end do
    end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,LLL,size(LLL),MPI_INTEGER,MPI_SUM,comm,i)
  END SUBROUTINE construct_map_3d_to_1d_grid


  SUBROUTINE get_map_3d_to_1d_grid( rgrid, LLL )
    implicit none
    type(grid),intent(IN) :: rgrid
    integer,allocatable,intent(INOUT) :: LLL(:,:,:)
    integer :: i,i1,i2,i3,m1,m2,m3,n1,n2,n3
    if ( .not.allocated(LLL) ) then
       m1=rgrid%g3%x%head_global
       n1=rgrid%g3%x%tail_global
       m2=rgrid%g3%y%head_global
       n2=rgrid%g3%y%tail_global
       m3=rgrid%g3%z%head_global
       n3=rgrid%g3%z%tail_global
       allocate( LLL(m1:n1,m2:n2,m3:n3) )
    end if
    LLL=0
    i=rgrid%g1%head-1
    do i3=rgrid%g3%z%head,rgrid%g3%z%tail
    do i2=rgrid%g3%y%head,rgrid%g3%y%tail
    do i1=rgrid%g3%x%head,rgrid%g3%x%tail
       i=i+1
       LLL(i1,i2,i3)=i
    end do
    end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,LLL,size(LLL),MPI_INTEGER,MPI_SUM,rgrid%comm,i)
  END SUBROUTINE get_map_3d_to_1d_grid


  SUBROUTINE construct_map_1d_to_3d_grid( Ngrid, Igrid, comm, LL )
    implicit none
    integer,intent(IN) :: Ngrid(0:3), Igrid(2,0:3), comm
    integer,allocatable,intent(INOUT) :: LL(:,:)
    integer :: i,i1,i2,i3
    if ( .not.allocated(LL) ) then
       allocate( LL(3,Ngrid(0)) )
    end if
    LL=0
    i=Igrid(1,0)-1
    do i3=Igrid(1,3),Igrid(2,3)
    do i2=Igrid(1,2),Igrid(2,2)
    do i1=Igrid(1,1),Igrid(2,1)
       i=i+1
       LL(1,i)=i1
       LL(2,i)=i2
       LL(3,i)=i3
    end do
    end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,LL,size(LL),MPI_INTEGER,MPI_SUM,comm,i)
  END SUBROUTINE construct_map_1d_to_3d_grid


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
