MODULE rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: rsdft_allreduce_sum
  PUBLIC :: rsdft_bcast

  INTERFACE rsdft_allreduce_sum
     MODULE PROCEDURE d_rsdft_allreduce_sum_0, z_rsdft_allreduce_sum_0, &
                      d_rsdft_allreduce_sum_1, z_rsdft_allreduce_sum_1, &
                      d_rsdft_allreduce_sum_2, z_rsdft_allreduce_sum_2
  END INTERFACE

  INTERFACE rsdft_bcast
     MODULE PROCEDURE d_rsdft_bcast2, z_rsdft_bcast2
  END INTERFACE

CONTAINS

!-------------------------------------------------- rsdft_allreduce_sum

  SUBROUTINE d_rsdft_allreduce_sum_0( a, comm )
    implicit none
    real(8),intent(INOUT) :: a
    integer,intent(IN) :: comm
    integer :: ierr
    real(8) :: a0
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    a0=a
    call MPI_ALLREDUCE( a0, a, 1, MPI_REAL8, MPI_SUM, comm, ierr )
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_0

  SUBROUTINE z_rsdft_allreduce_sum_0( a, comm )
    implicit none
    complex(8),intent(INOUT) :: a
    integer,intent(IN) :: comm
    integer :: ierr
    complex(8) :: a0
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    a0=a
    call MPI_ALLREDUCE( a0, a, 1, MPI_COMPLEX16, MPI_SUM, comm, ierr )
#endif
  END SUBROUTINE z_rsdft_allreduce_sum_0


  SUBROUTINE d_rsdft_allreduce_sum_1( a, comm )
    implicit none
    real(8),intent(INOUT) :: a(:)
    integer,intent(IN) :: comm
    integer :: n, ierr
    real(8),allocatable :: a0(:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    n=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, n, MPI_REAL8, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,n,MPI_REAL8,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_1

  SUBROUTINE z_rsdft_allreduce_sum_1( a, comm )
    implicit none
    complex(8),intent(INOUT) :: a(:)
    integer,intent(IN) :: comm
    integer :: n, ierr
    complex(8),allocatable :: a0(:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    n=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, n, MPI_COMPLEX16, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,n,MPI_COMPLEX16,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE z_rsdft_allreduce_sum_1


  SUBROUTINE d_rsdft_allreduce_sum_2( a, comm )
    implicit none
    real(8),intent(INOUT) :: a(:,:)
    integer,intent(IN) :: comm
    integer :: m, n, ierr
    real(8),allocatable :: a0(:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    m=size(a,1)
    n=size(a,2)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, m*n, MPI_REAL8, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,m*n,MPI_REAL8,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_2

  SUBROUTINE z_rsdft_allreduce_sum_2( a, comm )
    implicit none
    complex(8),intent(INOUT) :: a(:,:)
    integer,intent(IN) :: comm
    integer :: m, n, ierr
    complex(8),allocatable :: a0(:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    m=size(a,1)
    n=size(a,2)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, m*n, MPI_COMPLEX16, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,m*n,MPI_COMPLEX16,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE z_rsdft_allreduce_sum_2

!-------------------------------------------------- rsdft_bcast

  SUBROUTINE d_rsdft_bcast2( a, src, comm )
    implicit none
    real(8),intent(INOUT) :: a(:,:)
    integer,intent(IN) :: src, comm
    integer :: ierr
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    call MPI_BCAST( a, size(a), MPI_REAL8, src, comm, ierr )
#endif
  END SUBROUTINE d_rsdft_bcast2

  SUBROUTINE z_rsdft_bcast2( a, src, comm )
    implicit none
    complex(8),intent(INOUT) :: a(:,:)
    integer,intent(IN) :: src, comm
    integer :: ierr
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    call MPI_BCAST( a, size(a), MPI_COMPLEX16, src, comm, ierr )
#endif
  END SUBROUTINE z_rsdft_bcast2


END MODULE rsdft_mpi_module
