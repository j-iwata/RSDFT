MODULE rsdft_mpi_module

  implicit none

  PRIVATE
  PUBLIC :: rsdft_allreduce_sum
  PUBLIC :: d_rsdft_allreduce_sum_5, d_rsdft_allreduce_sum_6
  PUBLIC :: rsdft_bcast
  PUBLIC :: rsdft_allgather
  PUBLIC :: rsdft_allgatherv

  INTERFACE rsdft_allreduce_sum
     MODULE PROCEDURE d_rsdft_allreduce_sum_0, z_rsdft_allreduce_sum_0, &
                      i_rsdft_allreduce_sum_1,                          &
                      d_rsdft_allreduce_sum_1, z_rsdft_allreduce_sum_1, &
                      d_rsdft_allreduce_sum_2, z_rsdft_allreduce_sum_2, &
                      i_rsdft_allreduce_sum_3,                          &
                      d_rsdft_allreduce_sum_3, z_rsdft_allreduce_sum_3
  END INTERFACE

  INTERFACE rsdft_bcast
     MODULE PROCEDURE d_rsdft_bcast2, z_rsdft_bcast2
  END INTERFACE

  INTERFACE rsdft_allgather
     MODULE PROCEDURE l_rsdft_allgather12, &
                      d_rsdft_allgather12, &
                      d_rsdft_allgather23, &
                      d_rsdft_allgather34, &
                      d_rsdft_allgather11, &
                      d_rsdft_allgather22, &
                      z_rsdft_allgather23
  END INTERFACE

  INTERFACE rsdft_allgatherv
     MODULE PROCEDURE i_rsdft_allgatherv22, &
                      d_rsdft_allgatherv11, &
                      d_rsdft_allgatherv22, &
                      d_rsdft_allgatherv33, &
                      z_rsdft_allgatherv22, &
                      z_rsdft_allgatherv33
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


  SUBROUTINE i_rsdft_allreduce_sum_1( a, comm )
    implicit none
    integer,intent(INOUT) :: a(:)
    integer,intent(IN) :: comm
    integer :: n, ierr
    integer,allocatable :: a0(:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    n=size(a)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, n, MPI_INTEGER, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,n,MPI_INTEGER,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE i_rsdft_allreduce_sum_1

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


  SUBROUTINE i_rsdft_allreduce_sum_3( a, comm )
    implicit none
    integer,intent(INOUT) :: a(:,:,:)
    integer,intent(IN) :: comm
    integer :: l, m, n, ierr
    integer,allocatable :: a0(:,:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    l=size(a,1)
    m=size(a,2)
    n=size(a,3)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(l,m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, l*m*n, MPI_INTEGER, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,l*m*n,MPI_INTEGER,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE i_rsdft_allreduce_sum_3

  SUBROUTINE d_rsdft_allreduce_sum_3( a, comm )
    implicit none
    real(8),intent(INOUT) :: a(:,:,:)
    integer,intent(IN) :: comm
    integer :: l, m, n, ierr
    real(8),allocatable :: a0(:,:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    l=size(a,1)
    m=size(a,2)
    n=size(a,3)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(l,m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, l*m*n, MPI_REAL8, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,l*m*n,MPI_REAL8,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_3

  SUBROUTINE z_rsdft_allreduce_sum_3( a, comm )
    implicit none
    complex(8),intent(INOUT) :: a(:,:,:)
    integer,intent(IN) :: comm
    integer :: l, m, n, ierr
    complex(8),allocatable :: a0(:,:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    l=size(a,1)
    m=size(a,2)
    n=size(a,3)
#ifdef _NO_MPI_INPLACE_
    allocate( a0(l,m,n) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, l*m*n, MPI_COMPLEX16, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,l*m*n,MPI_COMPLEX16,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE z_rsdft_allreduce_sum_3


  SUBROUTINE d_rsdft_allreduce_sum_5( a, comm )
    implicit none
    real(8),intent(INOUT) :: a(:,:,:,:,:)
    integer,intent(IN) :: comm
    integer :: i, m(6), mm, ierr
    real(8),allocatable :: a0(:,:,:,:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    mm=1
    do i=1,5
       m(i)=size(a,i)
       mm=mm*m(i)
    end do
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m(1),m(2),m(3),m(4),m(5)) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, mm, MPI_REAL8, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,mm,MPI_REAL8,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_5


  SUBROUTINE d_rsdft_allreduce_sum_6( a, comm )
    implicit none
    real(8),intent(INOUT) :: a(:,:,:,:,:,:)
    integer,intent(IN) :: comm
    integer :: i, m(6), mm, ierr
    real(8),allocatable :: a0(:,:,:,:,:,:)
#ifdef _NOMPI_
    return
#else
    include 'mpif.h'
    mm=1
    do i=1,6
       m(i)=size(a,i)
       mm=mm*m(i)
    end do
#ifdef _NO_MPI_INPLACE_
    allocate( a0(m(1),m(2),m(3),m(4),m(5),m(6)) ) ; a0=a
    call MPI_ALLREDUCE( a0, a, mm, MPI_REAL8, MPI_SUM, comm, ierr )
    deallocate( a0 )
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,a,mm,MPI_REAL8,MPI_SUM,comm,ierr)
#endif
#endif
  END SUBROUTINE d_rsdft_allreduce_sum_6

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

!-------------------------------------------------- rsdft_allgather

  SUBROUTINE l_rsdft_allgather12( f, g, comm )
    implicit none
    logical,intent(INOUT) :: f(:), g(:,0:)
    integer,intent(IN) :: comm
    integer :: i
    logical,allocatable :: t(:)
#ifdef _NO_MPI_
    g(:,0) = f(:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f)) ) ; t=f
    call MPI_ALLGATHER( t,size(f),MPI_LOGICAL, g,size(f),MPI_LOGICAL, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f,size(f),MPI_LOGICAL, g,size(f),MPI_LOGICAL, comm, i )
#endif
#endif
  END SUBROUTINE l_rsdft_allgather12

  SUBROUTINE d_rsdft_allgather12( f, g, comm )
    implicit none
    real(8),intent(INOUT) :: f(:), g(:,0:)
    integer,intent(IN) :: comm
    integer :: i
    real(8),allocatable :: t(:)
#ifdef _NO_MPI_
    g(:,0) = f(:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f)) ) ; t=f
    call MPI_ALLGATHER( t, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgather12

  SUBROUTINE d_rsdft_allgather23( f, g, comm )
    implicit none
    real(8),intent(INOUT) :: f(:,:), g(:,:,0:)
    integer,intent(IN) :: comm
    integer :: i
    real(8),allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:,0) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHER( t, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgather23

  SUBROUTINE d_rsdft_allgather34( f, g, comm )
    implicit none
    real(8),intent(INOUT) :: f(:,:,:), g(:,:,:,0:)
    integer,intent(IN) :: comm
    integer :: i
    real(8),allocatable :: t(:,:,:)
#ifdef _NO_MPI_
    g(:,:,:,0) = f(:,:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2),size(f,3)) ) ; t=f
    call MPI_ALLGATHER( t, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgather34

  SUBROUTINE d_rsdft_allgather11( f, g, comm )
    implicit none
    real(8),intent(INOUT) :: f(:), g(:)
    integer,intent(IN) :: comm
    integer :: i
    real(8),allocatable :: t(:)
#ifdef _NO_MPI_
    g(:) = f(:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1)) ) ; t=f
    call MPI_ALLGATHER( t, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgather11

  SUBROUTINE d_rsdft_allgather22( f, g, comm )
    implicit none
    real(8),intent(INOUT) :: f(:,:), g(:,:)
    integer,intent(IN) :: comm
    integer :: i
    real(8),allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHER( t, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHER( f, size(f), MPI_REAL8, g, size(f), MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgather22

  SUBROUTINE z_rsdft_allgather23( f, g, comm )
    implicit none
    complex(8),intent(INOUT) :: f(:,:), g(:,:,0:)
    integer,intent(IN) :: comm
    integer :: i
    complex(8),allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:,0) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHER(t,size(f),MPI_COMPLEX16,g,size(f),MPI_COMPLEX16,comm,i)
    deallocate( t )
#else
    call MPI_ALLGATHER(f,size(f),MPI_COMPLEX16,g,size(f),MPI_COMPLEX16,comm,i)
#endif
#endif
  END SUBROUTINE z_rsdft_allgather23

!-------------------------------------------------- rsdft_allgatherv

  SUBROUTINE i_rsdft_allgatherv22( f, g, ir, id, comm )
    implicit none
    integer,intent(INOUT) :: f(:,:), g(:,:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    integer,allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHERV( t,size(f),MPI_INTEGER, g,ir,id,MPI_INTEGER, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHERV( f,size(f),MPI_INTEGER, g,ir,id,MPI_INTEGER, comm, i )
#endif
#endif
  END SUBROUTINE i_rsdft_allgatherv22

  SUBROUTINE d_rsdft_allgatherv11( f, g, ir, id, comm )
    implicit none
    real(8),intent(INOUT) :: f(:), g(:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    real(8),allocatable :: t(:)
#ifdef _NO_MPI_
    g(:) = f(:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1)) ) ; t=f
    call MPI_ALLGATHERV( t, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHERV( f, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgatherv11

  SUBROUTINE d_rsdft_allgatherv22( f, g, ir, id, comm )
    implicit none
    real(8),intent(INOUT) :: f(:,:), g(:,:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    real(8),allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHERV( t, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHERV( f, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgatherv22

  SUBROUTINE d_rsdft_allgatherv33( f, g, ir, id, comm )
    implicit none
    real(8),intent(INOUT) :: f(:,:,:), g(:,:,:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    real(8),allocatable :: t(:,:,:)
#ifdef _NO_MPI_
    g(:,:,:) = f(:,:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2),size(f,3)) ) ; t=f
    call MPI_ALLGATHERV( t, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
    deallocate( t )
#else
    call MPI_ALLGATHERV( f, size(f), MPI_REAL8, g, ir, id, MPI_REAL8, comm, i )
#endif
#endif
  END SUBROUTINE d_rsdft_allgatherv33

  SUBROUTINE z_rsdft_allgatherv22( f, g, ir, id, comm )
    implicit none
    complex(8),intent(INOUT) :: f(:,:), g(:,:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    complex(8),allocatable :: t(:,:)
#ifdef _NO_MPI_
    g(:,:) = f(:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2)) ) ; t=f
    call MPI_ALLGATHERV(t,size(f),MPI_COMPLEX16,g,ir,id,MPI_COMPLEX16,comm,i)
    deallocate( t )
#else
    call MPI_ALLGATHERV(f,size(f),MPI_COMPLEX16,g,ir,id,MPI_COMPLEX16,comm,i)
#endif
#endif
  END SUBROUTINE z_rsdft_allgatherv22

  SUBROUTINE z_rsdft_allgatherv33( f, g, ir, id, comm )
    implicit none
    complex(8),intent(INOUT) :: f(:,:,:), g(:,:,:)
    integer,intent(IN) :: ir(:),id(:),comm
    integer :: i
    complex(8),allocatable :: t(:,:,:)
#ifdef _NO_MPI_
    g(:,:,:) = f(:,:,:)
    return
#else
    include 'mpif.h'
#ifdef _NO_MPI_INPLACE_
    allocate( t(size(f,1),size(f,2),size(f,3)) ) ; t=f
    call MPI_ALLGATHERV(t,size(f),MPI_COMPLEX16,g,ir,id,MPI_COMPLEX16,comm,i)
    deallocate( t )
#else
    call MPI_ALLGATHERV(f,size(f),MPI_COMPLEX16,g,ir,id,MPI_COMPLEX16,comm,i)
#endif
#endif
  END SUBROUTINE z_rsdft_allgatherv33


END MODULE rsdft_mpi_module
