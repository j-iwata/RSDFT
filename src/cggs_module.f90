MODULE cggs_module

  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: cggs

CONTAINS


  SUBROUTINE cggs( iswitch, ms, ns, nt, dv, x, g )
    integer,intent(IN) :: iswitch,ms,ns,nt
    real(8),intent(IN) :: dv
#ifdef _DRSDFT_
    real(8),intent(IN) :: x(ms,ns)
    real(8),intent(INOUT) :: g(ms)
#else
    complex(8),intent(IN) :: x(ms,ns)
    complex(8),intent(INOUT) :: g(ms)
#endif
    select case( iswitch )
    case default
       return
    case( 1 )
       call cggs_1( ms, nt-1, dv, x, g )
    case( 2 )
       call cggs_2( ms, ns, nt, dv, x, g )
    end select
  END SUBROUTINE cggs


  SUBROUTINE cggs_1( ms, ns, dv, x, g )
    implicit none
    integer,intent(IN) :: ms,ns
    real(8),intent(IN) :: dv
#ifdef _DRSDFT_
    real(8),intent(IN) :: x(ms,ns)
    real(8),intent(INOUT) :: g(ms)
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: xg(:)
#else
    complex(8),intent(IN) :: x(ms,ns)
    complex(8),intent(INOUT) :: g(ms)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: xg(:)
#endif
    integer :: i,n

    if ( ns == 0 ) return

    allocate( xg(ns) ) ; xg=zero

    do n=1,ns
#ifdef _DRSDFT_
       do i=1,ms
          xg(n) = xg(n) + x(i,n)*g(i) 
       end do
#else
       do i=1,ms
          xg(n) = xg(n) + conjg( x(i,n) )*g(i) 
       end do
#endif
    end do
    xg(:)=xg(:)*dv
    call mpi_allreduce( MPI_IN_PLACE, xg, ns, MPI_REAL8,MPI_SUM,comm_grid,i )

    do n=1,ns
       do i=1,ms
          g(i) = g(i) - x(i,n)*xg(n)
       end do
    end do

    deallocate( xg )

  END SUBROUTINE cggs_1


  SUBROUTINE cggs_2( ms, ns, nt, dv, x, g )
    implicit none
    integer,intent(IN) :: ms,ns,nt
    real(8),intent(IN) :: dv
#ifdef _DRSDFT_
    real(8),intent(IN) :: x(ms,ns)
    real(8),intent(INOUT) :: g(ms)
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: xg(:)
#else
    complex(8),intent(IN) :: x(ms,ns)
    complex(8),intent(INOUT) :: g(ms)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: xg(:)
#endif
    integer :: i,n

    allocate( xg(ns) ) ; xg=zero

    do n=1,ns
       if ( n == nt ) cycle
#ifdef _DRSDFT_
       do i=1,ms
          xg(n) = xg(n) + x(i,n)*g(i) 
       end do
#else
       do i=1,ms
          xg(n) = xg(n) + conjg( x(i,n) )*g(i) 
       end do
#endif
       xg(n)=xg(n)*dv
    end do
    call mpi_allreduce( MPI_IN_PLACE, xg, ns, MPI_REAL8,MPI_SUM,comm_grid,i )

    do n=1,ns
       if ( n == nt ) cycle
       do i=1,ms
          g(i) = g(i) - x(i,n)*xg(n)
       end do
    end do

    deallocate( xg )

  END SUBROUTINE cggs_2


END MODULE cggs_module
