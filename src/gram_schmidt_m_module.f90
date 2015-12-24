MODULE gram_schmidt_m_module

  use rgrid_module
  use wf_module
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: gram_schmidt_m

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN=MPI_REAL8
  real(8) :: uu,vv
#else
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
  complex(8) :: uu,vv
#endif

CONTAINS

  SUBROUTINE gram_schmidt_m(n0,n1,k,s)
    implicit none
    integer,intent(IN) :: n0,n1,k,s
    integer :: n,m,ierr
    real(8) :: c,d

    call gather_b_wf( k,s )

    do n=n0,n1
       do m=n0,n-1
#ifdef _DRSDFT_
          uu=sum( unk(:,m,k,s)*unk(:,n,k,s) )*dV
#else
          uu=sum( conjg(unk(:,m,k,s))*unk(:,n,k,s) )*dV
#endif
          call mpi_allreduce(uu,vv,1,TYPE_MAIN,MPI_SUM,comm_grid,ierr)
          unk(:,n,k,s)=unk(:,n,k,s)-unk(:,m,k,s)*vv
       end do
       c=sum(abs(unk(:,n,k,s))**2)*dV
       call mpi_allreduce(c,d,1,MPI_REAL8,MPI_SUM,comm_grid,ierr)
       c=1.d0/sqrt(d)
       unk(:,n,k,s)=c*unk(:,n,k,s)
    end do

  END SUBROUTINE gram_schmidt_m

END MODULE gram_schmidt_m_module
