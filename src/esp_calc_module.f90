MODULE esp_calc_module

  use hamiltonian_module
  use rgrid_module, only: dV
  use parallel_module
  use wf_module, only: hunk, iflag_hunk

  implicit none

  PRIVATE
  PUBLIC :: esp_calc

CONTAINS

  SUBROUTINE esp_calc(k,s,wf,n1,n2,ns,ne,e)
    implicit none
    integer,intent(IN) :: k,s,n1,n2,ns,ne
#ifdef _DRSDFT_
    real(8),intent(IN) :: wf(n1:n2,ns:ne)
    real(8),allocatable :: hwf(:,:)
#else
    complex(8),intent(IN) :: wf(n1:n2,ns:ne)
    complex(8),allocatable :: hwf(:,:)
#endif
    real(8),intent(OUT) :: e(ns:ne)
    integer :: n,ierr

    e(:)=0.0d0

    allocate( hwf(n1:n2,1) ) ; hwf=0.0d0

    do n=ns,ne
       if ( iflag_hunk >= 1 ) then
          hwf(:,1)=hunk(:,n,k,s)
       else
          call hamiltonian(k,s,wf(n1,n),hwf,n1,n2,n,n)
       end if
#ifdef _DRSDFT_
       e(n) = sum( wf(:,n)*hwf(:,1) )*dV
#else
       e(n) = sum( conjg(wf(:,n))*hwf(:,1) )*dV
#endif
    end do

    deallocate( hwf )

    call mpi_allreduce(MPI_IN_PLACE,e,ne-ns+1,MPI_REAL8,MPI_SUM,comm_grid,ierr)

  END SUBROUTINE esp_calc

END MODULE esp_calc_module
