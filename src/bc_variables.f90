MODULE bc_variables

  implicit none

  PRIVATE
  PUBLIC :: n_neighbor, fdinfo_send, fdinfo_recv
  PUBLIC :: www, sbuf, rbuf
  PUBLIC :: TYPE_MAIN, zero, Md

  include 'mpif.h'

  integer :: n_neighbor(6)
  integer,allocatable :: fdinfo_send(:,:,:),fdinfo_recv(:,:,:)

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN=MPI_REAL8
  real(8),allocatable :: www(:,:,:,:)
  real(8),allocatable :: sbuf(:,:,:),rbuf(:,:,:)
  real(8),parameter :: zero=0.0d0
#else
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
  complex(8),allocatable :: www(:,:,:,:)
  complex(8),allocatable :: sbuf(:,:,:),rbuf(:,:,:)
  complex(8),parameter :: zero=(0.d0,0.d0)
#endif

  integer :: Md

END MODULE bc_variables
