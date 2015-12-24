MODULE subspace_diag_variables

  implicit none

  PRIVATE
  PUBLIC :: MB_diag, Hsub, Vsub, mat_block &
           ,NBLK1, NBLK2, zero, one, TYPE_MAIN

  include 'mpif.h'

  integer,allocatable :: mat_block(:,:)

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN = MPI_REAL8
  real(8),allocatable :: Hsub(:,:), Vsub(:,:)
  real(8),parameter :: zero=0.d0,one=1.d0
#else
  integer,parameter :: TYPE_MAIN = MPI_COMPLEX16
  complex(8),allocatable :: Hsub(:,:), Vsub(:,:)
  complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
#endif
  integer :: MB_diag,NBLK1,NBLK2

END MODULE subspace_diag_variables
