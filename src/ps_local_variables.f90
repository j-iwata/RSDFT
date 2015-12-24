MODULE ps_local_variables

  implicit none

  real(8),allocatable :: vqlg(:,:)

  type pslocal
     real(8),allocatable :: vqlg(:,:)
     real(8) :: const_ps_local
     real(8),allocatable :: Vion(:)
  end type pslocal

  type(pslocal),PUBLIC :: psloc

END MODULE ps_local_variables
