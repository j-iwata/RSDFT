MODULE hartree_variables

  implicit none

  PRIVATE
  PUBLIC :: E_hartree, Vh

  real(8) :: E_hartree
  real(8),allocatable :: Vh(:)

END MODULE hartree_variables
