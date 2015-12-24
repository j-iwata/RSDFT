MODULE ewald_variables

  implicit none

  PRIVATE
  PUBLIC :: eta, mg, mr, LG, LR, ipair, mpair

  real(8) :: eta
  integer :: mg,mr
  integer,allocatable :: LG(:,:),LR(:,:)
  integer :: mpair
  integer,allocatable :: ipair(:,:)

END MODULE ewald_variables
