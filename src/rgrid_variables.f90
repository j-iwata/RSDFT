MODULE rgrid_variables

  implicit none

  PRIVATE
  PUBLIC :: Ngrid, Hgrid, dV, zdV, Igrid

  integer :: Ngrid(0:3),Igrid(2,0:3)
  real(8) :: Hgrid(3)
  real(8) :: dV
#ifdef _DRSDFT_
  real(8) :: zdV
#else
  complex(8) :: zdV
#endif

END MODULE rgrid_variables
