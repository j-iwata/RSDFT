MODULE xc_variables

  use basic_type_factory, only: GSArray

  implicit none

  PRIVATE
  PUBLIC :: xcpot, xcene

  type xcpot
     type( GSArray ) :: xc
     type( GSArray ) :: x
     type( GSArray ) :: c
  end type xcpot

  type xcene
     real(8) :: Exc
     real(8) :: Ex
     real(8) :: Ec
  end type xcene

END MODULE xc_variables
