MODULE eion_module

  use ewald_module
  use eion_mol_module

  implicit none

  PRIVATE
  PUBLIC :: Eewald, calc_eion, init_eion

  real(8) :: Eewald

  integer :: SYStype = 0
  logical :: iswitch_test_ewald = .false.

CONTAINS


  SUBROUTINE init_eion( SYStype_in, disp_switch )
    implicit none
    integer,intent(IN) :: SYStype_in
    logical,intent(IN) :: disp_switch
    call write_border( 80, " init_eion(start)" )
    SYStype = SYStype_in
    if ( SYStype == 0 .and. iswitch_test_ewald ) call test_ewald(Eewald)
    call calc_eion
    if ( disp_switch ) then
       write(*,'(1x,"Eewald(SYSTYPE=",i1,")=",f25.15)'),SYStype,Eewald
    end if
    call write_border( 80, " init_eion(end)" )
  END SUBROUTINE init_eion


  SUBROUTINE calc_eion
    implicit none

    select case(SYStype)
    case default

       call calc_ewald( Eewald )

    case( 1 )

       call calc_eion_mol( Eewald )

    end select

  END SUBROUTINE calc_eion


END MODULE eion_module
