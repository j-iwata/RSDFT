MODULE ps_nloc_initiate_module

  use ps_nloc2_init_module, only: ps_nloc2_init,rcfac,qcfac,etafac
  use ps_nloc2_module, only: prep_ps_nloc2
  use ps_nloc3_module, only: init_ps_nloc3, prep_ps_nloc3
  use ps_nloc_mr_module, only: prep_ps_nloc_mr
  use pseudopot_module, only: pselect
  use ps_q_init_module, only: initKtoKPSQ_ps_q_init, ps_q_init
  use ps_qrij_prep_module, only: prepQRijp102
  use ps_prepNzqr_g_module, only: prepNzqr
  use var_sys_parameter, only: pp_kind
  use var_ps_member, only: ps_type

  implicit none

  PRIVATE
  PUBLIC :: ps_nloc_initiate

CONTAINS

  SUBROUTINE ps_nloc_initiate( Gcut )

    implicit none
    real(8),intent(IN) :: Gcut

    call write_border( 80, " ps_nloc_initiate(start)" )

    select case( pselect )
    case default

       pp_kind="NCPP"

       call ps_nloc2_init( Gcut )
       if ( ps_type == 0 ) then
          call prep_ps_nloc2
       else if ( ps_type == 1 ) then
          call prep_ps_nloc_mr
       end if

    case( 3 )

       pp_kind="NCPP"

       call init_ps_nloc3
       call prep_ps_nloc3

    case( 102 )

       pp_kind="USPP"

       call initKtoKPSQ_ps_q_init
       call ps_nloc2_init( Gcut )
       call ps_q_init( Gcut,rcfac,qcfac,etafac )
       call prep_ps_nloc2
       call prepNzqr
       call prepQRijp102

    end select

    call write_border( 80, " ps_nloc_initiate(end)" )

  END SUBROUTINE ps_nloc_initiate

END MODULE ps_nloc_initiate_module
