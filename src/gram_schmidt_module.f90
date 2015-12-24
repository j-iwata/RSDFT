MODULE gram_schmidt_module

  use gram_schmidt_m_module
  use gram_schmidt_t_module
 !use gram_schmidt_u_module
  use gram_schmidt_g_module
  use var_sys_parameter
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: gram_schmidt

  integer :: iswitch_algorithm = 0
  include 'mpif.h'

  logical :: flag_init_read = .true.

CONTAINS

  SUBROUTINE read_gram_schmidt
    implicit none
    call IOTools_readIntegerKeyword( "GS", iswitch_algorithm )
    flag_init_read = .false.
  END SUBROUTINE read_gram_schmidt

  SUBROUTINE gram_schmidt(n0,n1,k,s)

    implicit none
    integer,intent(IN) :: n0,n1,k,s

    call write_border( 1, " gram_schmidt(start)" )

    if ( flag_init_read ) call read_gram_schmidt

    if ( pp_kind == "USPP" ) then

       call gram_schmidt_g( n0,n1,k,s )

    else

       select case( iswitch_algorithm )
       case default
          call gram_schmidt_t(n0,n1,k,s)
       case( 1 )
          call gram_schmidt_m(n0,n1,k,s)
      !case( 2 )
      !   call gram_schmidt_u(n0,n1,k,s)
      end select

    end if

    call write_border( 1, " gram_schmidt(end)" )

  END SUBROUTINE gram_schmidt

END MODULE gram_schmidt_module
