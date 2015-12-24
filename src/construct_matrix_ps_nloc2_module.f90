MODULE construct_matrix_ps_nloc2_module

  use ps_nloc2_variables
  use rgrid_module

  implicit none

  PRIVATE
  PUBLIC :: construct_matrix_ps_nloc2

CONTAINS

  SUBROUTINE construct_matrix_ps_nloc2( k, ML, Hmat )
    implicit none
    integer,intent(IN) :: k, ML
#ifdef _DRSDFT_
    real(8),intent(INOUT) :: Hmat(ML,ML)
#else
    complex(8),intent(INOUT) :: Hmat(ML,ML)
#endif
    integer :: i1,i2,j1,j2,lma
    real(8) :: c

    do lma=1,nzlma

       c = iuV(lma)*dV

       do j2=1,MJJ(lma)

          i2 = JJP(j2,lma)

          do j1=1,MJJ(lma)

             i1 = JJP(j1,lma)

#ifdef _DRSDFT_
             Hmat(i1,i2) = Hmat(i1,i2) + c*uVk(j1,lma,k)*uVk(j2,lma,k)
#else
             Hmat(i1,i2) = Hmat(i1,i2) + c*uVk(j1,lma,k)*conjg(uVk(j2,lma,k))
#endif

          end do ! j2

       end do ! j1

    end do ! lma

  END SUBROUTINE construct_matrix_ps_nloc2


END MODULE construct_matrix_ps_nloc2_module
