MODULE kinetic_mol_module

  use rgrid_mol_module, only: LL
  use bc_module, only: www, bcset, bcset_1, bcset_3
  use kinetic_variables, only: coef_lap0, coef_lap, coef_kin, Md

  implicit none

  PRIVATE
  PUBLIC :: op_kinetic_mol

CONTAINS

  SUBROUTINE op_kinetic_mol(n1,n2,ib1,ib2,tpsi,htpsi)
    implicit none
    integer,intent(IN) :: n1,n2,ib1,ib2
#ifdef _DRSDFT_
    real(8),intent(IN)    ::  tpsi(n1:n2,ib1:ib2)
    real(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    real(8) :: d
#else
    complex(8),intent(IN)    ::  tpsi(n1:n2,ib1:ib2)
    complex(8),intent(INOUT) :: htpsi(n1:n2,ib1:ib2)
    complex(8) :: d
#endif
    integer :: i,ib,i1,i2,i3,m,n
    real(8) :: c

    do ib=ib1,ib2
!$OMP do
       do i=n1,n2
          www( LL(1,i),LL(2,i),LL(3,i),ib-ib1+1 ) = tpsi(i,ib)
       end do
!$OMP end do
    end do

    call bcset_3(1,ib2-ib1+1,Md,0)
!$OMP barrier

    do ib=ib1,ib2
!$OMP do
       do i=n1,n2
          htpsi(i,ib) = htpsi(i,ib) + coef_lap0*tpsi(i,ib)
       end do
!$OMP end do
    end do

    do ib=ib1,ib2
       n=ib-ib1+1
!$OMP do
       do i=n1,n2
          i1 = LL(1,i)
          i2 = LL(2,i)
          i3 = LL(3,i)
          d  = htpsi(i,ib)
          do m=1,Md
             c = coef_kin(m)
             d = d + c*( www(i1-m,i2,i3,n)+www(i1+m,i2,i3,n) &
                       + www(i1,i2-m,i3,n)+www(i1,i2+m,i3,n) &
                       + www(i1,i2,i3-m,n)+www(i1,i2,i3+m,n) )
          end do
          htpsi(i,ib) = d
       end do
!$OMP end do
    end do ! ib

  END SUBROUTINE op_kinetic_mol

END MODULE kinetic_mol_module
