MODULE strfac_module

  use ggrid_module
  use atom_module

  implicit none

  PRIVATE
  PUBLIC :: SGK, construct_strfac, destruct_strfac &
           ,construct_strfac_2

  complex(8),allocatable :: SGK(:,:)

CONTAINS

  SUBROUTINE construct_strfac

    integer :: a,i,i1,i2,i3,ik,j,ierr,MG
    real(8) :: Gr,pi2,a1,a2,a3

    call write_border( 80, " construct_strfac(start)" )

    pi2 = 2.d0*acos(-1.d0)
    MG  = NGgrid(0)

    allocate( SGK(MG,Nelement) ) ; SGK=(0.d0,0.d0)

    call construct_Ggrid(1)

    do a=1,Natom
       ik=ki_atom(a)
       a1=pi2*aa_atom(1,a)
       a2=pi2*aa_atom(2,a)
       a3=pi2*aa_atom(3,a)
       do i=MG_0,MG_1
          Gr=LLG(1,i)*a1+LLG(2,i)*a2+LLG(3,i)*a3
          SGK(i,ik)=SGK(i,ik)+dcmplx(cos(Gr),-sin(Gr))
       end do
    end do

    call destruct_Ggrid

    call write_border( 80, " construct_strfac(end)" )

  END SUBROUTINE construct_strfac

  SUBROUTINE destruct_strfac
    call write_border( 80, " destruct_strfac(start)" )
    deallocate( SGK )
    call write_border( 80, " destruct_strfac(end)" )
  END SUBROUTINE destruct_strfac


  SUBROUTINE construct_strfac_2(MG2_0,MG2_1)
    implicit none
    integer,intent(IN) :: MG2_0,MG2_1
    integer :: a,ik,i
    real(8) :: a1,a2,a3,Gr,pi2
    pi2=2.0d0*acos(-1.0d0)
    allocate( SGK(MG2_0:MG2_1,Nelement) )
    SGK=(0.0d0,0.0d0)
    do a=1,Natom
       ik=ki_atom(a)
       a1=pi2*aa_atom(1,a)
       a2=pi2*aa_atom(2,a)
       a3=pi2*aa_atom(3,a)
       do i=MG2_0,MG2_1
          Gr=LLG(1,i)*a1+LLG(2,i)*a2+LLG(3,i)*a3
          SGK(i,ik)=SGK(i,ik)+dcmplx(cos(Gr),-sin(Gr))
       end do
    end do
  END SUBROUTINE construct_strfac_2


END MODULE strfac_module
