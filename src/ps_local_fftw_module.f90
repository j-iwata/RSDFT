MODULE ps_local_fftw_module

  use rgrid_module, only: Ngrid
  use ggrid_module, only: NGgrid, MGL, MG_0,MG_1, LLG, allgatherv_ggrid &
                         ,construct_ggrid, destruct_ggrid
  use fftw_module, only: ML1_c, ML2_c, N_ML3_c, ML3_c0 &
       ,zwork3_ptr0, zwork3_ptr1, plan_backward, z3_to_d1_fftw
  use,intrinsic :: iso_c_binding

  implicit none

  PRIVATE
  PUBLIC :: construct_ps_local_fftw

CONTAINS


  SUBROUTINE construct_ps_local_fftw( vqlg, SGK, Vion )

    implicit none
    real(8),intent(IN) :: vqlg(:,:)
    complex(8),intent(IN) :: SGK(:,:)
    real(8),intent(OUT) :: Vion(:)
#ifdef _FFTW_
    integer :: i,i1,i2,i3,j1,j2,j3,ik,j,MG
    integer :: ML1,ML2,ML3,ML,Nelement
    complex(8),allocatable :: zwork3(:,:,:),vg(:)
    include "fftw3-mpi.f03"

    call write_border( 0, " construct_ps_local_fftw(start)" )

    MG  = NGgrid(0)
    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    Nelement = size( SGK, 2 )

    allocate( vg(MG) ) ; vg=(0.0d0,0.0d0)

    do ik=1,Nelement
       do i=MG_0,MG_1
          j=MGL(i)
          vg(i)=vg(i)+vqlg(j,ik)*SGK(i,ik)
       end do
    end do
    call allgatherv_Ggrid(vg)

    call construct_Ggrid(2)

    allocate( zwork3(0:ML1-1,0:ML2-1,0:ML3-1) )
    zwork3(:,:,:)=(0.0d0,0.0d0)

    do i=1,NGgrid(0)
       zwork3(LLG(1,i),LLG(2,i),LLG(3,i))=vg(i)
    end do

    call destruct_Ggrid

    deallocate( vg )

    do i3=1,N_ML3_c
       do i2=1,ML2_c
       do i1=1,ML1_c
          j1=i1-1
          j2=i2-1
          j3=i3-1+ML3_c0
          zwork3_ptr0(i1,i2,i3) = zwork3(j1,j2,j3)
       end do
       end do
    end do

    call fftw_mpi_execute_dft( plan_backward, zwork3_ptr0, zwork3_ptr1 )

    call z3_to_d1_fftw( zwork3_ptr1, Vion )

    deallocate( zwork3 )

    call write_border( 0, " construct_ps_local_fftw(end)" )
#else
    Vion=0.0d0
#endif
  END SUBROUTINE construct_ps_local_fftw


END MODULE ps_local_fftw_module
