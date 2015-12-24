MODULE hartree_sol_fftw_module

  use hartree_variables, only: E_hartree, Vh
  use bb_module, only: bb
  use rgrid_module, only: Ngrid, Igrid, dV
  use ggrid_module, only: NGgrid, LLG, construct_ggrid, destruct_ggrid
  use grid_module, only: inner_product_grid
  use fftw_module, only: ML1_c, ML2_c, ML3_c0, N_ML3_c &
       , zwork3_ptr0, zwork3_ptr1, z3_to_d1_fftw &
       , plan_forward, plan_backward
  use fft_module, only: d1_to_z3_fft, init_fft
  use,intrinsic :: iso_c_binding

  implicit none

  PRIVATE
  PUBLIC :: calc_hartree_sol_fftw

  logical :: first_time=.true.
  integer :: NGHT
  integer,allocatable :: LGHT(:,:)
  integer,allocatable :: IGHT(:,:)
  real(8),allocatable :: GGHT(:)

CONTAINS


  SUBROUTINE calc_hartree_sol_fftw(n1,n2,n3,rho)
    implicit none
    integer,intent(IN) :: n1,n2,n3
    real(8),intent(IN) :: rho(n1:n2,n3)
#ifdef _FFTW_
    integer :: i,i1,i2,i3,j1,j2,j3,n
    real(8) :: Eh0,g2,const
    complex(8),parameter :: z0=(0.0d0,0.0d0)
    real(8),allocatable :: work(:)
    complex(8),allocatable :: zwork3(:,:,:)
    include 'fftw3-mpi.f03'

    call write_border( 1, " calc_hartree_sol_fftw(start)" )

    if ( first_time ) then
       call construct_Ggrid(0)
       n=0
       do i=1,NGgrid(0)
          i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
          i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
          i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
          if ( all(LLG(1:3,i)==0) ) cycle
          if ( ML3_c0 <= i3 .and. i3 <= ML3_c0+N_ML3_c-1 ) then
             n=n+1
          end if
       end do
       allocate( LGHT(3,n) ) ; LGHT=0
       allocate( IGHT(3,n) ) ; IGHT=0
       allocate( GGHT(n)   ) ; GGHT=0.0d0
       const = 4.0d0*acos(-1.d0)/Ngrid(0)
       n=0
       do i=1,NGgrid(0)
          i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
          i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
          i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
          if ( all(LLG(1:3,i)==0) ) cycle
          if ( ML3_c0 <= i3 .and. i3 <= ML3_c0+N_ML3_c-1 ) then
             n=n+1
             LGHT(1,n)=i1+1
             LGHT(2,n)=i2+1
             LGHT(3,n)=i3+1-ML3_c0
             g2=( bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i) )**2 &
               +( bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i) )**2 &
               +( bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i) )**2
             GGHT(n)=const/g2
          end if
       end do
       NGHT=n
       call destruct_Ggrid
       call init_fft( .true. )
       first_time=.false.
    end if

    allocate( work(size(rho,1)) ) ; work(:)=rho(:,1)
    do n=2,n3
       work(:)=work(:)+rho(:,n)
    end do
    call d1_to_z3_fft( work, zwork3 )

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

    call fftw_mpi_execute_dft( plan_forward, zwork3_ptr0, zwork3_ptr1 )

    zwork3_ptr0(:,:,:)=z0

    do i=1,NGHT
       i1=LGHT(1,i)
       i2=LGHT(2,i)
       i3=LGHT(3,i)
       zwork3_ptr0(i1,i2,i3)=zwork3_ptr1(i1,i2,i3)*GGHT(i)
    end do

    call fftw_mpi_execute_dft( plan_backward, zwork3_ptr0, zwork3_ptr1 )

    call z3_to_d1_fftw( zwork3_ptr1, Vh )

    call inner_product_grid( work, Vh, 0.5d0*dV, E_hartree ) 

    deallocate( zwork3 )
    deallocate( work )

    call write_border( 1, " calc_hartree_sol_fftw(end)" )
#endif
  END SUBROUTINE calc_hartree_sol_fftw


END MODULE hartree_sol_fftw_module
