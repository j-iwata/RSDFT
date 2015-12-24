MODULE fock_fftw_module

  use parallel_module
  use rgrid_module, only: Ngrid, Igrid
  use ggrid_module, only: NGgrid, Ecut, LLG, construct_ggrid, destruct_ggrid
  use bb_module, only: bb
  use xc_hybrid_module, only: q_fock, R_hf, omega, gamma_hf &
                             ,iflag_lcwpbe, iflag_hse, iflag_hf, iflag_pbe0 &
                             ,n_kq_fock,i_kq_fock,kq_fock
  use bz_module, only: kbb
  use omp_variables
  use fft_module, only: init_fft, d1_to_z3_fft, z1_to_z3_fft
  use fftw_module, only: ML1_c, ML2_c, N_ML3_c, ML3_c0 &
       ,zwork3_ptr0, zwork3_ptr1, plan_forward, plan_backward &
       ,z3_to_d1_fftw, z3_to_z1_fftw
  use,intrinsic :: iso_c_binding

  implicit none

  PRIVATE
  PUBLIC :: fock_fftw
  PUBLIC :: fock_fftw_double
  PUBLIC :: init_fock_fftw

#ifdef _FFTW_
  include "fftw3-mpi.f03"
#endif

  logical :: first_time=.true.
  integer :: NGHT
  integer,allocatable :: LGHT(:,:)
  real(8),allocatable :: GGHT(:,:),GXYZ(:,:)

  complex(8),allocatable :: zwork3(:,:,:)

CONTAINS


  SUBROUTINE init_fock_fftw

    implicit none
    integer :: n,i,i1,i2,i3,j1,j2,j3,a2b,b2b,a3b,b3b
    real(8) :: g2,pi,c,d

    call write_border( 0, "init_fock_fftw(start)" )

    pi  = acos(-1.0d0)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)

    call construct_Ggrid(0)

    n=0
    do i=1,NGgrid(0)
       i1=mod( Ngrid(1)+LLG(1,i), Ngrid(1) )
       i2=mod( Ngrid(2)+LLG(2,i), Ngrid(2) )
       i3=mod( Ngrid(3)+LLG(3,i), Ngrid(3) )
       if ( ML3_c0 <= i3 .and. i3 <= ML3_c0+N_ML3_c-1 ) then
          n=n+1
       end if
    end do ! i

    if ( allocated(GGHT) ) deallocate( GGHT )
    if ( allocated(GXYZ) ) deallocate( GXYZ )
    if ( allocated(LGHT) ) deallocate( LGHT )
    allocate( LGHT(3,n) ) ; LGHT=0
    allocate( GXYZ(3,n) ) ; GXYZ=0.0d0
    allocate( GGHT(n,n_kq_fock) ) ; GGHT=0.0d0

    n=0
    do i=1,NGgrid(0)
       i1=LLG(1,i)
       i2=LLG(2,i)
       i3=LLG(3,i)
       j1=mod( Ngrid(1)+i1, Ngrid(1) )
       j2=mod( Ngrid(2)+i2, Ngrid(2) )
       j3=mod( Ngrid(3)+i3, Ngrid(3) )
       if ( ML3_c0 <= j3 .and. j3 <= ML3_c0+N_ML3_c-1 ) then
          n=n+1
          LGHT(1,n)=j1+1
          LGHT(2,n)=j2+1
          LGHT(3,n)=j3+1-ML3_c0
          GXYZ(1,n)=bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i)
          GXYZ(2,n)=bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i)
          GXYZ(3,n)=bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i)
       end if
    end do ! i
    NGHT=n

    d=1.0d0/Ngrid(0)

    do i=1,n_kq_fock

       do n=1,NGHT

          g2=( GXYZ(1,n) + kq_fock(1,i) )**2 &
            +( GXYZ(2,n) + kq_fock(2,i) )**2 &
            +( GXYZ(3,n) + kq_fock(3,i) )**2

          if ( g2 < 0.0d0 .or. g2 > Ecut ) cycle

          if ( iflag_hf > 0 .or. iflag_pbe0 > 0 ) then

             if ( g2 <= 1.d-10 ) then
                GGHT(n,i) = 2.0d0*pi*R_hf**2 * d
             else
                GGHT(n,i) = ( 1.0d0 - cos(sqrt(g2)*R_hf) )*4.0d0*pi/g2 * d
             end if

          else if ( iflag_hse > 0 ) then

             if ( g2 <= 1.d-10 ) then
                GGHT(n,i) = pi/(omega*omega) * d
             else
                c=0.25d0/(omega*omega)
                GGHT(n,i) = ( 1.0d0 - exp(-g2*c) )*4.0d0*pi/g2 * d
             end if

          end if

       end do ! n
    end do ! i

    call destruct_Ggrid

    call init_fft( .true. )

    first_time=.false.

    call write_border( 0, "init_fock_fftw(end)" )

  END SUBROUTINE init_fock_fftw


  SUBROUTINE fock_fftw( n1, n2, k, q, trho, tVh, t )

!$  use omp_lib
    implicit none
    integer,intent(IN) :: n1,n2,k,q,t
    integer :: i,i1,i2,i3,j,j1,j2,j3
#ifdef _DRSDFT_
    real(8),intent(IN)    :: trho(n1:n2)
    real(8),intent(INOUT) :: tVh(n1:n2)
#else
    complex(8),intent(IN)    :: trho(n1:n2)
    complex(8),intent(INOUT) :: tVh(n1:n2)
#endif

#ifdef _FFTW_

    if ( first_time ) call init_fock_fftw

#ifdef _DRSDFT_
    call d1_to_z3_fft( trho, zwork3 )
#else
    call z1_to_z3_fft( trho, zwork3 )
#endif

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

    zwork3_ptr0=(0.0d0,0.0d0)

    if ( gamma_hf == 1 ) then

!$OMP parallel do private( i1,i2,i3 )
       do i=1,NGHT
          i1=LGHT(1,i)
          i2=LGHT(2,i)
          i3=LGHT(3,i)
          zwork3_ptr0(i1,i2,i3)=zwork3_ptr1(i1,i2,i3)*GGHT(i,1)
       end do
!$OMP end parallel do

    else if ( gamma_hf == 0 ) then

       j=i_kq_fock(k,q,t)
!$OMP parallel do private( i1,i2,i3 )
       do i=1,NGHT
          i1=LGHT(1,i)
          i2=LGHT(2,i)
          i3=LGHT(3,i)
          zwork3_ptr0(i1,i2,i3)=zwork3_ptr1(i1,i2,i3)*GGHT(i,j)
       end do
!$OMP end parallel do

    end if ! [ gamma_hf ]

    call fftw_mpi_execute_dft( plan_backward, zwork3_ptr0, zwork3_ptr1 )

#ifdef _DRSDFT_
    call z3_to_d1_fftw( zwork3_ptr1, tVh )
#else
    call z3_to_z1_fftw( zwork3_ptr1, tVh )
#endif

#endif

  END SUBROUTINE fock_fftw


  SUBROUTINE fock_fftw_double( n1, n2, trho, tVh )

!$  use omp_lib
    implicit none
    integer,intent(IN)     :: n1,n2
    complex(8),intent(IN)  :: trho(n1:n2)
    complex(8),intent(OUT) :: tVh(n1:n2)
    integer :: i,i1,i2,i3,j1,j2,j3
    complex(8),parameter :: z0=(0.0d0,0.0d0)

#ifdef _FFTW_

    if ( first_time ) call init_fock_fftw

    call z1_to_z3_fft( trho, zwork3 )

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

    zwork3_ptr0=(0.0d0,0.0d0)

!$OMP parallel do private( i1,i2,i3 )
    do i=1,NGHT
       i1=LGHT(1,i)
       i2=LGHT(2,i)
       i3=LGHT(3,i)
       zwork3_ptr0(i1,i2,i3)=zwork3_ptr1(i1,i2,i3)*GGHT(i,1)
    end do
!$OMP end parallel do

    call fftw_mpi_execute_dft( plan_backward, zwork3_ptr0, zwork3_ptr1 )

    call z3_to_z1_fftw( zwork3_ptr1, tVh )

#else
    tVh=z0
#endif

  END SUBROUTINE fock_fftw_double


END MODULE fock_fftw_module
