MODULE fftw_module

  use,intrinsic :: iso_c_binding
  use grid_module, only: grid, get_range_rgrid

  implicit none

  PRIVATE
  PUBLIC :: init_fftw
  PUBLIC :: finalize_fftw
  PUBLIC :: plan_forward, plan_backward
  PUBLIC :: ML1_c, ML2_c, N_ML3_c, ML3_c0
  PUBLIC :: zwork3_ptr0, zwork3_ptr1
  PUBLIC :: z3_to_d1_fftw, z3_to_z1_fftw

  integer :: comm_fftw
  type(c_ptr) :: plan_forward,plan_backward
  type(c_ptr) :: zwork3_cptr0,zwork3_cptr1
  integer(c_intptr_t) :: ML1_c,ML2_c,ML3_c,ML3_c0,N_ML3_c
  complex(c_double_complex), pointer :: zwork3_ptr0(:,:,:)=>null()
  complex(c_double_complex), pointer :: zwork3_ptr1(:,:,:)=>null()

CONTAINS


  SUBROUTINE init_fftw( Ngrid, Np, comm_grid, myrank_g )

    implicit none
    integer,intent(IN) :: Ngrid(3),Np(3),comm_grid, myrank_g
#ifdef _FFTW_
    integer :: ML1,Ml2,ML3,np1,np2,np3
    integer :: i1,i2,i3,irank,icolor,ierr
    integer(c_intptr_t) :: alloc_ML3_c
    include "fftw3-mpi.f03"

    call write_border( 0, " init_fftw(start)" )

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    np1 = Np(1)
    np2 = Np(2)
    np3 = Np(3)

    if ( mod(ML3,np3) /= 0 ) then
       call stop_program( "fftw_init: mod(ML3,np3)/=0 is not supported." )
    end if

    irank=-1
    do i3=1,np3
    do i2=1,np2
    do i1=1,np1
       irank=irank+1
       if ( irank == myrank_g ) icolor=i1+(i2-1)*np1
    end do
    end do
    end do

    call mpi_comm_split(comm_grid,icolor,0,comm_fftw,ierr)

    call fftw_mpi_init()

    ML1_c = ML1
    ML2_c = ML2
    ML3_c = ML3

    alloc_ML3_c = &
         fftw_mpi_local_size_3d(ML3_c,ML2_c,ML1_c,comm_fftw,N_ML3_c,ML3_c0)
    zwork3_cptr0=fftw_alloc_complex(alloc_ML3_c)
    zwork3_cptr1=fftw_alloc_complex(alloc_ML3_c)
    call c_f_pointer(zwork3_cptr0, zwork3_ptr0, [ML1_c,ML2_c,N_ML3_c])
    call c_f_pointer(zwork3_cptr1, zwork3_ptr1, [ML1_c,ML2_c,N_ML3_c])
    zwork3_ptr0=(0.0d0,0.0d0)
    zwork3_ptr1=(0.0d0,0.0d0)

    plan_forward  = fftw_mpi_plan_dft_3d( ML3_c,ML2_c,ML1_c &
         ,zwork3_ptr0,zwork3_ptr1,comm_fftw,fftw_forward,fftw_measure )

    plan_backward = fftw_mpi_plan_dft_3d( ML3_c,ML2_c,ML1_c &
         ,zwork3_ptr0,zwork3_ptr1,comm_fftw,fftw_backward,fftw_measure )

    call write_border( 0, " init_fftw(end)" )
#endif
  END SUBROUTINE init_fftw


  SUBROUTINE finalize_fftw
#ifdef _FFTW_
    implicit none
    integer :: ierr
    include "fftw3-mpi.f03"
    call mpi_comm_free(comm_fftw,ierr)
    call fftw_destroy_plan(plan_forward)
    call fftw_destroy_plan(plan_backward)
    call fftw_free(zwork3_cptr0)
    call fftw_free(zwork3_cptr1)
    call fftw_mpi_cleanup()
#endif
  END SUBROUTINE finalize_fftw


  SUBROUTINE z3_to_d1_fftw( z3, d1 )
    implicit none
    complex(8),intent(IN) :: z3(:,:,:)
    real(8),intent(OUT) :: d1(:)
    integer :: i1,i2,i3,i
    type(grid) :: rgrid
    call get_range_rgrid( rgrid )
    i=0
    do i3=1,N_ML3_c
    do i2=rgrid%g3%y%head+1,rgrid%g3%y%tail+1
    do i1=rgrid%g3%x%head+1,rgrid%g3%x%tail+1
       i=i+1
       d1(i) = dble( z3(i1,i2,i3) )
    end do
    end do
    end do
  END SUBROUTINE z3_to_d1_fftw

  SUBROUTINE z3_to_z1_fftw( z3, z1 )
    implicit none
    complex(8),intent(IN) :: z3(:,:,:)
    complex(8),intent(OUT) :: z1(:)
    integer :: i1,i2,i3,i
    type(grid) :: rgrid
    call get_range_rgrid( rgrid )
    i=0
    do i3=1,N_ML3_c
    do i2=rgrid%g3%y%head+1,rgrid%g3%y%tail+1
    do i1=rgrid%g3%x%head+1,rgrid%g3%x%tail+1
       i=i+1
       z1(i) = z3(i1,i2,i3)
    end do
    end do
    end do
  END SUBROUTINE z3_to_z1_fftw


END MODULE fftw_module
