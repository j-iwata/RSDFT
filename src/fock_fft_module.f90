MODULE fock_fft_module

  use parallel_module
  use rgrid_module, only: Ngrid, Igrid
  use ggrid_module, only: NGgrid, Ecut
  use bb_module, only: bb
  use xc_hybrid_module, only: q_fock, R_hf, omega &
                             ,iflag_lcwpbe, iflag_hse, iflag_hf, iflag_pbe0
  use watch_module
  use bz_module, only: kbb
  use fock_ffte_module, only: fock_ffte_double, fock_ffte, init_fock_ffte &
                             ,ct_fock_ffte,et_fock_ffte
  use fock_fftw_module, only: fock_fftw_double, fock_fftw, init_fock_fftw
  use fft_module

  implicit none

  PRIVATE
  PUBLIC :: ct_fock_fft, et_focK_fft
  PUBLIC :: Fock_FFT
  PUBLIC :: Fock_FFT_Double
  PUBLIC :: init_fock_fft

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN=MPI_REAL8
#else
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
#endif

  real(8) :: ct_fock_fft(10),et_fock_fft(10)

CONTAINS


  SUBROUTINE init_fock_fft
    implicit none
#ifdef _FFTE_
    call init_fock_ffte
#elif _FFTW_
    call init_fock_fftw
#endif
  END SUBROUTINE init_fock_fft


  SUBROUTINE Fock_fft( n1, n2, k, q, trho, tVh, t )

    implicit none

    integer,intent(IN) :: n1,n2,k,q,t
    integer :: i,i1,i2,i3,j1,j2,j3,ierr,irank,a1,a2,a3,b1,b2,b3
    integer :: ML1,ML2,ML3,ML
    real(8) :: pi,pi4,g2,const1,const2,k_fock(3)
    complex(8),parameter :: z0=(0.d0,0.d0)
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:)
#ifdef _DRSDFT_
    real(8),allocatable :: work(:)
    real(8),intent(IN)    :: trho(n1:n2)
    real(8),intent(INOUT) :: tVh(n1:n2)
#else
    complex(8),allocatable :: work(:)
    complex(8),intent(IN)    :: trho(n1:n2)
    complex(8),intent(INOUT) :: tVh(n1:n2)
#endif
    real(8) :: ctt(0:5),ett(0:5)

#ifdef _FFTE_

    ct_fock_ffte(:)=0.0d0
    et_fock_ffte(:)=0.0d0
    call Fock_FFTE( n1, n2, k, q, trho, tVh, t )
    ct_fock_fft(:)=ct_fock_fft(:)+ct_fock_ffte(:)
    et_fock_fft(:)=et_fock_fft(:)+et_fock_ffte(:)
    return

#elif _FFTW_

    call FocK_FFTW( n1, n2, k, q, trho, tVh, t )
    return

#endif

    pi  = acos(-1.0d0)
    pi4 = 4.0d0*pi

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    k_fock(1) = bb(1,1)*kbb(1,k) + bb(1,2)*kbb(2,k) + bb(1,3)*kbb(3,k)
    k_fock(2) = bb(2,1)*kbb(1,k) + bb(2,2)*kbb(2,k) + bb(2,3)*kbb(3,k)
    k_fock(3) = bb(3,1)*kbb(1,k) + bb(3,2)*kbb(2,k) + bb(3,3)*kbb(3,k)

    call watch(ctt(0),ett(0))

    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork0=z0

    allocate( work(ML) ) ; work=z0

    call mpi_allgatherv(trho(n1),ir_grid(myrank_g),TYPE_MAIN &
         ,work,ir_grid,id_grid,TYPE_MAIN,comm_grid,ierr)

    i=0
    irank=-1
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       a1=pinfo_grid(1,irank) ; b1=pinfo_grid(2,irank)+a1-1
       a2=pinfo_grid(3,irank) ; b2=pinfo_grid(4,irank)+a2-1
       a3=pinfo_grid(5,irank) ; b3=pinfo_grid(6,irank)+a3-1
       do j3=a3,b3
       do j2=a2,b2
       do j1=a1,b1
          i=i+1
          zwork0(j1,j2,j3)=work(i)
       end do
       end do
       end do
    end do
    end do
    end do

    deallocate( work )

!
! ---
!

    allocate( zwork1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork1=z0

    call init_fft

    call watch(ctt(1),ett(1))

    call forward_fft( zwork0, zwork1 )

    call watch(ctt(2),ett(2))

    zwork1(:,:,:)=z0

    if ( iflag_hf > 0 .or. iflag_pbe0 > 0 ) then
       do i3=-NGgrid(3),NGgrid(3)
       do i2=-NGgrid(2),NGgrid(2)
       do i1=-NGgrid(1),NGgrid(1)
          g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1,q,t))**2 &
            +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2,q,t))**2 &
            +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3,q,t))**2
          if ( g2 < 0.0d0 .or. g2 > Ecut ) cycle
          j1=mod(i1+ML1,ML1)
          j2=mod(i2+ML2,ML2)
          j3=mod(i3+ML3,ML3)
          if ( g2 <= 1.d-10 ) then
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*2.d0*Pi*R_hf**2
          else
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*pi4*(1.d0-cos(sqrt(g2)*R_hf))/g2
          end if
       end do
       end do
       end do
    end if

    if ( iflag_hse > 0 ) then
       const1 = 0.25d0/(omega*omega)
       const2 = pi/(omega*omega)
       do i3=-NGgrid(3),NGgrid(3)
       do i2=-NGgrid(2),NGgrid(2)
       do i1=-NGgrid(1),NGgrid(1)
          g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1,q,t))**2 &
            +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2,q,t))**2 &
            +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3,q,t))**2
          if ( g2 < 0.0d0 .or. g2 > Ecut ) cycle
          j1=mod(i1+ML1,ML1)
          j2=mod(i2+ML2,ML2)
          j3=mod(i3+ML3,ML3)
          if ( g2 <= 1.d-10 ) then
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*const2
          else
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*pi4*(1.0d0-exp(-g2*const1))/g2
          end if
       end do
       end do
       end do
    end if

    if ( iflag_lcwpbe > 0 ) then
       do i3=-NGgrid(3),NGgrid(3)
       do i2=-NGgrid(2),NGgrid(2)
       do i1=-NGgrid(1),NGgrid(1)
          g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1,q,t))**2 &
            +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2,q,t))**2 &
            +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3,q,t))**2
          if ( g2 < 0.0d0 .or. g2 > Ecut ) cycle
          j1=mod(i1+ML1,ML1)
          j2=mod(i2+ML2,ML2)
          j3=mod(i3+ML3,ML3)
          if ( g2 <= 1.d-10 ) then
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*(2.d0*Pi*R_hf**2.d0 &
                  -Pi/(omega**2.d0))
          else
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3) &
                  *pi4*(exp(-0.25d0*g2/(omega**2.d0))-cos(sqrt(g2)*R_hf))/g2
          end if
       end do
       end do
       end do
    end if

    call watch(ctt(3),ett(3))

    call backward_fft( zwork1, zwork0 )

    call watch(ctt(4),ett(4))

    call finalize_fft

    deallocate( zwork0 )

#ifdef _DRSDFT_        
    call z3_to_d1_fft( zwork1, tVh )
#else
    call z3_to_z1_fft( zwork1, tVh )
#endif

    deallocate( zwork1 )

    call watch(ctt(5),ett(5))

    ct_fock_fft(1) = ct_fock_fft(1) + ctt(1) - ctt(0)
    et_fock_fft(1) = et_fock_fft(1) + ett(1) - ett(0)
    ct_fock_fft(2) = ct_fock_fft(2) + ctt(2) - ctt(1)
    et_fock_fft(2) = et_fock_fft(2) + ett(2) - ett(1)
    ct_fock_fft(3) = ct_fock_fft(3) + ctt(3) - ctt(2)
    et_fock_fft(3) = et_fock_fft(3) + ett(3) - ett(2)
    ct_fock_fft(4) = ct_fock_fft(4) + ctt(4) - ctt(3)
    et_fock_fft(4) = et_fock_fft(4) + ett(4) - ett(3)
    ct_fock_fft(5) = ct_fock_fft(5) + ctt(5) - ctt(4)
    et_fock_fft(5) = et_fock_fft(5) + ett(5) - ett(4)

    return
  END SUBROUTINE Fock_FFT


  SUBROUTINE Fock_FFT_Double( n1, n2, trho, tVh )

    implicit none
    integer,intent(IN)     :: n1,n2
    complex(8),intent(IN)  :: trho(n1:n2)
    complex(8),intent(OUT) :: tVh(n1:n2)
    integer :: i,i1,i2,i3,j1,j2,j3,ierr,irank,a1,a2,a3,b1,b2,b3
    integer :: ML1,ML2,ML3,ML,k,q,t
    real(8) :: pi,pi4,g2,const1,const2,k_fock(3)
    complex(8),parameter :: z0=(0.d0,0.d0)
    complex(8),allocatable :: zwork0(:,:,:),zwork1(:,:,:)
    complex(8),allocatable :: work(:)
    real(8) :: ctt(0:5),ett(0:5)

#ifdef _FFTE_

    ct_fock_ffte(:)=0.0d0
    et_fock_ffte(:)=0.0d0
    call Fock_FFTE_Double( n1, n2, trho, tVh )
    ct_fock_fft(:)=ct_fock_fft(:)+ct_fock_ffte(:)
    et_fock_fft(:)=et_fock_fft(:)+et_fock_ffte(:)
    return

#elif _FFTW_

    call Fock_FFTW_Double( n1, n2, trho, tVh )
    return

#endif

    pi  = acos(-1.0d0)
    pi4 = 4.0d0*pi

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    k = 1
    q = 1
    t = 1

    k_fock(1) = bb(1,1)*kbb(1,k) + bb(1,2)*kbb(2,k) + bb(1,3)*kbb(3,k)
    k_fock(2) = bb(2,1)*kbb(1,k) + bb(2,2)*kbb(2,k) + bb(2,3)*kbb(3,k)
    k_fock(3) = bb(3,1)*kbb(1,k) + bb(3,2)*kbb(2,k) + bb(3,3)*kbb(3,k)

    call watch(ctt(0),ett(0))

    allocate( zwork0(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork0=z0

    allocate( work(ML) ) ; work=z0

    call mpi_allgatherv(trho(n1),ir_grid(myrank_g),MPI_COMPLEX16 &
         ,work,ir_grid,id_grid,MPI_COMPLEX16,comm_grid,ierr)

    i=0
    irank=-1
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       a1=pinfo_grid(1,irank) ; b1=pinfo_grid(2,irank)+a1-1
       a2=pinfo_grid(3,irank) ; b2=pinfo_grid(4,irank)+a2-1
       a3=pinfo_grid(5,irank) ; b3=pinfo_grid(6,irank)+a3-1
       do j3=a3,b3
       do j2=a2,b2
       do j1=a1,b1
          i=i+1
          zwork0(j1,j2,j3)=work(i)
       end do
       end do
       end do
    end do
    end do
    end do

    deallocate( work )

! ---

    allocate( zwork1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; zwork1=z0

    call init_fft

    call watch(ctt(1),ett(1))

    call forward_fft( zwork0, zwork1 )

    call watch(ctt(2),ett(2))

    zwork1(:,:,:)=z0

    if ( iflag_hf > 0 .or. iflag_pbe0 > 0 ) then
       do i3=-NGgrid(3),NGgrid(3)
       do i2=-NGgrid(2),NGgrid(2)
       do i1=-NGgrid(1),NGgrid(1)
          g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1,q,t))**2 &
            +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2,q,t))**2 &
            +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3,q,t))**2
          if ( g2 < 0.0d0 .or. g2 > Ecut ) cycle
          j1=mod(i1+ML1,ML1)
          j2=mod(i2+ML2,ML2)
          j3=mod(i3+ML3,ML3)
          if ( g2 <= 1.d-10 ) then
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*2.d0*Pi*R_hf**2
          else
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*pi4*(1.d0-cos(sqrt(g2)*R_hf))/g2
          end if
       end do
       end do
       end do
    end if

    if ( iflag_hse > 0 ) then
       const1 = 0.25d0/(omega*omega)
       const2 = pi/(omega*omega)
       do i3=-NGgrid(3),NGgrid(3)
       do i2=-NGgrid(2),NGgrid(2)
       do i1=-NGgrid(1),NGgrid(1)
          g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1,q,t))**2 &
            +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2,q,t))**2 &
            +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3,q,t))**2
          if ( g2 < 0.0d0 .or. g2 > Ecut ) cycle
          j1=mod(i1+ML1,ML1)
          j2=mod(i2+ML2,ML2)
          j3=mod(i3+ML3,ML3)
          if ( g2 <= 1.d-10 ) then
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*const2
          else
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*pi4*(1.0d0-exp(-g2*const1))/g2
          end if
       end do
       end do
       end do
    end if

    if ( iflag_lcwpbe > 0 ) then
       do i3=-NGgrid(3),NGgrid(3)
       do i2=-NGgrid(2),NGgrid(2)
       do i1=-NGgrid(1),NGgrid(1)
          g2=(bb(1,1)*i1+bb(1,2)*i2+bb(1,3)*i3+k_fock(1)-q_fock(1,q,t))**2 &
            +(bb(2,1)*i1+bb(2,2)*i2+bb(2,3)*i3+k_fock(2)-q_fock(2,q,t))**2 &
            +(bb(3,1)*i1+bb(3,2)*i2+bb(3,3)*i3+k_fock(3)-q_fock(3,q,t))**2
          if ( g2 < 0.0d0 .or. g2 > Ecut ) cycle
          j1=mod(i1+ML1,ML1)
          j2=mod(i2+ML2,ML2)
          j3=mod(i3+ML3,ML3)
          if ( g2 <= 1.d-10 ) then
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3)*(2.d0*Pi*R_hf**2.d0 &
                  -Pi/(omega**2.d0))
          else
             zwork1(j1,j2,j3)=zwork0(j1,j2,j3) &
                  *pi4*(exp(-0.25d0*g2/(omega**2.d0))-cos(sqrt(g2)*R_hf))/g2
          end if
       end do
       end do
       end do
    end if

    call watch(ctt(3),ett(3))

    call backward_fft( zwork1, zwork0 )

    call watch(ctt(4),ett(4))

    call finalize_fft

    deallocate( zwork0 )

    call z3_to_z1_fft( zwork1, tVh )

    deallocate( zwork1 )

    call watch(ctt(5),ett(5))

    ct_fock_fft(1) = ct_fock_fft(1) + ctt(1) - ctt(0)
    et_fock_fft(1) = et_fock_fft(1) + ett(1) - ett(0)
    ct_fock_fft(2) = ct_fock_fft(2) + ctt(2) - ctt(1)
    et_fock_fft(2) = et_fock_fft(2) + ett(2) - ett(1)
    ct_fock_fft(3) = ct_fock_fft(3) + ctt(3) - ctt(2)
    et_fock_fft(3) = et_fock_fft(3) + ett(3) - ett(2)
    ct_fock_fft(4) = ct_fock_fft(4) + ctt(4) - ctt(3)
    et_fock_fft(4) = et_fock_fft(4) + ett(4) - ett(3)
    ct_fock_fft(5) = ct_fock_fft(5) + ctt(5) - ctt(4)
    et_fock_fft(5) = et_fock_fft(5) + ett(5) - ett(4)

    return

  END SUBROUTINE Fock_FFT_Double


END MODULE fock_fft_module
