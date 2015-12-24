MODULE fock_ffte_module

  use parallel_module
  use rgrid_module, only: Ngrid, Igrid
  use ggrid_module, only: NGgrid, Ecut, LLG, construct_ggrid, destruct_ggrid
  use bb_module, only: bb
  use xc_hybrid_module, only: q_fock, R_hf, omega, gamma_hf &
                             ,iflag_lcwpbe, iflag_hse, iflag_hf, iflag_pbe0 &
                             ,n_kq_fock,i_kq_fock,kq_fock
  use watch_module
  use bz_module, only: kbb
  use ffte_sub_module, only: zwork1_ffte, zwork2_ffte, npuz, npuy, npux &
                            ,comm_fftx, comm_ffty, comm_fftz
  use omp_variables

  implicit none

  PRIVATE
  PUBLIC :: ct_fock_ffte,et_fock_ffte, fock_ffte, init_fock_ffte
  PUBLIC :: fock_ffte_double

  real(8) :: ct_fock_ffte(10),et_fock_ffte(10)

  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
  logical :: first_time=.true.
  integer :: NGHT
  integer,allocatable :: LGHT(:,:)
  real(8),allocatable :: GGHT(:,:),GXYZ(:,:)

CONTAINS


  SUBROUTINE init_fock_ffte

    implicit none
    integer :: n,i,i1,i2,i3,j1,j2,j3,a2b,b2b,a3b,b3b
    real(8) :: g2,pi,c

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
       if ( a2b <= i2 .and. i2 <= b2b .and. a3b <= i3 .and. i3 <= b3b ) then
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
       if ( a2b <= j2 .and. j2 <= b2b .and. a3b <= j3 .and. j3 <= b3b ) then
          n=n+1
          LGHT(1,n)=j1
          LGHT(2,n)=j2
          LGHT(3,n)=j3
          GXYZ(1,n)=bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i)
          GXYZ(2,n)=bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i)
          GXYZ(3,n)=bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i)
       end if
    end do ! i
    NGHT=n

    do i=1,n_kq_fock

       do n=1,NGHT

          g2=( GXYZ(1,n) + kq_fock(1,i) )**2 &
            +( GXYZ(2,n) + kq_fock(2,i) )**2 &
            +( GXYZ(3,n) + kq_fock(3,i) )**2

          if ( g2 < 0.0d0 .or. g2 > Ecut ) cycle

          if ( iflag_hf > 0 .or. iflag_pbe0 > 0 ) then

             if ( g2 <= 1.d-10 ) then
                GGHT(n,i) = 2.0d0*pi*R_hf**2
             else
                GGHT(n,i) = ( 1.0d0 - cos(sqrt(g2)*R_hf) )*4.0d0*pi/g2
             end if

          else if ( iflag_hse > 0 ) then

             if ( g2 <= 1.d-10 ) then
                GGHT(n,i) = pi/(omega*omega)
             else
                c=0.25d0/(omega*omega)
                GGHT(n,i) = ( 1.0d0 - exp(-g2*c) )*4.0d0*pi/g2
             end if

          end if

       end do ! n
    end do ! i

    call destruct_Ggrid

    first_time=.false.

  END SUBROUTINE init_fock_ffte


  SUBROUTINE fock_ffte( n1, n2, k, q, trho, tVh, t )

!$  use omp_lib
    implicit none

    integer,intent(IN) :: n1,n2,k,q,t
    integer :: i,i1,i2,i3,j1,j2,j3,ierr,irank,a1,a2,a3,b1,b2,b3
    integer :: ML1,ML2,ML3,ML,j
    real(8) :: pi,pi4,g2,const1,const2,k_fock(3)
    complex(8),parameter :: z0=(0.d0,0.d0)
#ifdef _DRSDFT_
    real(8),allocatable :: work(:)
    real(8),intent(IN)    :: trho(n1:n2)
    real(8),intent(INOUT) :: tVh(n1:n2)
#else
    complex(8),allocatable :: work(:)
    complex(8),intent(IN)    :: trho(n1:n2)
    complex(8),intent(INOUT) :: tVh(n1:n2)
#endif
    integer :: ispin,n,mm
    real(8) :: ctt(0:5),ett(0:5)
    integer :: MG,ML_0,ML_1,a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12
    integer :: MG1,MG2,MG3,NG1,NG2,NG3,np1,np2,np3

    pi  = acos(-1.0d0)
    pi4 = 4.0d0*pi

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    k_fock(1) = bb(1,1)*kbb(1,k) + bb(1,2)*kbb(2,k) + bb(1,3)*kbb(3,k)
    k_fock(2) = bb(2,1)*kbb(1,k) + bb(2,2)*kbb(2,k) + bb(2,3)*kbb(3,k)
    k_fock(3) = bb(3,1)*kbb(1,k) + bb(3,2)*kbb(2,k) + bb(3,3)*kbb(3,k)

    MG  = NGgrid(0)
    ML_0= Igrid(1,0)
    ML_1= Igrid(2,0)
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = (b1b-a1b+1)
    ab12= (b1b-a1b+1)*(b2b-a2b+1)

    mm = 0

    if ( first_time ) call init_fock_ffte
    
    ctt(:)=0.0d0
    ett(:)=0.0d0

    call watch(ctt(0),ett(0))

!$OMP parallel private( mm,i1,i2,i3,i )
!$  mm=omp_get_thread_num()
!$OMP workshare
    zwork1_ffte(:,:,:)=z0
!$OMP end workshare
    do i3=Igrid_omp(1,3,mm),Igrid_omp(2,3,mm)
    do i2=Igrid_omp(1,2,mm),Igrid_omp(2,2,mm)
    do i1=Igrid_omp(1,1,mm),Igrid_omp(2,1,mm)
       i=ML_0+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
       zwork1_ffte(i1,i2,i3) = trho(i)
    end do
    end do
    end do
!$OMP end parallel

    call mpi_allreduce(zwork1_ffte,zwork2_ffte,ML1*(b2b-a2b+1)*(b3b-a3b+1) &
         ,mpi_complex16,mpi_sum,comm_fftx,ierr)

    call watch(ctt(1),ett(1))

    call pzfft3dv(zwork2_ffte,zwork1_ffte,ML1,ML2,ML3 &
         ,comm_ffty,comm_fftz,npuy,npuz,-1)

    call watch(ctt(2),ett(2))

    zwork2_ffte(:,:,:)=z0

    if ( gamma_hf == 1 ) then

!$OMP parallel do private( i1,i2,i3 )
       do i=1,NGHT
          i1=LGHT(1,i)
          i2=LGHT(2,i)
          i3=LGHT(3,i)
          zwork2_ffte(i1,i2,i3)=zwork1_ffte(i1,i2,i3)*GGHT(i,1)
       end do
!$OMP end parallel do

    else if ( gamma_hf == 0 ) then

       j=i_kq_fock(k,q,t)
!$OMP parallel do private( i1,i2,i3 )
       do i=1,NGHT
          i1=LGHT(1,i)
          i2=LGHT(2,i)
          i3=LGHT(3,i)
          zwork2_ffte(i1,i2,i3)=zwork1_ffte(i1,i2,i3)*GGHT(i,j)
       end do
!$OMP end parallel do

    end if ! [ gamma_hf ]

    call watch(ctt(3),ett(3))

    call pzfft3dv(zwork2_ffte,zwork1_ffte,ML1,ML2,ML3 &
         ,comm_ffty,comm_fftz,npuy,npuz,1)

    call watch(ctt(4),ett(4))

!$OMP parallel private( mm,i,i1,i2,i3 )
!$  mm=omp_get_thread_num()
    do i3=Igrid_omp(1,3,mm),Igrid_omp(2,3,mm)
    do i2=Igrid_omp(1,2,mm),Igrid_omp(2,2,mm)
    do i1=Igrid_omp(1,1,mm),Igrid_omp(2,1,mm)
       i=ML_0+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
#ifdef _DRSDFT_
       tVh(i) = real( zwork1_ffte(i1,i2,i3) )
#else
       tVh(i) = zwork1_ffte(i1,i2,i3)
#endif
    end do
    end do
    end do
!$OMP end parallel

    call watch(ctt(5),ett(5))

    ct_fock_ffte(1) = ctt(1) - ctt(0)
    et_fock_ffte(1) = ett(1) - ett(0)
    ct_fock_ffte(2) = ctt(2) - ctt(1)
    et_fock_ffte(2) = ett(2) - ett(1)
    ct_fock_ffte(3) = ctt(3) - ctt(2)
    et_fock_ffte(3) = ett(3) - ett(2)
    ct_fock_ffte(4) = ctt(4) - ctt(3)
    et_fock_ffte(4) = ett(4) - ett(3)
    ct_fock_ffte(5) = ctt(5) - ctt(4)
    et_fock_ffte(5) = ett(5) - ett(4)

  END SUBROUTINE fock_ffte


  SUBROUTINE fock_ffte_double( n1, n2, trho, tVh )
!$  use omp_lib
    implicit none
    integer,intent(IN)     :: n1,n2
    complex(8),intent(IN)  :: trho(n1:n2)
    complex(8),intent(OUT) :: tVh(n1:n2)
    complex(8),parameter :: z0=(0.d0,0.d0)
    real(8) :: ctt(0:5),ett(0:5)
    integer :: ML1,ML2,ML3,mm,i,i1,i2,i3,ierr
    integer :: ML_0,ML_1,a1b,b1b,a2b,b2b,a3b,b3b,ab1,ab12

    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

    ML_0= Igrid(1,0)
    ML_1= Igrid(2,0)
    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)
    ab1 = (b1b-a1b+1)
    ab12= (b1b-a1b+1)*(b2b-a2b+1)

    mm = 0

    if ( first_time ) call init_fock_ffte
    
    ctt(:)=0.0d0
    ett(:)=0.0d0

    call watch(ctt(0),ett(0))

!$OMP parallel private( mm,i1,i2,i3,i )
!$  mm=omp_get_thread_num()
!$OMP workshare
    zwork1_ffte(:,:,:)=z0
!$OMP end workshare
    do i3=Igrid_omp(1,3,mm),Igrid_omp(2,3,mm)
    do i2=Igrid_omp(1,2,mm),Igrid_omp(2,2,mm)
    do i1=Igrid_omp(1,1,mm),Igrid_omp(2,1,mm)
       i=ML_0+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
       zwork1_ffte(i1,i2,i3) = trho(i)
    end do
    end do
    end do
!$OMP end parallel

    call mpi_allreduce(zwork1_ffte,zwork2_ffte,ML1*(b2b-a2b+1)*(b3b-a3b+1) &
         ,mpi_complex16,mpi_sum,comm_fftx,ierr)

    call watch(ctt(1),ett(1))

    call pzfft3dv(zwork2_ffte,zwork1_ffte,ML1,ML2,ML3 &
         ,comm_ffty,comm_fftz,npuy,npuz,-1)

    call watch(ctt(2),ett(2))

    zwork2_ffte(:,:,:)=z0

!$OMP parallel do private( i1,i2,i3 )
    do i=1,NGHT
       i1=LGHT(1,i)
       i2=LGHT(2,i)
       i3=LGHT(3,i)
       zwork2_ffte(i1,i2,i3)=zwork1_ffte(i1,i2,i3)*GGHT(i,1)
    end do
!$OMP end parallel do

    call watch(ctt(3),ett(3))

    call pzfft3dv(zwork2_ffte,zwork1_ffte,ML1,ML2,ML3 &
         ,comm_ffty,comm_fftz,npuy,npuz,1)

    call watch(ctt(4),ett(4))

!$OMP parallel private( mm,i,i1,i2,i3 )
!$  mm=omp_get_thread_num()
    do i3=Igrid_omp(1,3,mm),Igrid_omp(2,3,mm)
    do i2=Igrid_omp(1,2,mm),Igrid_omp(2,2,mm)
    do i1=Igrid_omp(1,1,mm),Igrid_omp(2,1,mm)
       i=ML_0+(i1-a1b)+(i2-a2b)*ab1+(i3-a3b)*ab12
       tVh(i) = zwork1_ffte(i1,i2,i3)
    end do
    end do
    end do
!$OMP end parallel

    call watch(ctt(5),ett(5))

    ct_fock_ffte(1) = ctt(1) - ctt(0)
    et_fock_ffte(1) = ett(1) - ett(0)
    ct_fock_ffte(2) = ctt(2) - ctt(1)
    et_fock_ffte(2) = ett(2) - ett(1)
    ct_fock_ffte(3) = ctt(3) - ctt(2)
    et_fock_ffte(3) = ett(3) - ett(2)
    ct_fock_ffte(4) = ctt(4) - ctt(3)
    et_fock_ffte(4) = ett(4) - ett(3)
    ct_fock_ffte(5) = ctt(5) - ctt(4)
    et_fock_ffte(5) = ett(5) - ett(4)

  END SUBROUTINE fock_ffte_double


END MODULE fock_ffte_module
