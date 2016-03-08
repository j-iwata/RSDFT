PROGRAM dos_simple

  implicit none

  integer,parameter :: u5=5, u99=99, u10=10
  integer,parameter :: max_loop=10000000
  integer :: mbv,mbc,ne,nbk,natom,ie,mb,msp,s
  integer :: func_type, n_MP
  real(8),parameter :: HT=27.2116d0
  real(8),allocatable :: eval(:,:,:),occp(:,:,:)
  real(8),allocatable :: factor(:)
  real(8) :: gamma,evb,ecb,e1,e2,e,de,emin,emax,f(2),c1,c2

  natom = 0
  mbv   = 0
  mb    = 0
  nbk   = 0
  gamma = 0.0d0

!------------------------------

  write(*,*) "# of atoms: natom="
  read(u5,*) natom
  if ( natom <= 0 ) natom=1
  write(*,*) natom

  write(*,*) "band index of energy origin: mbv="
  read(u5,*) mbv
  write(*,*) mbv

  write(*,*) "width of smearing (eV): gamma="
  read(u5,*) gamma
  if ( gamma < 1.d-5 ) then
     gamma = 1.d-4
     write(*,*) "gamma(eV) is replaced to ",gamma
  end if
  write(*,*) gamma

!------------------------------

  call read_from_fort99(u99)
!  call read_from_band_eigv("band_eigv")

!------------------------------

  emin = minval( eval )
  emax = maxval( eval )

  write(*,*) "emin     (HT,eV)=",emin,emin*HT
  write(*,*) "emax     (HT,eV)=",emax,emax*HT
  write(*,*) "emax-emin(HT,eV)=",emax-emin,(emax-emin)*HT

!------------------------------

  if ( 1 <= mbv .and. mbv <= mb ) then

     mbc = mbv + 1
     evb = maxval( eval(1:mbv,1:nbk,1:msp) )
     ecb = minval( eval(mbc:mb,1:nbk,1:msp) )
     eval=eval-evb

     emin = minval( eval )
     emax = maxval( eval )

     write(*,*) "Energy is shifted by ", evb
     write(*,*) "emin   (HT,eV)=",emin,emin*HT
     write(*,*) "emax   (HT,eV)=",emax,emax*HT
     write(*,*) "ecb-evb(HT,eV)=",ecb-evb,(ecb-evb)*HT

  end if

!------------------------------

  write(*,*) "# of energy grid: ne= ( 0: default value is set )"
  read(u5,*) ne
  write(*,*) ne
  write(*,*) "energy range: e1,e2= ( 0,0: default values are set )"
  read(u5,*) e1,e2
  write(*,*) e1,e2
  write(*,*) "func type [0:gaussian, 1:Lorentzian, 2:Methfessel-Paxton]"
  read(u5,*) func_type
  if ( func_type < 0 .or. 2 < func_type ) then
     func_type = 0
     write(*,*) "func_type is replaced to ",func_type
  end if
  write(*,*) func_type
  if ( func_type == 2 ) then
     write(*,*) "Order of MP =?"
     read(u5,*) n_MP
     write(*,*) n_MP
     call init_gg0( n_MP )
  end if

  if ( ne <= 0 .or. e1 >= e2 ) then
     ne = 2000
     e1 = emin
     e2 = emax
  end if

  de = (e2-e1)/ne

  gamma = gamma/HT

  c1 = HT
  c2 = 1.0d0/(HT*natom)

  rewind u10

  do ie=1,ne+1

     e=e1+(ie-1)*de

     select case( func_type )
     case( 0 )
        do s=1,msp
           call dos0(s,e,f(s))
        end do
     case( 1 )
        do s=1,msp
           call dos1(s,e,f(s))
        end do
     case( 2 )
        do s=1,msp
           call dos2(s,e,f(s),n_MP)
        end do
     end select

     write(u10,'(1x,4f20.15)') e*c1, f(1:msp)*c2, sum(f(1:msp))*c2

  end do ! ie

CONTAINS


  SUBROUTINE dos0(s,e,f)
    implicit none
    integer,intent(IN)  :: s
    real(8),intent(IN)  :: e
    real(8),intent(OUT) :: f
    real(8) :: a,b,c,x
    integer :: i,k

    a = ( 1.d0/gamma )**2
    b = sqrt( a/acos(-1.0d0) )

    f=0.0d0
    do k=1,nbk
       do i=1,mb

          x=e-eval(i,k,s)
          c=b*occp(i,k,s)
          f=f+c*exp(-a*x*x)

       end do
    end do

  END SUBROUTINE dos0


  SUBROUTINE dos1(s,e,f)
    implicit none
    integer,intent(IN)  :: s
    real(8),intent(IN)  :: e
    real(8),intent(OUT) :: f
    real(8) :: b,c,x
    integer :: i,k

    b = 1.0d0/acos(-1.0d0)

    f=0.0d0
    do k=1,nbk
       do i=1,mb

          x=e-eval(i,k,s)
          c=b*occp(i,k,s)
          f=f+c*gamma/(x*x+gamma*gamma)

       end do
    end do

  END SUBROUTINE dos1


  SUBROUTINE read_from_fort99(u)
    implicit none
    integer,intent(IN) :: u
    integer :: n,k,i,j,s
    real(8) :: dummy(3)
    character(11) :: cbuf

    rewind u
    do i=1,max_loop
       read(u,*,end=999) cbuf
       if ( cbuf == "Eigenvalues" ) exit
    end do
    read(u,*)
    do i=1,max_loop
       read(u,*,end=99) k,n,dummy(1:3)
    end do
    if ( i > max_loop ) stop "max_loop < # of data"
99  mb=n
    nbk=k
    j=i

    rewind u
    do i=1,max_loop
       read(u,*,end=999) cbuf
       if ( cbuf == "Eigenvalues" ) exit
    end do
    read(u,*)
    do i=1,max_loop
       read(u,*,end=98) k,n,dummy(1:6)
    end do
    if ( i > max_loop ) stop "max_loop < # of data"
98  continue
    msp=1
    if ( i == j ) msp=2 

    write(*,*) "mb,nbk,msp=",mb,nbk,msp

    allocate( eval(mb,nbk,msp) ) ; eval=0.0d0
    allocate( occp(mb,nbk,msp) ) ; occp=0.0d0

    rewind u
    do i=1,max_loop
       read(u,*,end=999) cbuf
       if ( cbuf == "Eigenvalues" ) exit
    end do
    read(u,*)
    do k=1,nbk
       do n=1,mb
          read(u,*) i,j,(eval(n,k,s),dummy(1),occp(n,k,s),s=1,msp)
       end do
       do s=1,msp
          occp(2:mb,k,s)=occp(1,k,s)
       end do
    end do

    return

999 continue
    stop "The format of fort.99 is strange"

  END SUBROUTINE read_from_fort99

  SUBROUTINE dos2(s,e,f,N_MP)
    implicit none
    integer,intent(IN)  :: s
    real(8),intent(IN)  :: e
    real(8),intent(OUT) :: f
    integer,intent(IN)  :: N_MP
    real(8) :: c,x
    integer :: i,k

    f=0.0d0
    do k=1,nbk
       do i=1,mb

          x=e-eval(i,k,s)
          c=occp(i,k,s)/gamma
          f=f+c*gg0(N_MP,x/gamma)

       end do
    end do

  END SUBROUTINE dos2

  FUNCTION gg0(n,x)
    implicit none
    integer :: n,i
    real(8) :: x,gg,gg0,hp0,hp1,hp2,hp3
    gg0 = factor(0)*exp(-x*x)
    if ( n <= 0 ) return
    hp0 = 1.0d0
    hp1 = 2.0d0*x
    gg  = 0.0d0
    do i=1,n
       hp2 = 2.d0*x*hp1 - 2.d0*(2*i-1)*hp0
       hp3 = 2.d0*x*hp2 - 2.d0*(2*i  )*hp1
       gg  = gg + factor(i)*hp2
       hp0 = hp2
       hp1 = hp3
    end do
    gg0 = gg0 + gg*exp(-x*x)
    return
  END FUNCTION gg0

  SUBROUTINE init_gg0( nmax )
    implicit none
    integer,intent(IN) :: nmax
    integer :: n
    allocate( factor(0:nmax) ) ; factor=0.0d0
    factor(0)= 1.0d0/sqrt(acos(-1.0d0))
    if ( nmax == 0 ) return
    factor(1)=-1.0d0/(4.0d0*sqrt(acos(-1.0d0)))
    do n=2,nmax
       factor(n)=-factor(n-1)/(4.0d0*n)
    end do
  END SUBROUTINE init_gg0



END PROGRAM dos_simple
