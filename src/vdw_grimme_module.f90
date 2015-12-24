!
! --- Semiempirical van der Waals functional ---
!
! REFERENCE
! Stefan Grimme, "Semiempirical GGA-Type Density Functional Constructed
! with a Long-Range Dispersion Correction",
! Journal of Computational Chemistry 27,1787 (2006).
!
MODULE vdw_grimme_module

  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: init_vdw_grimme, calc_E_vdw_grimme, calc_F_vdw_grimme &
       ,read_vdw_grimme, get_E_vdw_grimme

  real(8) :: C6(54),R0(54),s6

  real(8),parameter :: d        = 20.0d0
  real(8),parameter :: Rcut     = 30.0d0/0.529177d0
  real(8),parameter :: s6_b97   = 1.25d0
  real(8),parameter :: s6_pbe   = 0.75d0
  real(8),parameter :: s6_blyp  = 1.20d0
  real(8),parameter :: s6_bp86  = 1.05d0
  real(8),parameter :: s6_tpss  = 1.00d0
  real(8),parameter :: s6_b3lyp = 1.05d0

  real(8) :: Edisp0 = 0.0d0
  real(8) :: aa(3,3)
  integer :: MI_0, MI_1
  integer :: lmax
  integer,allocatable :: atom_num(:), Kion(:)
  logical :: flag_init = .true.
  logical :: iswitch_vdw = .false.
  character(8) :: xctype

  include 'mpif.h'

CONTAINS


  SUBROUTINE read_vdw_grimme
    implicit none
    call write_border( 0, " read_vdw_grimme(start)" )
    call IOTools_findKeyword( "VDW", iswitch_vdw, flag_bcast=.true. )
    call IOTools_readStringKeyword( "XCTYPE", xctype )
    call write_border( 0, " read_vdw_grimme(end)" )
  END SUBROUTINE read_vdw_grimme


  SUBROUTINE init_vdw_grimme( aa_in, ki_in, z_in )
    implicit none
    real(8),intent(IN) :: aa_in(3,3)
    integer,intent(IN) :: ki_in(:), z_in(:)
    integer :: Natom

    if ( .not.iswitch_vdw ) return

    call write_border( 0, " init_vdw_grimme(start)" )

    aa(:,:) = aa_in(:,:)

    Natom = size(ki_in)

    allocate( Kion(Natom) ) ; Kion=0
    allocate( atom_num(size(z_in)) ) ; atom_num=0
    Kion(:) = ki_in(:)
    atom_num(:) = z_in(:)

    if ( any( atom_num == 0 ) ) then
       write(*,*) "atomic number is necessary for vdW"
       stop "stop@init_vdw_grimme"
    end if

    if ( xctype(1:8) == "GGAPBE96" ) then
       s6 = s6_pbe
    else
       s6 = s6_B97
    end if

    call get_lattice_max( Rcut, aa, lmax )
    call prep_parallel( Natom, MI_0, MI_1 )

!-----------------------------------------
! C6 parameters
!
    C6( 1: 6) = (/ 0.14d0, 0.08d0, 1.61d0, 1.61d0, 3.13d0, 1.75d0/)
    C6( 7:12) = (/ 1.23d0, 0.70d0, 0.75d0, 0.63d0, 5.71d0, 5.71d0/)
    C6(13:18) = (/10.79d0, 9.23d0, 7.84d0, 5.57d0, 5.07d0, 4.61d0/)
    C6(19:20) = (/10.80d0,10.80d0 /)
    C6(21:30) = 10.80d0
    C6(31:36) = (/16.99d0,17.10d0,16.37d0,12.64d0,12.47d0,12.01d0/)
    C6(37:38) = (/24.67d0,24.67d0/)
    C6(39:48) = 24.67d0
    C6(49:54) = (/37.32d0,38.71d0,38.44d0,31.74d0,31.50d0,29.99d0/)
!
! van der Waals Radii
!
    R0( 1: 6) = (/ 1.001d0, 1.012d0, 0.825d0, 1.408d0, 1.485d0, 1.452d0 /)
    R0( 7:12) = (/ 1.397d0, 1.342d0, 1.287d0, 1.243d0, 1.144d0, 1.364d0 /)
    R0(13:18) = (/ 1.639d0, 1.716d0, 1.705d0, 1.683d0, 1.639d0, 1.595d0 /)
    R0(19:20) = (/ 1.485d0, 1.474d0 /)
    R0(21:30) = 1.562d0
    R0(31:36) = (/ 1.650d0, 1.727d0, 1.760d0, 1.771d0, 1.749d0, 1.727d0 /)
    R0(37:38) = (/ 1.628d0, 1.606d0 /)
    R0(39:48) = 1.639d0
    R0(49:54) = (/ 1.672d0, 1.804d0, 1.881d0, 1.892d0, 1.892d0, 1.881d0 /)

    C6(:)=17.338d0*C6(:)   !( J/mol*nm^6 -> a.u. )
    R0(:)=R0(:)/0.529177d0 !( \AA -> a.u. )
!
!-----------------------------------------

    flag_init = .false.

    call write_border( 0, " init_vdw_grimme(end)" )

  END SUBROUTINE init_vdw_grimme


  SUBROUTINE calc_E_vdw_grimme( asi, Edisp_out )
    implicit none
    real(8),intent(IN)  :: asi(:,:)
    real(8),optional,intent(OUT) :: Edisp_out
    integer :: i,j,m1,m2,m3,l,ierr,MI
    real(8) :: Edisp
    real(8) :: tx,ty,tz,err,err0,chk_sum
    real(8) :: Rix,Riy,Riz,Rjx,Rjy,Rjz,Rij
    real(8) :: C6i,C6j,C6ij,R0i,R0j,fdmp,Rr
    real(8),save :: chk_sum0=1.d10

    if ( .not.iswitch_vdw ) return

    call write_border( 1, " calc_vdw_grimme(start)" )

    if ( flag_init ) then
       write(*,*) "You should call init_vdw_grimme first"
       stop "stop@calc_Edisp_vdw_grimme"
    end if

!    chk_sum  = sum( asi**2 )
!    err      = abs( chk_sum - chk_sum0 )
!    chk_sum0 = chk_sum
!    if ( err < 1.d-14 ) then
!       Edisp = Edisp0
!       return
!    end if

    Edisp  = 0.0d0
!    Edisp0 = 1.0d10
    Edisp0 = 0.0d0

    MI = size( asi, 2 )

    do l=lmax,lmax

       do m3=-l,l
       do m2=-l,l
       do m1=-l,l

!        if ( abs(m1)<l .and. abs(m2)<l .and. abs(m3)<l ) cycle

          tx = aa(1,1)*m1 + aa(1,2)*m2 + aa(1,3)*m3
          ty = aa(2,1)*m1 + aa(2,2)*m2 + aa(2,3)*m3
          tz = aa(3,1)*m1 + aa(3,2)*m2 + aa(3,3)*m3

          do i=MI_0,MI_1

             Rix = aa(1,1)*asi(1,i) + aa(1,2)*asi(2,i) + aa(1,3)*asi(3,i) + tx
             Riy = aa(2,1)*asi(1,i) + aa(2,2)*asi(2,i) + aa(2,3)*asi(3,i) + ty
             Riz = aa(3,1)*asi(1,i) + aa(3,2)*asi(2,i) + aa(3,3)*asi(3,i) + tz

             C6i = C6( atom_num(Kion(i)) )
             R0i = R0( atom_num(Kion(i)) )

             do j=1,MI

                Rjx = aa(1,1)*asi(1,j) + aa(1,2)*asi(2,j) + aa(1,3)*asi(3,j)
                Rjy = aa(2,1)*asi(1,j) + aa(2,2)*asi(2,j) + aa(2,3)*asi(3,j)
                Rjz = aa(3,1)*asi(1,j) + aa(3,2)*asi(2,j) + aa(3,3)*asi(3,j)

                Rij = sqrt( (Rix-Rjx)**2 + (Riy-Rjy)**2 + (Riz-Rjz)**2 )

                if ( Rij < 1.d-10 ) cycle

                C6j = C6( atom_num(Kion(j)) )
                R0j = R0( atom_num(Kion(j)) )

                C6ij = sqrt( C6i*C6j )

                Rr = R0i + R0j

                fdmp =1.d0/( 1.d0 + exp( -d*(Rij/Rr-1.d0) ) )

                Edisp = Edisp + C6ij/Rij**6 * fdmp

             end do ! j

          end do ! i

          err0   = abs( Edisp - Edisp0 )
          Edisp0 = Edisp

       end do ! m3
       end do ! m2
       end do ! m1

    end do ! l

    call mpi_allreduce(Edisp0,Edisp,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(err0,err,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

    Edisp  = -s6*Edisp/2.0d0
    Edisp0 = Edisp

    if ( present(Edisp_out) ) Edisp_out = Edisp

    if ( err > 1.d-6 ) then
       write(*,*) "Edisp was not converged"
       write(*,*) "Edisp,err =",Edisp,err
    end if

    call write_border( 1, " calc_vdw_grimme(end)" )

  END SUBROUTINE calc_E_vdw_grimme


  SUBROUTINE get_E_vdw_grimme( Edisp_out )
    implicit none
    real(8),intent(OUT) :: Edisp_out
    Edisp_out = Edisp0
  END SUBROUTINE get_E_vdw_grimme


  SUBROUTINE calc_F_vdw_grimme( asi, Fdisp )
    implicit none
    real(8),intent(IN)  :: asi(:,:)
    real(8),intent(INOUT) :: Fdisp(:,:)
    integer :: i,j,m1,m2,m3,l,MI
    real(8) :: Rix,Riy,Riz,Rjx,Rjy,Rjz,Rij
    real(8) :: C6i,C6j,C6ij,R0i,R0j,fdmp,Rr
    real(8) :: tx,ty,tz,err,err0,c0,c1,ierr
    real(8),allocatable :: Fdisp0(:,:), Fdisp1(:,:)

    if ( .not.iswitch_vdw ) return

    if ( flag_init ) then
       write(*,*) "You should call init_vdw_grimme first"
       stop "stop@calc_F_vdw_grimme"
    end if

!    Fdisp(:,:) = 0.d0

    MI = size( asi, 2 )

    allocate( Fdisp0(3,MI) ) ; Fdisp0=1.0d10
    allocate( Fdisp1(3,MI) ) ; Fdisp1=0.0d0

    c0 = s6*6.d0
    c1 = s6*d

    do l=lmax,lmax

       do m3=-l,l
       do m2=-l,l
       do m1=-l,l

!        if ( abs(m1)<l .and. abs(m2)<l .and. abs(m3)<l ) cycle

          tx = aa(1,1)*m1 + aa(1,2)*m2 + aa(1,3)*m3
          ty = aa(2,1)*m1 + aa(2,2)*m2 + aa(2,3)*m3
          tz = aa(3,1)*m1 + aa(3,2)*m2 + aa(3,3)*m3

          do i=MI_0,MI_1

             Rix = aa(1,1)*asi(1,i) + aa(1,2)*asi(2,i) + aa(1,3)*asi(3,i) + tx
             Riy = aa(2,1)*asi(1,i) + aa(2,2)*asi(2,i) + aa(2,3)*asi(3,i) + ty
             Riz = aa(3,1)*asi(1,i) + aa(3,2)*asi(2,i) + aa(3,3)*asi(3,i) + tz

             C6i = C6( atom_num(Kion(i)) )
             R0i = R0( atom_num(Kion(i)) )

             do j=1,MI

                Rjx = aa(1,1)*asi(1,j) + aa(1,2)*asi(2,j) + aa(1,3)*asi(3,j)
                Rjy = aa(2,1)*asi(1,j) + aa(2,2)*asi(2,j) + aa(2,3)*asi(3,j)
                Rjz = aa(3,1)*asi(1,j) + aa(3,2)*asi(2,j) + aa(3,3)*asi(3,j)

                Rij = sqrt( (Rix-Rjx)**2 + (Riy-Rjy)**2 + (Riz-Rjz)**2 )

                if ( Rij < 1.d-10 ) cycle

                C6j = C6( atom_num(Kion(j)) )
                R0j = R0( atom_num(Kion(j)) )

                C6ij = sqrt( C6i*C6j )

                Rr = R0i + R0j

                fdmp =1.d0/( 1.d0 + exp( -d*(Rij/Rr-1.d0) ) )

                Fdisp1(1,j)=Fdisp1(1,j)-c0*C6ij*fdmp*(Rix-Rjx)/Rij**8
                Fdisp1(2,j)=Fdisp1(2,j)-c0*C6ij*fdmp*(Riy-Rjy)/Rij**8
                Fdisp1(3,j)=Fdisp1(3,j)-c0*C6ij*fdmp*(Riz-Rjz)/Rij**8

                Fdisp1(1,j)=Fdisp1(1,j)-c1*C6ij/Rr*fdmp*(1-fdmp)*(Rix-Rjx)/Rij**7
                Fdisp1(2,j)=Fdisp1(2,j)-c1*C6ij/Rr*fdmp*(1-fdmp)*(Riy-Rjy)/Rij**7
                Fdisp1(3,j)=Fdisp1(3,j)-c1*C6ij/Rr*fdmp*(1-fdmp)*(Riz-Rjz)/Rij**7

             end do ! j

          end do ! i

          err0        = sum( (Fdisp1-Fdisp0)**2 )
          Fdisp0(:,:) = Fdisp1(:,:)

       end do ! m3
       end do ! m2
       end do ! m1

    end do ! l

    call mpi_allreduce(Fdisp0,Fdisp1,3*MI,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(err0,err,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

    Fdisp(:,:) = Fdisp(:,:) + Fdisp1(:,:)

    deallocate( Fdisp1 )
    deallocate( Fdisp0 )

    if ( err > 1.0d-10 ) then
       write(*,*) "Fdisp was not converged"
       write(*,*) "err =",err
    end if

  END SUBROUTINE calc_F_vdw_grimme


  SUBROUTINE get_lattice_max(rcut,aa,lmax)
    implicit none
    real(8),intent(IN)  :: rcut,aa(3,3)
    integer,intent(OUT) :: lmax
    real(8) :: taaaa(3,3),taa(3,3),x(3,3),y(3,3),e0(3),e(3),c,r2
    integer :: loop,i,j,k,imax,jmax,kmax,m

    taa(:,:) = transpose( aa(:,:) )
    taaaa(:,:) = matmul( taa(:,:),aa(:,:) )

!- solve eigenvalue problem by power method

    call random_number(x)
    e=0.d0
    do loop=1,1000
       y(:,:) = matmul( taaaa(:,:),x(:,:) )
       e0(:) = e(:)
       do i=1,3
          e(i) = sum( x(:,i)*y(:,i) )
       end do
       do i=1,3
          do j=1,i-1
             c=sum(y(:,j)*y(:,i))
             y(:,i)=y(:,i)-c*y(:,j)
          end do
          c = 1.d0/sqrt( sum(y(:,i)*y(:,i)) )
          y(:,i) = c*y(:,i)
       end do
       x(:,:)=y(:,:)
    end do
    c=sum(abs(e-e0)/e)
    if ( c > 1.d-12 ) then
       write(*,*) "e,e0,e-e0=",e,e0,e-e0
       stop "not converged (get_rcut)"
    end if

    r2   = rcut*rcut
    c    = sqrt(r2/minval(e))
    lmax = int(c)+1

    write(*,*) "Rcut=",Rcut
    write(*,*) "e=",e(1:3)
    write(*,*) "c=",c
    write(*,*) "lmax=",lmax

    m=0
    imax=0
    jmax=0
    kmax=0
    do k=-lmax,lmax
    do j=-lmax,lmax
    do i=-lmax,lmax
       x(1,1) = i*aa(1,1) + j*aa(1,2) + k*aa(1,3)
       x(2,1) = i*aa(2,1) + j*aa(2,2) + k*aa(2,3)
       x(3,1) = i*aa(3,1) + j*aa(3,2) + k*aa(3,3)
       c = sum( x(1:3,1)**2 )
       if ( c <= r2 ) then
          m=m+1
          imax = max( imax, abs(i) )
          jmax = max( jmax, abs(i) )
          kmax = max( kmax, abs(i) )
       end if
    end do
    end do
    end do

    write(*,*) "m=",m
    write(*,*) "imax,jmax,kmax=",imax,jmax,kmax

  END SUBROUTINE get_lattice_max


  SUBROUTINE prep_parallel( n, n0,n1 )
    implicit none
    integer,intent(IN)  :: n
    integer,intent(OUT) :: n0,n1
    integer :: i,j,np,myrnk
    integer,allocatable :: ircnt(:)

    call MPI_COMM_SIZE( MPI_COMM_WORLD, np   , i )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myrnk, i )

    allocate( ircnt(0:np-1) )

    ircnt(:)=0
    do i=1,n
       j=mod(i-1,np)
       ircnt(j)=ircnt(j)+1
    end do

    n0 = sum( ircnt(0:myrnk) ) - ircnt(myrnk) + 1
    n1 = n0 + ircnt(myrnk) - 1

    deallocate( ircnt )

  END SUBROUTINE prep_parallel


END MODULE vdw_grimme_module
