!--------------------------------------------
! TM format
! J.-L. Martins' code ( latest version is atom-5.7.0 )
! The unit of energy is in Rydberg, and converted to Hartree in this routine
!--------------------------------------------
MODULE ps_read_TM_module

  use var_ps_member, only: ps1d, ps_allocate_ps1d

  implicit none

  PRIVATE
  PUBLIC :: ps_read_TM

CONTAINS

  SUBROUTINE ps_read_TM( g, ps )
    implicit none
    integer,intent(IN) :: g
    type(ps1d),intent(INOUT) :: ps
    character(2)  :: icorr,nameat
    character(3)  :: irel
    character(4)  :: nicore
    character(10) :: iray(6),ititle(7)
    integer :: i,j,norb,numu,nrr0,nrr,Lref
    integer,allocatable :: tlo(:),tinorm(:)
    integer,parameter :: lmax=3, nrmax=5000
    real(8) :: a0,b0,viou(lmax,nrmax),pi4,Zps
    real(8),parameter :: eps=1.d-7
    real(8),allocatable :: tviod(:,:),tanorm(:)

    write(*,'(a40," ps_read_tm")') repeat("-",40)

    pi4 = 4.0d0*acos(-1.0d0)

    allocate( tlo(2*lmax)       ) ; tlo=0
    allocate( tviod(lmax,nrmax) ) ; tviod=0.0d0
    allocate( tanorm(lmax)      ) ; tanorm=0.0d0
    allocate( tinorm(lmax)      ) ; tinorm=0

    read(g) nameat,icorr,irel,nicore,iray,ititle,norb,numu,nrr,a0,b0,Zps

    nrr=nrr+1

    ps%Mr = nrr
    ps%norb = norb-1
    call ps_allocate_ps1d( ps )

    ps%Zps = Zps

    read(g) ( ps%rad(i) ,i=2,nrr )

    do i=1,norb
       read(g) tlo(i),(tviod(i,j),j=2,nrr)
    end do
    do i=1,numu
       read(g) tlo(i+norb),(viou(tlo(i+norb)+1,j),j=2,nrr)
    end do

    read(g) ( ps%cdc(i),i=2,nrr)
    read(g) ( ps%cdd(i),i=2,nrr)
    read(g) ( ps%vql(i),i=2,nrr)

    read(g) norb

    do i=1,norb
       read(g) tinorm(tlo(i)+1),tanorm(tlo(i)+1)
    end do

    j=0
    do i=1,norb
       if ( tinorm(tlo(i)+1) /= 0 ) then
          j=j+1
          ps%lo(j)=tlo(i)
          ps%inorm(j)=tinorm(tlo(i)+1)
          ps%anorm(j)=tanorm(tlo(i)+1)
          ps%viod(2:nrr,j)=tviod(i,2:nrr)
       else
          Lref=tlo(i)
       end if
    end do
    norb=j

    deallocate( tinorm, tanorm, tviod, tlo )

! Values at the origin
!
    ps%rad(1)    = 0.0d0
    ps%vql(1)    = ps%vql(2)
    ps%viod(1,:) = ps%viod(2,:)

    do j=1,norb
       ps%viod(:,j) = ps%viod(:,j)*ps%rad(:)
    end do
    do i=2,nrr
       ps%cdc(i) = ps%cdc(i)*0.5d0/(pi4*ps%rad(i)**2)
    end do

    ps%cdc(1) = ps%cdc(2)
    ps%cdd(1) = ps%cdd(2)

! dr/dx
!
    ps%rab(:) = b0*( ps%rad(:) + a0 )

! Convert unit from Rydberg to Hartree
!
    ps%viod(:,:) = ps%viod(:,:)*sqrt(0.5d0)
    ps%anorm(:)  = ps%anorm(:)*sqrt(0.5d0)
    ps%vql(:)    = ps%vql(:)*0.5d0

! Cut-off radius
!
    do j=1,norb
       nrr0=1
       do i=nrr,1,-1
          if ( abs( ps%viod(i,j) ) > 1.d-2 ) then
             nrr0=i
             exit
          end if
       end do ! i
       do i=nrr0,nrr
          if ( abs( ps%viod(i,j) ) < eps ) then
             ps%Rps(j)  = ps%rad(i)
             ps%NRps(j) = i
             exit
          end if
       end do ! i
    end do ! j

    write(*,*) "*** TM Format ***"
    write(*,*) "Znuc=",ps%Zps
    write(*,*) "# of radial mesh points =",ps%Mr
    write(*,*) "# of orbitals =",ps%norb
    write(*,*) "angular momentum =",ps%lo(1:norb)
    write(*,*) "cut off radius =",ps%Rps(1:norb)
    write(*,*) "# of grid points within cut off radius",ps%NRps(1:norb)
    write(*,*) "uVu integral (anorm) ="
    write(*,'(1x,8f10.5)') ( ps%inorm(i)*ps%anorm(i),i=1,norb )
!    write(*,*) "Dij ="
!    write(*,'(1x,9f10.5)') (( ps_tm%Dij(i,j),i=1,norb ),j=1,norb)
    write(*,*) "sum(rhov)=",sum(ps%cdd*ps%rab)
    write(*,*) "sum(rhoc)=",sum(ps%cdc*ps%rab*(ps%rad)**2)*pi4

    write(*,'(a40," ps_read_tm(end)")') repeat("-",40)

  END SUBROUTINE ps_read_TM

END MODULE ps_read_TM_module
