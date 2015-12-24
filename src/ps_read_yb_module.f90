MODULE ps_read_YB_module

  use var_ps_member, only: ps1d, ps_allocate_ps1d

  implicit none

  PRIVATE
  PUBLIC :: ps_read_YB

CONTAINS

  SUBROUTINE ps_read_YB( g, ps )
    implicit none
    integer,intent(IN) :: g
    type(ps1d),intent(INOUT) :: ps
    real(8),parameter :: E2=14.39965d0, H2M=3.80998d0
    real(8),parameter :: a_B=0.529177249d0
    real(8),parameter :: Ry=13.6056981d0, Pi=3.141592653589793d0
    integer :: i,j,L,Mlps0,Lref,nrr,norb
    real(8) :: rPC,dr,c,znuc_tmp,Zps
    real(8),allocatable :: vtmp(:,:),utmp(:,:),rtmp(:)

! Read

     read(g,*) nrr, dr, Mlps0, Zps

     nrr = nrr + 1
     norb = Mlps0 + 1

     allocate( vtmp(nrr,norb) ) ; vtmp=0.0d0
     allocate( utmp(nrr,norb) ) ; utmp=0.0d0
     allocate( rtmp(norb)     ) ; rtmp=0.0d0

     ps%Mr = nrr
     ps%norb = norb - 1

     call ps_allocate_ps1d( ps )

     ps%Zps = Zps

     read(g,*) rPC, ( rtmp(L), L=1,Mlps0+1 )

     do i=1,nrr
        read(g,*) ps%rad(i),ps%cdc(i),( vtmp(i,L), L=1,Mlps0+1 )
     end do
     do i=1,nrr
        read(g,*) ps%rad(i),( utmp(i,L), L=1,Mlps0+1 )
     end do

     ps%rab(1:nrr) = dr

     Lref = Mlps0

     i=0
     do L=0,Mlps0
        if ( L == Lref ) then
           ps%vql(1:nrr) = vtmp(1:nrr,L+1)
           cycle
        end if
        i=i+1
        ps%lo(i)=L
        ps%Rps(i)=rtmp(L+1)
        ps%viod(1:nrr,i) = vtmp(1:nrr,L+1)
        ps%ups(1:nrr,i) = utmp(1:nrr,L+1)
     end do
     norb=i

     do j=1,norb
        do i=1,nrr
           if ( ps%rad(i) >= ps%Rps(j) ) then
              ps%NRps(j)=i
              exit
           end if
        end do
        ps%Rps(j) = ps%rad( ps%NRps(j) )
     end do

     do j=1,norb
        ps%viod(1:nrr,j)=( ps%viod(1:nrr,j) - ps%vql(1:nrr) )*ps%ups(1:nrr,j)
        ps%anorm(j) = sum( ps%viod(1:nrr,j)*ps%ups(1:nrr,j) )*dr
        ps%inorm(j) = sign( 1.0d0, ps%anorm(j) )
        ps%anorm(j) = abs( ps%anorm(j) )
        ps%viod(1:nrr,j) = ps%viod(1:nrr,j)/sqrt( ps%anorm(j) )
     end do

     ps%rad(:)    = ps%rad(:)/a_B
     ps%rab(:)    = ps%rab(:)/a_B
     ps%Rps(:)    = ps%Rps(:)/a_B
     ps%vql(:)    = ps%vql(:)/(2.0d0*Ry)
     ps%viod(:,:) = ps%viod(:,:)/sqrt(2.0d0*Ry/a_B)
     ps%anorm(:)  = ps%anorm(:)/(2.0d0*Ry)

     ps%cdd(:)=0.0d0
     znuc_tmp=0.0d0
     do j=1,norb
        c=2.0d0*(2*ps%lo(j)+1)
        znuc_tmp=znuc_tmp+c
        if ( znuc_tmp <= ps%Zps ) then
           ps%cdd(1:nrr) = ps%cdd(1:nrr) + c*ps%ups(1:nrr,j)**2
        end if
     end do

     write(*,*) "*** KY format ***"
     write(*,*) "Znuc=",ps%Zps
     write(*,*) "# of radial mesh points =",nrr
     write(*,*) "# of orbitals =",ps%norb
     write(*,*) "angular momentum =",ps%lo(1:norb)
     write(*,*) "reference L =",Lref
     write(*,*) "cut off radius =",ps%Rps(1:norb)
     write(*,*) "uVu integral (anorm) ="
     write(*,'(1x,8f10.5)') ( ps%inorm(i)*ps%anorm(i),i=1,norb )
     write(*,*) "sum(rhov)=",sum(ps%cdd*ps%rab)

     deallocate( rtmp )
     deallocate( utmp )
     deallocate( vtmp )

   END SUBROUTINE ps_read_YB

 END MODULE ps_read_YB_module
