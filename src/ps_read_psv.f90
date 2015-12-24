MODULE ps_read_PSV

  use var_ps_member, only: ps1d, file_ps, ippform, ps_allocate_ps1d
  use ps_read_PSV_g

  implicit none

  PRIVATE
  PUBLIC :: read_PSV

CONTAINS
  
  SUBROUTINE read_PSV( unit_ps, ielm, psp )

    implicit none
    integer,intent(IN) :: unit_ps,ielm
    type(ps1d),intent(INOUT) :: psp
    integer :: nsmpl,ndlc,ndlc_1,ndata,iost,iorb,jorb
    integer :: nl,l,i,j,m,ifchrg,verpot,ngauss
    integer,allocatable :: nr(:)
    real(8) :: temp,h,znuc,r1,r2,rend,dummy
    real(8) :: zatom,Rc_in,a1,a2,a3,a4
    real(8),allocatable :: rr(:),vl(:),cc(:),r(:)
    real(8),allocatable :: psi_(:,:,:),phi_(:,:,:),bet_(:,:,:)
    real(8),allocatable :: ddi_(:,:,:),qqr_(:,:,:)
    real(8),allocatable :: a0(:),b0(:),c0(:),cdd_coef_0(:,:,:)
    character(8) :: xc_pot
    character(18) :: inbuf18
    character(30) :: file_name
    character(80) :: inbuf
    logical :: multref, uspp

! Check pseudopotential type

    read(unit_ps,'(A)') inbuf

    if ( index(inbuf,'#PSV')==1 ) then

! Ver. 2.0 or more
!
       verpot=2
       read(unit_ps,*) xc_pot, zatom, ifchrg

! COMMENT LINE SKIP('#....')
!
       do i=1,10000
          read(unit_ps,'(A)') inbuf
          if ( index(inbuf,'#')==1 ) cycle
          backspace unit_ps
          exit
       end do
       if ( i>10000 ) stop "stop at psv"

       read(unit_ps,*) znuc

    else if ( index(inbuf,'#xPSV')==1 ) then
! xTAPP
!
       verpot=3
       ifchrg=0
       do i=1,10000
          read(unit_ps,'(A)') inbuf
          if ( index(inbuf,'#')==1 ) cycle
          backspace unit_ps
          exit
       end do
       if ( i>10000 ) stop "stop at psv"

       read(unit_ps,*) xc_pot
       read(unit_ps,*) znuc,zatom

    else

! Ver. 1.0 
!
       verpot=1
       ifchrg=0
       znuc=0
       write(6,'("WARNING! THE POTENTIAL TYPE IS OLD."/ &
                ,"XC-POTENTIAL TYPE: LDAPZ81 IS ASSUMED"/   &
                ,"We recommend rewrite or remake this potential data."/)')

       xc_pot='LDAPZ81'

! COMMENT LINE SKIP
!
       read(unit_ps,*)

    end if

    read(unit_ps,*) nl
    allocate( nr(nl) )
    read(unit_ps,*) nr(1:nl)
    read(unit_ps,*) nsmpl,Rc_in
    read(unit_ps,*) a1,a2,a3,a4

    temp=abs(a1+a3-1.0d0)
    if ( temp>1.0d-10 ) then
       write(6,*) ' LOCAL PARAMETER FITTING ERROR'
       write(6,*) 'A1,A2,A3,A4',a1,a2,a3,a4
       stop
    end if

! Read local potential
!
    read(unit_ps,*) ndlc

    allocate( rr(ndlc),vl(ndlc),cc(ndlc) )

    read(unit_ps,*,IOSTAT=iost) rr(1),vl(1),cc(1)
    backspace unit_ps

    read(unit_ps,*,IOSTAT=iost) rr(2),vl(2),cc(2)
    if ( rr(1)==rr(2) ) then
       iost=0
    else
       iost=1
       backspace unit_ps
       backspace unit_ps
    end if

    if ( iost==0 ) then
       do i=2,ndlc
          read(unit_ps,*) rr(i),vl(i),cc(i)
       end do
    else
       cc(:)=0.d0
       do i=2,ndlc
          read(unit_ps,*) rr(i),vl(i)
       end do
    end if

    if ( rr(1)<=0.d0 ) then
       write(6,*) ' ERROR IN LOACAL POT DATA rr(1)=',rr(1)
       stop 'stop at psv_format'
    end if

    h=exp( log(rr(ndlc)/rr(1))/(ndlc-1) )
    temp=h*rr(1)

    if ( abs(temp-rr(2))>1.d-12 ) then
       write(6,*) ' LOCAL POTENTIAL ERROR NDLC=',ndlc, 'H=',h
       write(6,*) ' rr(1),rr(2),r(3)=',rr(1),rr(2),rr(3)
       write(6,*) '          r(ndlc)=',rr(ndlc)
       write(6,*) ' temp=',temp,log(10.d0),log(exp(1.d0))
       stop
    end if

! Read nonlocal potential
!
    read(unit_ps,*)
    read(unit_ps,*) ndata,r1,rend,r2

    allocate( r(ndata) )

! MAKE R(I)
!
    if ( ndata/=nsmpl ) then
       write(6,*)' ERROR ON # OF DATA  NSMPL,NDATA=',nsmpl,ndata
       stop' MKQG0'
    end if
    h=exp( log(rend/r1)/(ndata-1) )
    r(1)=r1
    do i=2,ndata
       r(i)=r(i-1)*h
    end do
    if ( (abs(r(2)-r2)>1.d-6).or.(abs(rr(ndata)-rend)>1.d-8) ) then
       write(6,*)'ERROR  NDATA=',ndata
       write(6,'(" R(1),X1=",2e15.7," R(2),X2=",2e15.7/," R(NDATA),XEND=" &
                 ,2e15.7)') rr(1),r1,rr(2),r2,rr(ndata),rend
       stop
    end if
    temp=sqrt( sum( (rr(1:ndata)-r(1:ndata))**2 )/ndata )
    if ( temp>1.d-7 ) then
       write(*,*) "temp=",temp
       stop "Radial grid is inconsistent."
    end if

    m=maxval(nr(1:nl))
    allocate( psi_(1:nsmpl,m,nl) ) ; psi_=0.0d0
    allocate( phi_(1:nsmpl,m,nl) ) ; phi_=0.0d0
    allocate( bet_(1:nsmpl,m,nl) ) ; bet_=0.0d0
    allocate( ddi_(m,m,nl)       ) ; ddi_=0.0d0
    allocate( qqr_(m,m,nl)       ) ; qqr_=0.0d0

    do l=1,nl
       read(unit_ps,*)
       m=nr(l)
       read(unit_ps,'(3e20.12)') ( (ddi_(i,j,l),i=1,j), j=1,m )
       do j=1,m
       do i=j+1,m
          ddi_(i,j,l)=ddi_(j,i,l)
       end do
       end do
       read(unit_ps,*) ( (qqr_(i,j,l),i=1,j),j=1,m )
       do j=1,m
       do i=j+1,m
          qqr_(i,j,l)=qqr_(j,i,l)
       end do
       end do
       do j=1,m
          read(unit_ps,*)
          do i=1,nsmpl
             read(unit_ps,*) psi_(i,j,l),phi_(i,j,l),bet_(i,j,l)
          end do
       end do
    end do

    if ( verpot /= 3 ) then
       do j=1,10000
          read(unit_ps,'(A)') inbuf18
          if ( inbuf18=='### initial charge' ) then
             read(unit_ps,*) ngauss
             allocate( a0(ngauss),b0(ngauss),c0(ngauss) )
             do i=1,ngauss
                read(unit_ps,*) a0(i),b0(i),c0(i)
             end do
             exit
          end if
       end do
       if ( j>10000 ) stop "read_psv"
    end if

    psp%Mr = max( ndlc,nsmpl,ndata ) + 1
    psp%norb = max( 1, sum(nr(1:nl)) )
    call ps_allocate_ps1d( psp )

    psp%nlf      = nl
    psp%Zps      = znuc
    psp%Zelement = zatom

    iorb=0
    do l=1,nl
       psp%nrf(l) = nr(l)
       do j=1,nr(l)
          iorb=iorb+1
          psp%lo(iorb)    = l-1
          psp%no(iorb)    = j
          psp%anorm(iorb) = abs( ddi_(j,j,l) )
          psp%inorm(iorb) = sign( 1.0d0, ddi_(j,j,l) )
          psp%NRps(iorb)  = nsmpl+1
          psp%Rps(iorb)   = Rc_in
       end do
    end do

    do jorb=1,psp%norb
    do iorb=1,psp%norb
       if ( psp%lo(iorb) /= psp%lo(jorb) ) cycle
       l=psp%lo(iorb)
       i=psp%no(iorb)
       j=psp%no(jorb)
       psp%Dij(iorb,jorb) = ddi_(i,j,l+1)
    end do
    end do

    multref = .false.
    if ( count(psp%anorm/=0.0d0) < count(psp%Dij/=0.0d0) ) multref = .true.

    uspp = .false.
    if ( ippform(ielm) > 100 ) uspp = .true.

! r=0
    do i=1,ndlc
       psp%vql(i+1) = vl(i)
       psp%cdc(i+1) = cc(i)
       psp%rad(i+1) = rr(i)
    end do
    psp%rad(1) = 0.0d0
    psp%vql(1) = psp%vql(2)
!
    iorb=0
    do l=1,nl
       do j=1,nr(l)
          iorb=iorb+1
          temp=sqrt( psp%anorm(iorb) )
          if ( uspp .or. multref ) temp=1.0d0 ! uspp & multireference
          do i=1,nsmpl
             psp%viod(i+1,iorb) = bet_(i,j,l)*temp
          end do
          psp%viod(1,iorb) = 0.0d0
       end do ! j
    end do ! l

! dr/dx
    h=log( rr(ndlc)/rr(1) )/(ndlc-1)
    psp%rab(1) = 0.0d0
    do i=1,ndlc
       psp%rab(i+1) = rr(i)*h
    end do

    psp%parloc(1) = a1
    psp%parloc(2) = a2
    psp%parloc(3) = a3
    psp%parloc(4) = a4
!
! initial charge
!
    select case( verpot )
    case default

       temp=16.d0*atan(1.d0)

       psp%cdd(:)=0.0d0
       do l=1,ngauss
          do i=1,ndlc
             psp%cdd(i+1)=psp%cdd(i+1) &
                  +rr(i)**2*temp*(a0(l)+b0(l)*rr(i)**2)*exp(-c0(l)*rr(i)**2)
          end do
       end do

       allocate( psp%cdd_coef(3,ngauss) )
       psp%cdd_coef=0.0d0
       psp%cdd_coef(1,1:ngauss) = a0(1:ngauss)
       psp%cdd_coef(2,1:ngauss) = b0(1:ngauss)
       psp%cdd_coef(3,1:ngauss) = c0(1:ngauss)
       psp%ngauss = ngauss

    case( 3 )

       file_name=trim(file_ps(ielm))//".ichr"
       open(unit_ps+1,file=file_name,status='old')
       read(unit_ps+1,*) ndlc_1
       if ( ndlc /= ndlc_1 ) stop "ndlc/=ndlc_1!!!"

       temp=16.d0*atan(1.d0)

       psp%cdd(:)=0.0d0
       do i=1,ndlc
          read(unit_ps+1,*) dummy, r1
          psp%cdd(i+1)=temp*r1*rr(i)**2
       end do
       close(unit_ps+1)

    end select

    write(*,*) "*** PSV format ***"
    write(*,*) "zatom                   =",psp%Zelement
    write(*,*) "# of radial mesh points =",psp%Mr
    write(*,*) "r(1),...,r(end)         =",psp%rad(1),psp%rad(psp%Mr)
    write(*,*) "NRps                    =",psp%NRps(1:psp%norb)
    write(*,*) "# of orbitals           =",psp%norb
    write(*,*) "angular momentum        =",psp%lo(1:psp%norb)
    write(*,*) "cut off radius          =",psp%Rps(1:psp%norb)
    write(*,*) "Zps                     =",psp%Zps
    write(*,*) "xc_pot                  =   ",xc_pot
    write(*,*) "a1,a2,a3,a4             =",a1,a2
    write(*,*) "                         ",a3,a4
    write(*,*) "Reference               =",psp%nrf(1:psp%nlf)
    write(*,*) "anorm                   =",psp%anorm(1:psp%norb)
    write(*,*) "inorm                   =",psp%inorm(1:psp%norb)
    write(*,*) "sum(cdd)                =",sum( psp%cdd(:)*psp%rab(:) )
    write(*,*) "sum(cdc)                =",sum( psp%cdc(:)*psp%rab(:) )*temp

    if ( uspp ) call readPSVG( unit_ps,ddi_,qqr_,psi_,phi_,bet_,psp )

    if ( verpot /= 3 ) deallocate( c0,b0,a0 )
    deallocate( qqr_ )
    deallocate( ddi_ )
    deallocate( bet_ )
    deallocate( phi_ )
    deallocate( psi_ )
    deallocate( r )
    deallocate( cc,vl,rr )
    deallocate( nr )

  END SUBROUTINE read_PSV

END MODULE ps_read_PSV
