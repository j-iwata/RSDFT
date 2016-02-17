MODULE bz_module

  use io_tools_module
  use symmetry_module, only: isymmetry, nsym, rgb

  implicit none

  PRIVATE
  PUBLIC :: generate_bz
  PUBLIC :: read_bz

  integer,PUBLIC :: nk
  integer,PUBLIC :: mmm(3,2)
  integer,PUBLIC :: Nbzsm
  real(8),allocatable,PUBLIC :: kbb(:,:)
  real(8),allocatable,PUBLIC :: weight_bz(:)

  integer,PUBLIC :: MMBZ
  integer :: npbz
  real(8),allocatable :: wbz(:)
  integer :: ndata_read_k=0
  real(8) :: kbb0(3)
  data kbb0/0.d0,0.d0,0.d0/
  integer :: use_inversion=1

CONTAINS


  SUBROUTINE read_bz
    implicit none
    nk=2
    mmm(1:3,1)=(/ 2,2,2 /)
    mmm(1:3,2)=(/ 2,2,2 /)
    ndata_read_k=0
    kbb0(1:3)=0.0d0
    npbz=0
    use_inversion=1
    call IOTools_readIntegerKeyword( "NK", nk )
    call IOTools_readIntegerKeyword( "MMM1", mmm(:,1) )
    call IOTools_readIntegerKeyword( "MMM2", mmm(:,2) )
    call IOTools_readIntegerKeyword( "NPBZ", npbz )
    call IOTools_readIntegerKeyword( "INVBZ", use_inversion )
  END SUBROUTINE read_bz


  SUBROUTINE generate_bz
    implicit none
    logical :: disp_switch
    integer :: i,k,k1,iw,iw0,iw1,m1,m2,m3,mm1,mm2,mm3,i1,i2,i3,p1(3),p2(3)
    integer,allocatable :: mm(:,:),m(:,:),w(:)

    if ( isymmetry > 0 ) then
       call generate_bz_sym
       return
    end if

    call write_border( 0," generate_bz(start)" )

    m1 =mmm(1,1) ; m2 =mmm(2,1) ; m3 =mmm(3,1)
    mm1=mmm(1,2) ; mm2=mmm(2,2) ; mm3=mmm(3,2)

    k=(2*m1+1)*(2*m2+1)*(2*m3+1)*2
    allocate( mm(3,k),m(3,k),w(k) )
    mm=0 ; m=0 ; w=0
    k=0 ; k1=0

    iw1= 1
    iw0=-1 ; if ( use_inversion < 1 ) iw0=1

    do i1=-m1,m1,mm1
    do i2=-m2,m2,mm2
    loop_A : do i3=-m3,m3,mm3

       do iw=iw1,iw0,-2

          p1(1)=i1*iw ; p1(2)=i2*iw ; p1(3)=i3*iw ; p2(1:3)=p1(1:3)
          do i=1,3
             p1(i)=mod(p2(i),nk)
             if ( p1(i)>  nk/2 ) p1(i)=p1(i)-nk
             if ( p1(i)<=-nk/2 ) p1(i)=p1(i)+nk
          end do
          if ( k1>0 ) then
             do i=1,k1
                if ( mm(1,i)==p1(1) .and. &
                     mm(2,i)==p1(2) .and. &
                     mm(3,i)==p1(3)         ) cycle loop_A
             end do
          end if

          if ( iw==1 ) then
             k=k+1
             m(1:3,k)=p1(1:3)
             w(k)=1
          else
             w(k)=2
          end if

          k1=k1+1
          mm(1:3,k1)=p1(1:3)

       end do ! iw

    end do loop_A
    end do
    end do

    Nbzsm = k
    if ( npbz > k ) Nbzsm=npbz

    MMBZ  = k1

    allocate( weight_bz(Nbzsm) ) ; weight_bz=0.d0
    allocate( kbb(3,Nbzsm)     ) ; kbb=0.d0
    allocate( wbz(Nbzsm)       ) ; wbz=0.d0

    kbb(1:3,1:k)=real(m(1:3,1:k),8)/nk
    weight_bz(1:k)=real(w(1:k),8)/MMBZ
    wbz(1:Nbzsm) = weight_bz(1:Nbzsm)

    if ( Nbzsm == ndata_read_k ) then
       kbb(1:3,1) = kbb0(1:3)
       ndata_read_k = 0
    end if

    deallocate( w,m,mm )

    call check_disp_switch( disp_switch, 0 )
    if ( disp_switch ) then
       write(*,*) "Nbzsm, MMBZ =",Nbzsm,MMBZ
       write(*,'(1x,a4,a30,a12)') "","kbb","weight_bz"
       do k=1,Nbzsm
          write(*,'(1x,i4,3f10.5,f12.5)') k,kbb(:,k),weight_bz(k)
       end do
    end if

    call write_border( 0," generate_bz(end)" )

  END SUBROUTINE generate_bz


  SUBROUTINE generate_bz_sym
    implicit none
    logical :: disp_switch
    integer,allocatable :: mm(:,:),m(:,:),w(:),w1(:),w2(:)
    integer :: m1,m2,m3,mm1,mm2,mm3,i1,i2,i3,p1(3),p2(3)
    integer :: i,k,k1,iw,iw0,iw1,nkmax,p3(3),ns,ni,is,nni,ig
    real(8) :: c,tmp(3)

    call write_border( 0," generate_bz_sym(start)" )

    m1 =mmm(1,1) ; m2 =mmm(2,1) ; m3 =mmm(3,1)
    mm1=mmm(1,2) ; mm2=mmm(2,2) ; mm3=mmm(3,2)

    nkmax=(2*m1+1)*(2*m2+1)*(2*m3+1)*2
    allocate( mm(3,nkmax),m(3,nkmax),w(nkmax) ) ; mm=0 ; m=0 ; w=0
    allocate( w1(nkmax),w2(nkmax) ) ; w1=0 ; w2=0

    iw1= 1
    iw0=-1 ; if ( use_inversion < 1 ) iw0=1

    ni=1
    is=1

    do i1=-m1,m1,mm1
    do i2=-m2,m2,mm2
       loop_3 : do i3=-m3,m3,mm3

          p1(1)=i1 ; p1(2)=i2 ; p1(3)=i3 ; p2(1:3)=p1(1:3)

          do i=1,3
             p1(i)=mod(p2(i),nk)
             if ( p1(i) >   nk/2 ) p1(i)=p1(i)-nk
             if ( p1(i) <= -nk/2 ) p1(i)=p1(i)+nk
          end do

          do i=1,ni-1
             if ( all(p1(1:3)==mm(1:3,i)) ) cycle loop_3
          end do

          ns =0
          nni=ni

          do iw=iw1,iw0,-2
             loop_sym : do ig=1,nsym

                if ( ni>nkmax ) stop "generate_bz_sym(1)"

                tmp(:) = matmul( rgb(:,:,ig),p1(:) )*iw
                p3(:) = nint( tmp(:) )

                do i=1,3
                   p2(i)=mod(p3(i),nk)
                   if ( p2(i) >  nk/2 ) p2(i)=p2(i)-nk
                   if ( p2(i) <=-nk/2 ) p2(i)=p2(i)+nk
                end do

                do i=nni,ni-1
                   if ( all(p2(:)==mm(:,i)) ) cycle loop_sym
                end do
                ns=ns+1
                mm(:,ni)=p2(:)
                ni=ni+1

             end do loop_sym
          end do ! iw

          w(is)=ns

          m(1:3,is)=mm(1:3,nni)

          is=is+1

       end do loop_3
    end do ! i2
    end do ! i1

    is=is-1
    ni=ni-1

    do k=1,ni
       w1(k)=1
       do i=1,3
          i1=mod(i,3)+1
          i2=mod(i+1,3)+1
          if ( abs(mm(i,k)) == nk/2 ) then
             do k1=1,ni
                if ( k == k1 ) cycle
                if ( mm(i ,k1) ==-mm(i ,k) .and. &
                     mm(i1,k1) == mm(i1,k) .and. &
                     mm(i2,k1) == mm(i2,k) ) w1(k)=w1(k)*2
             end do
          end if
       end do ! i
    end do ! k

    do k=1,is
       do k1=1,ni
          if ( all(m(1:3,k)==mm(1:3,k1)) ) then
             w2(k)=w1(k1)
             exit
          end if
       end do
    end do

    Nbzsm = is
    if ( npbz > is ) Nbzsm = npbz

    MMBZ = ni

    allocate( weight_bz(Nbzsm) ) ; weight_bz=0.0d0
    allocate( kbb(3,Nbzsm)     ) ; kbb=0.0d0
    allocate( wbz(Nbzsm)       ) ; wbz=0.0d0

    do k=1,Nbzsm
       kbb(1:3,k)   = real( m(1:3,k), 8 )/real( nk, 8 )
       weight_bz(k) = real( w(k), 8 )/real( w2(k), 8 )
    end do

    c=sum( weight_bz(:) )
    weight_bz(:)=weight_bz(:)/c

    wbz(1:Nbzsm) = weight_bz(1:Nbzsm)

    call check_disp_switch( disp_switch, 0 )
    if ( disp_switch ) then
       write(*,*) "Nbzsm, MMBZ =",Nbzsm,MMBZ
       write(*,*) "KBB"
       do k=1,Nbzsm
          write(*,'(1x,i4,3f10.5,f12.5,2i5)') k,kbb(:,k),wbz(k),w(k),w2(k)
       end do
       write(*,*) "sum(w)  =",sum(w),sum(w2),sum(w1)
       write(*,*) "sum(wbz)=",sum(wbz)
    end if

    deallocate( w2,w1,w,m,mm )

    call write_border( 0," generate_bz_sym(end)" )

    return

  END SUBROUTINE generate_bz_sym


END MODULE bz_module
