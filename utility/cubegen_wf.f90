PROGRAM cubegen_wf
  implicit none

  integer,parameter :: u1=1,u2=970,u3=5,iu=2
  integer :: ML,ML1,ML2,ML3,zatom(9)
  integer :: i1,i2,i3,i,n,k,s,MI,MKI
  integer,allocatable :: LL(:,:),Kion(:)
  real(8),allocatable :: rtmp(:),rho(:,:,:),occ(:,:,:)
  real(8),allocatable :: asi(:,:),rsi(:,:)
  real(8) :: ax,aa(3,3)
  character(30) :: cbuf,file_name
  complex(8),allocatable :: ztmp(:)

  integer :: MB,MB1,MB2
  integer :: MBZ,MSP
  integer :: ndata, iflag_abs
  integer,allocatable :: idata(:,:)
  logical :: flag_real8

! ---

  read(u3,*) cbuf, iflag_abs
  read(u3,*) file_name
  read(u3,*) MBZ
  read(u3,*) MSP
  read(u3,*) ndata

  if ( cbuf == "r" .or. cbuf == "R" ) then
     flag_real8 = .true.
  else if ( cbuf == "c" .or. cbuf == "C" ) then
     flag_real8 = .false.
     if ( iflag_abs == 0 ) then
        stop "complex wave function should be converted to a real quantity"
     end if
  end if

  write(*,*) "flag_real8=",flag_real8, cbuf
  write(*,*) "iflag_abs =",iflag_abs
  write(*,*) "file_name =",file_name
  write(*,*) "MBZ, MSP  =",MBZ, MSP
  write(*,*) "ndata     =",ndata

  allocate( idata(3,ndata) ) ; idata=0

  do i=1,ndata
     read(u3,*) idata(1:3,i)
  end do

! ---

  read(u2,*) cbuf,ax
  read(u2,*) cbuf,aa(1:3,1)
  read(u2,*) cbuf,aa(1:3,2)
  read(u2,*) cbuf,aa(1:3,3)
  read(u2,*)
  read(u2,*) MKI,MI,zatom(1:MKI)

  allocate( asi(3,MI) ) ; asi=0.0d0
  allocate( rsi(3,MI) ) ; asi=0.0d0
  allocate( Kion(MI)  ) ; Kion=0

  do i=1,MI
     read(u2,*) Kion(i),asi(1:3,i)
  end do

  aa(:,:)=ax*aa(:,:)

  rsi(:,:)=matmul( aa,asi )

! ---

  open(u1,file=file_name,status='old',form='unformatted')

  read(u1) ML,ML1,ML2,ML3
  read(u1) MB,MB1,MB2

  allocate( LL(3,ML)        ) ; LL=0.0d0
  allocate( occ(MB,MBZ,MSP) ) ; occ=0.0d0

  read(u1) LL
  read(u1) occ

  write(*,*) "ML=",ML
  write(*,*) "ML1,ML2,ML3=",ML1,ML2,ML3
  write(*,*) "MB=",MB
  write(*,*) "MB1,MB2=",MB1,MB2
  write(*,*) "sum(occ)=",sum(occ)
  write(*,*) "min(occ),max(occ)=",minval(occ),maxval(occ)

  if ( flag_real8 ) then
     allocate( rtmp(ML) ) ; rtmp=0.0d0
  else
     allocate( ztmp(ML) ) ; ztmp=(0.0d0,0.0d0)
  end if

  allocate( rho(0:ML1-1,0:ML2-1,0:ML3-1) ) ; rho=0.0d0

  do s=1,MSP
  do k=1,MBZ
  do n=MB1,MB2

     if ( flag_real8 ) then
        read(u1) rtmp
     else
        read(u1) ztmp
     end if

     do i=1,ndata
        if ( idata(1,i) == n .and. &
             idata(2,i) == k .and. &
             idata(3,i) == s ) then
           call gen_cube(n,k,s)
        end if
     end do

  end do ! n
  end do ! k
  end do ! s

  close(u1)

CONTAINS


  SUBROUTINE gen_cube(n,k,s)

    implicit none

    integer,intent(IN) :: n,k,s
    integer :: i
    real(8) :: aa_del(3,3), r0(3)
    character(30) :: name, file_out
    character(1) :: cs
    character(3) :: ck
    character(5) :: cn

    aa_del(:,1)=aa(:,1)/ML1
    aa_del(:,2)=aa(:,2)/ML2
    aa_del(:,3)=aa(:,3)/ML3

    name  = 'SYS1'
    r0(:) = 0.0d0

    write(cn,'(i5.5)') n
    write(ck,'(i3.3)') k
    write(cs,'(i1.1)') s

    file_out = "wf_"//cn//"_"//ck//"_"//cs//".cub"

    select case( iflag_abs )
    case( 0 )
       if ( flag_real8 ) then
          do i=1,ML
             rho(LL(1,i),LL(2,i),LL(3,i)) = rtmp(i)
          end do
       else
          stop "complex wf can not be plot as it is"
       end if
    case( 1 )
       if ( flag_real8 ) then
          do i=1,ML
             rho(LL(1,i),LL(2,i),LL(3,i)) = abs( rtmp(i) )
          end do
       else
          do i=1,ML
             rho(LL(1,i),LL(2,i),LL(3,i)) = abs( ztmp(i) )
          end do
       end if
    case( 2 )
       if ( flag_real8 ) then
          do i=1,ML
             rho(LL(1,i),LL(2,i),LL(3,i)) = rtmp(i)**2
          end do
       else
          do i=1,ML
             rho(LL(1,i),LL(2,i),LL(3,i)) = abs( ztmp(i) )**2
          end do
       end if
    end select

    open(iu,file=file_out)
    write(iu,100) name
    write(iu,100) name
    write(iu,110) MI,r0
    write(iu,110) ML1,aa_del(:,1)
    write(iu,110) ML2,aa_del(:,2)
    write(iu,110) ML3,aa_del(:,3)
    do i=1,MI
       write(iu,110) zatom(Kion(i)),real(zatom(Kion(i)),8),rsi(:,i)
    end do
    do i1=0,ML1-1
       do i2=0,ML2-1
          write(iu,*) rho(i1,i2,:)
       end do
    end do
    close(iu)

100 format(a4)
110 format(i5,4f12.6)

  END SUBROUTINE gen_cube


END PROGRAM cubegen_wf
