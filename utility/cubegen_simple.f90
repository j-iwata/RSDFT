PROGRAM cubegen_simple

  implicit none

  integer,parameter :: u1=5, u2=970, iu=10
  real(8),parameter :: ab=0.529177d0
  integer :: idum,MB,MBZ,i,j,k,n,ML,ML1,ML2,ML3,MI,MKI
  integer :: NBLK1,NBLK2,NBLK,np1,np2,np3,i1,i2,i3,n1,n2,n3
  integer :: MBwr1,MBwr2,MBZ0,MB0
  integer :: zatom(5)
  real(8) :: aa(3,3),ax,r0(3),aL(3),aa_del(3,3),va,fac
  real(8) :: c,d,a,b
  integer,allocatable :: Kion(:),mkd(:)
  real(8),allocatable :: asi(:,:),rsi(:,:),rho(:,:,:),occ(:,:)
  real(8),allocatable :: rtmp(:)
  complex(8) ,parameter :: zi=(0.d0,1.d0)
  complex(8),allocatable :: unk(:,:,:,:,:),utmp(:)
  character(4) :: name,cbuf,ckey
  integer,allocatable :: ip1(:),ip2(:),ip3(:)
  integer,allocatable :: nip1(:),nip2(:),nip3(:)
  integer,allocatable :: LL2(:,:),LLL2(:,:,:)

! --- read input data ---

  read(u1,*) cbuf,ax
  do i=1,3
     read(u1,*) cbuf,aa(1:3,i)
  end do
  read(u1,*) ckey
  read(u1,*) MKI,MI,zatom(1:MKI)

  if ( MKI > 5 ) then
     write(*,*) "MKI=",MKI,size(zatom)
     stop
  end if

  allocate( asi(3,MI),Kion(MI),rsi(3,MI) )

  do i=1,MI
     read(u1,*) Kion(i),asi(:,i)
  end do

  aa=ax*aa
  do i=1,3
     aL(i)=sqrt( sum(aa(:,i)**2) )
  end do

  if ( ckey == "XYZ" ) then
     if ( ax == 1.0d0 ) then
        call convert_xyz2aa( MI, aa, asi )
        rsi=matmul( aa, asi )
     else
        rsi=asi
     end if
  else
     rsi=matmul( aa,asi )
  end if

  ML1=1.d0
  ML2=1.d0
  ML3=1.d0

  aa_del(:,1)=aa(:,1)/ML1
  aa_del(:,2)=aa(:,2)/ML2
  aa_del(:,3)=aa(:,3)/ML3

! --- cube file

  name='SYS1'
  r0(:)=0.d0

  open(iu,file='tmp.cube')
  write(iu,100) name
  write(iu,100) name
  write(iu,110) MI,r0
  write(iu,110) ML1,aa_del(:,1)
  write(iu,110) ML2,aa_del(:,2)
  write(iu,110) ML3,aa_del(:,3)
  do i=1,MI
     write(iu,110) zatom(Kion(i)),real(zatom(Kion(i)),8),rsi(:,i)
  end do
  close(iu)

100 format(a4)
110 format(i5,4f12.6)

  deallocate( rsi,Kion,asi )


CONTAINS


  SUBROUTINE convert_xyz2aa( MI, aa, asi )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(INOUT) :: aa(3,3),asi(3,MI)
    real(8) :: LL(1:2,3),RG(3)

    aa(:,:) = 0.0d0

    do i=1,3
       RG(i) = sum( asi(i,:) )
    end do

    write(*,'(1x,"RG=",3f10.5)') RG

    do i=1,3
       LL(1,i) = minval( asi(i,:) )
       LL(2,i) = maxval( asi(i,:) )
       aa(i,i) = LL(2,i) - LL(1,i) + 5.0d0
    end do

    do i=1,MI
       asi(1,i) = ( asi(1,i) + 0.5d0*aa(1,1) )/aa(1,1)
       asi(2,i) = ( asi(2,i) + 0.5d0*aa(2,2) )/aa(2,2)
       asi(3,i) = ( asi(3,i) + 0.5d0*aa(3,3) )/aa(3,3)
    end do

  END SUBROUTINE convert_xyz2aa


END PROGRAM cubegen_simple
