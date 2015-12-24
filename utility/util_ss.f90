PROGRAM band_plot

  implicit none

  integer,parameter :: u1=5,u2=110,u3=120,max_loop=100000
  integer,parameter :: jyuu=10
  integer,allocatable :: mbi(:),maxloc10(:,:,:),indx(:)
  integer,allocatable :: itmp(:),itmp1(:),itmp2(:),jyuu_tmp(:,:)
  integer,allocatable :: itmp4(:),jtmp4(:),itmp3(:),jtmp3(:)
  integer,allocatable :: iiii(:),jjjj(:),kkkk(:),itp(:),jtp(:)
  integer :: n,mb_max,mb_ovlp,mb_min,mbv,lll,ll,loop,ll0,lll0
  integer :: i,j,k,l,nbk,mbc,mb1,mb2,mb,jj0,j1,i1,ndata
  integer :: i0,iwork(8),jwork(8),idummy,idummy2,mbk,iblk,nblk
  integer :: jj,jj1,ii1,i_fixed,l0,ierr,m1,m0,m0_max,m1_max
  real(8),parameter :: HT=27.2116d0
  real(8),allocatable :: eval(:,:),esp(:,:,:),dk_bz(:,:),dval(:,:)
  real(8),allocatable :: sqovlp10(:,:,:),rtmp3(:),rtmp4(:)
  real(8) :: rdummy(3),rr0,dummy,evb,ecb,r1,bb(3,3),r

!  open(u1,file="band_eigv_tot",status='old')
!  open(u2,file="band_ovlp_tot",status='old',form="unformatted")
!  open(u3,file="band_dedk_tot",status='old')

  read(u1,*) bb(1:3,1)
  read(u1,*) bb(1:3,2)
  read(u1,*) bb(1:3,3)

  mb_max=0
  do loop=1,max_loop
     read(u1,*) mbv,mb
     write(*,*) "mb,mbv =",mb,mbv
     mb_max=max(mb_max,mb)
     do i=1,mb
        read(u1,*)
     end do
     read(u1,*,END=10)
     backspace(u1)
  end do
  stop "stop!!! too much data"
10 continue
  nbk=loop

  write(*,*) "mb_max =",mb_max
  write(*,*) "nbk    =",nbk

  allocate( mbi(nbk)            ) ; mbi=0
  allocate( dk_bz(3,nbk)        ) ; dk_bz=0.d0
  allocate( eval(mb_max,nbk)    ) ; eval=0.d0
  allocate( dval(mb_max,nbk)    ) ; dval=0.d0
  allocate( esp(nbk,mb_max,nbk) ) ; esp=0.d0

  rewind u1
  read(u1,*)
  read(u1,*)
  read(u1,*)
  do k=1,nbk
     read(u1,*) idummy,mbi(k),dk_bz(1:3,k)
     do i=1,mbi(k)
        read(u1,*) idummy,eval(i,k),dummy
     end do
  end do

  goto 9
  m1_max=0
  m0_max=0
  rewind u2
  do k=2,nbk
     read(u2) m1,m0
     write(*,*) k,m1,m0
     if ( m1 > m1_max .or. m0 > m0_max ) then
        if ( allocated(maxloc10) ) deallocate(maxloc10)
        if ( allocated(sqovlp10) ) deallocate(sqovlp10)
        m1_max=max(m1,m1_max)
        m0_max=max(m0,m0_max)
        allocate( maxloc10(m1_max,m0_max,nbk) )
        allocate( sqovlp10(m1_max,m0_max,nbk) )
     end if
     read(u2) maxloc10(1:m1,1:m0,1)
     read(u2) sqovlp10(1:m1,1:m0,1)
  end do
  maxloc10(:,:,:)=0
  sqovlp10(:,:,:)=0.d0
  rewind u2
  do k=2,nbk
     read(u2) m1,m0
     write(*,*) k,m1,m0
     read(u2) maxloc10(1:m1,1:m0,k)
     read(u2) sqovlp10(1:m1,1:m0,k)
  end do

  rewind u3
  do k=1,nbk
     read(u3,*)
     do i=1,mbi(k)
        read(u3,*) idummy,dummy,rdummy(1:2),dval(i,k)
     end do
  end do
9 continue
!  close(u1)
!  close(u2)
!  close(u3)

!------------------------------

!  call plot_eval(mb1,mb2,ecb,HT)

!------------------------------

!  do k=1,1
!     do i=1,mbi(k)
!        write(*,'(1x,2i6,f10.5,5(1x,i4,f10.5))') k,i,eval(i,k)
!     end do
!  end do
!  do k=2,nbk
!     do i=1,mbi(k)
!        write(*,'(1x,2i6,f10.5,5(1x,i4,f10.5))') &
!             k,i,eval(i,k),(maxloc10(i,j,k),sqovlp10(i,j,k),j=1,3)
!     end do
!  end do
!  stop

  rewind 10
  n=minval(mbi)
  do i=1,maxval(mbi)
  do k=1,nbk
     if ( eval(i,k) == 0.d0 ) cycle
     write(10,'(1x,i8,50f10.5)') k,eval(i,k)
  end do
  write(10,*)
  write(10,*)
  end do
  stop

  esp(:,:,:)=1.d10

  iblk=0

  do mbk=nbk,1,-1

     iblk=iblk+1
     write(*,*) "mbk=",mbk

     do i=1,mbi(mbk)
        if ( eval(i,mbk) == 1.d10 ) cycle
        esp(mbk,i,iblk)=eval(i,mbk)
        eval(i,mbk)=1.d10
     end do

     do i=1,mbi(mbk)
        i0=i
        write(*,*) "i0=",mbk,i0
        do k=mbk,2,-1
           j=maxloc10(i0,1,k)
           r=sqovlp10(i0,1,k)
           if ( r < 0.3d0 ) exit
           if ( eval(j,k-1) == 1.d10 ) exit
           write(*,*) "i0=",k-1,j,r
           esp(k-1,i,iblk)=eval(j,k-1)
           eval(j,k-1)=1.d10
           i0=j
        end do
     end do

  end do ! mbk

  nblk=iblk

  j=0
  do iblk=1,nblk
     do i=1,maxval(mbi)
        if ( all(esp(:,i,iblk)==1.d10) ) cycle
        do k=1,nbk
           if ( esp(k,i,iblk) == 1.d10 ) cycle
           write(10,'(1x,i8,50f10.5)') k,esp(k,i,iblk)
        end do
        write(10,*)
        write(10,*)
        j=j+1
        write(10,'("#",i6)') j
     end do
  end do


  write(*,*) count(eval==1.d10),count(eval==0.d0),count(eval/=0.d0),size(eval)

  rewind 11
  write(11,'(1x,"plot \")')
  write(11,'(1x,a16,i4," w lp \")') ' "fort.10" index',0
  do n=1,j
     write(11,'(1x,a16,i4," w lp \")') ',"fort.10" index',n
  end do
  stop

END PROGRAM band_plot
