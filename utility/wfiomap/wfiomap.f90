PROGRAM wfiomap

  implicit none

  integer :: u1=5, u2=10
  integer :: ML1,ML2,ML3
  integer :: MBN,MBZ,MSP,i1,i2,i3,i4,i5,i6,i,irank,jrank,n
  integer :: np_old(6),np_new(6),nprocs_old,nprocs_new
  integer,allocatable :: ntmp(:,:),Igrid(:,:,:),Jgrid(:,:,:)
  integer,allocatable :: itmp(:),iomap(:,:,:),nmap(:)
  integer :: a1b,b1b,a2b,b2b,a3b,b3b,a1a,b1a,a2a,b2a,a3a,b3a

  read(u1,*) ML1,ML2,ML3
  read(u1,*) MBN,MBZ,MSP
  write(*,*) "ML1,ML2,ML3",ML1,ML2,ML3
  write(*,*) "MBN,MBZ,MSP",MBN,MBZ,MSP

  np_old=1
  np_new=1
  read(u1,*) np_old(1:6)
  read(u1,*) np_new(1:6)
  write(*,*) "np_old",np_old
  write(*,*) "np_new",np_new

  nprocs_old=1
  nprocs_new=1
  do i=1,6
     nprocs_old=nprocs_old*np_old(i)
     nprocs_new=nprocs_new*np_new(i)
  end do
  write(*,*) "nprocs_old",nprocs_old
  write(*,*) "nprocs_new",nprocs_new

  i1=maxval( np_old(1:3) )
  i2=maxval( np_new(1:3) )
  i3=max( i1,i2 )
  allocate( ntmp(0:i3-1,3) )
!
  ntmp=0
  do i3=0,ML3-1
     i=mod(i3,np_old(3))
     ntmp(i,3)=ntmp(i,3)+1
  end do
  do i2=0,ML2-1
     i=mod(i2,np_old(2))
     ntmp(i,2)=ntmp(i,2)+1
  end do
  do i1=0,ML1-1
     i=mod(i1,np_old(1))
     ntmp(i,1)=ntmp(i,1)+1
  end do
!  write(*,'(1x,10i4)') ntmp(0:np_old(1)-1,1)
!  write(*,'(1x,10i4)') ntmp(0:np_old(2)-1,2)
!  write(*,'(1x,10i4)') ntmp(0:np_old(3)-1,3)
!
  allocate( Igrid(2,3,0:nprocs_old-1) ) ; Igrid=0

  irank=-1
  do i6=1,np_old(6)
  do i5=1,np_old(5)
  do i4=1,np_old(4)
     do i3=0,np_old(3)-1
     do i2=0,np_old(2)-1
     do i1=0,np_old(1)-1
        irank=irank+1
        Igrid(1,1,irank) = sum( ntmp(0:i1,1) ) - ntmp(i1,1)
        Igrid(2,1,irank) = Igrid(1,1,irank) + ntmp(i1,1) - 1
        Igrid(1,2,irank) = sum( ntmp(0:i2,2) ) - ntmp(i2,2)
        Igrid(2,2,irank) = Igrid(1,2,irank) + ntmp(i2,2) - 1
        Igrid(1,3,irank) = sum( ntmp(0:i3,3) ) - ntmp(i3,3)
        Igrid(2,3,irank) = Igrid(1,3,irank) + ntmp(i3,3) - 1
     end do
     end do
     end do
  end do
  end do
  end do

!  do irank=0,nprocs_old-1
!     write(*,'(1x,i4,2x,3(2i6,1x))') irank,Igrid(:,:,irank)
!  end do
!
  ntmp=0
  do i3=0,ML3-1
     i=mod(i3,np_new(3))
     ntmp(i,3)=ntmp(i,3)+1
  end do
  do i2=0,ML2-1
     i=mod(i2,np_new(2))
     ntmp(i,2)=ntmp(i,2)+1
  end do
  do i1=0,ML1-1
     i=mod(i1,np_new(1))
     ntmp(i,1)=ntmp(i,1)+1
  end do

!  write(*,'(1x,10i4)') ntmp(0:np_new(1)-1,1)
!  write(*,'(1x,10i4)') ntmp(0:np_new(2)-1,2)
!  write(*,'(1x,10i4)') ntmp(0:np_new(3)-1,3)

!---

  allocate( Jgrid(2,3,0:nprocs_new-1) ) ; Jgrid=0

  irank=-1
  do i6=1,np_new(6)
  do i5=1,np_new(5)
  do i4=1,np_new(4)
     do i3=0,np_new(3)-1
     do i2=0,np_new(2)-1
     do i1=0,np_new(1)-1
        irank=irank+1
        Jgrid(1,1,irank) = sum( ntmp(0:i1,1) ) - ntmp(i1,1)
        Jgrid(2,1,irank) = Jgrid(1,1,irank) + ntmp(i1,1) - 1
        Jgrid(1,2,irank) = sum( ntmp(0:i2,2) ) - ntmp(i2,2)
        Jgrid(2,2,irank) = Jgrid(1,2,irank) + ntmp(i2,2) - 1
        Jgrid(1,3,irank) = sum( ntmp(0:i3,3) ) - ntmp(i3,3)
        Jgrid(2,3,irank) = Jgrid(1,3,irank) + ntmp(i3,3) - 1
     end do
     end do
     end do
  end do
  end do
  end do

!  do irank=0,nprocs_new-1
!     write(*,'(1x,i4,2x,3(2i6,1x))') irank,Jgrid(:,:,irank)
!  end do
!---

  allocate( itmp(nprocs_old) ) ; itmp=0
  allocate( iomap(0:6,nprocs_old,nprocs_new) ) ; iomap=0
  allocate( nmap(nprocs_new) ) ; nmap=0

  open(u2,file="input_wfiodir")
  write(u2,*) nprocs_new

  do jrank=0,nprocs_new-1

     a1b = Jgrid(1,1,jrank)
     b1b = Jgrid(2,1,jrank)
     a2b = Jgrid(1,2,jrank)
     b2b = Jgrid(2,2,jrank)
     a3b = Jgrid(1,3,jrank)
     b3b = Jgrid(2,3,jrank)

     write(*,'(1x,"new rank: ",i4,1x,6i4)') jrank,a1b,b1b,a2b,b2b,a3b,b3b

     n=0
     do irank=0,nprocs_old-1

        a1a = Igrid(1,1,irank)
        b1a = Igrid(2,1,irank)
        a2a = Igrid(1,2,irank)
        b2a = Igrid(2,2,irank)
        a3a = Igrid(1,3,irank)
        b3a = Igrid(2,3,irank)

        if ( (a1a <= a1b .and. a1b <= b1a .or. a1a <= b1b .and. b1b <= b1a) &
        .and.(a2a <= a2b .and. a2b <= b2a .or. a2a <= b2b .and. b2b <= b2a) &
        .and.(a3a <= a3b .and. a3b <= b3a .or. a3a <= b3b .and. b3b <= b3a) &
           ) then
           write(*,'(1x,"    rank: ",i4,1x,6i4)') irank,a1a,b1a,a2a,b2a,a3a,b3a
           n=n+1
           itmp(n)=irank
           iomap(0,n,jrank+1) = irank
           iomap(1,n,jrank+1) = a1a
           iomap(2,n,jrank+1) = b1a
           iomap(3,n,jrank+1) = a2a
           iomap(4,n,jrank+1) = b2a
           iomap(5,n,jrank+1) = a3a
           iomap(6,n,jrank+1) = b3a
        end if

     end do ! irnak(old)

     nmap(jrank+1) = n

     write(u2,*) n
     do i=1,n
        write(u2,*) itmp(i)
     end do

  end do ! jrank(new)

  do jrank=0,nprocs_new-1
     write(u2,'(1x,"new_rank: ",i4,1x,6i4)') jrank,Jgrid(:,:,jrank)
     do i=1,nmap(jrank+1)
        irank=iomap(0,i,jrank+1)
        write(u2,'(1x,"    rank: ",i4,1x,6i4)') irank,Igrid(:,:,irank)
     end do
  end do
        
  close(u2)

  deallocate( nmap )
  deallocate( iomap )
  deallocate( itmp )

!---

  deallocate( Jgrid )
  deallocate( Igrid )
  deallocate( ntmp )

END PROGRAM wfiomap
