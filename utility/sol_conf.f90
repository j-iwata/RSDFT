
      PROGRAM MAIN

      implicit none
      integer,parameter :: MI0=8, NI=200000
      integer,parameter :: MKI=1, nmkd=1
      integer :: nc1,nc2,nc3,nprocs,u1,u2,nc
      integer :: a,a0,i,j,k,m,n,m0,n0,MI,MB,ML1,ML2,ML3,ML0,ML
      integer :: NPROW,NPCOL,MBSIZE,NBSIZE,np1,np2,np3,NPCOL0
      integer :: Kion(NI),mkd(NI),Kion0(MI0),mkd0(MI0),np_2d(1:4)
      real(8) :: asi0(3,MI0),asi(3,NI),RR(3),aa(3,3)
      real(8) :: ax,ax0,mem

!
! input
!
      write(*,*) "input size of the unit cell : nc1,nc2,nc3="
      write(*,*) "( # of atoms is 8*nc1*nc2*nc3 )"
      read(*,*) nc1,nc2,nc3

!      write(*,*) "input # of nodes: np_2d(1:4) ="
!      read(*,*) np_2d(1:4)
      np_2d(1:4)=1

      nprocs=np_2d(1)*np_2d(2)*np_2d(3)*np_2d(4)

!===========================================================================
! # of meshes along each direction
!( If you need more dense grid points, you should make larger the paramter.)
!
      ML0=16
!      ML0=32
!
! Lattice constant
!
      ax0 = 10.261d0
!      ax0 = 10.16d0
!
!( Grid spacings are defined as H1=H2=H3=ax0/ML0. )
!===========================================================================

      ax = ax0*max(nc1,nc2,nc3)

!
! ---
!
      np1=np_2d(1) ; np2=np_2d(2) ; np3=np_2d(3)
!
! # of atoms
!
      MI = MI0 * nc1*nc2*nc3

!
! # of grid points
!
      ML1=ML0*nc1 ; ML2=ML0*nc2 ; ML3=ML0*nc3

      ML=ML1*ML2*ML3

!
! # of bands
!
      MB=2*MI+nint(0.2d0*MI)

!
! parameters for SCALAPACK
!
      n      = np_2d(1)*np_2d(2)*np_2d(3)
      NPCOL0 = np_2d(1)*np_2d(2)*np_2d(3)*np_2d(4)
      NPCOL  = NPCOL0
      NPROW  = 1
      do i=2,n
         j=NPCOL0/i
         if ( i*j==NPCOL0 .and. i<=j .and. j-i<NPCOL-NPROW ) then
            NPCOL=j
            NPROW=i
         end if
      end do

      i=(MB+NPROW-1)/NPROW
      j=(MB+NPCOL-1)/NPCOL
      MBSIZE=min(i,j)
      NBSIZE=MBSIZE

      if ( NBSIZE*NPCOL/=MB ) then
         n=max( NBSIZE*NPCOL, MB )
         n=min( n, (NBSIZE+1)*NPCOL )
         MB=n
      end if

      NPROW=0
      NPCOL=0
      MBSIZE=0
      NBSIZE=0

!
! write parameters
!
      write(*,*) "# of nodes =",nprocs
      write(*,*) "# of atoms =",MI
      write(*,*) "# of grid points =",ML
      write(*,*) "# of eigenstates",MB

!
! memory estimation
!
      mem=16.d0*MB*dble(ML)/dble(nprocs)
      write(*,*) "memory estimation (GB) =",mem/1024.d0**3

!
! Lattice vectors
!
      aa=0.d0
      aa(1,1)=ax0*nc1/ax
      aa(2,2)=ax0*nc2/ax
      aa(3,3)=ax0*nc3/ax

!
! atomic coorinates
!
      asi0(:,1)=(/0.00d0, 0.00d0, 0.00d0/)
      asi0(:,2)=(/0.25d0, 0.25d0, 0.25d0/)
      asi0(:,3)=(/0.50d0, 0.00d0, 0.50d0/)
      asi0(:,4)=(/0.00d0, 0.50d0, 0.50d0/)
      asi0(:,5)=(/0.50d0, 0.50d0, 0.00d0/)
      asi0(:,6)=(/0.75d0, 0.25d0, 0.75d0/)
      asi0(:,7)=(/0.25d0, 0.75d0, 0.75d0/)
      asi0(:,8)=(/0.75d0, 0.75d0, 0.25d0/)
      Kion0(:)=1
      mkd0(:)=1

      a=0
      do k=1,nc3
      do j=1,nc2
      do i=1,nc1

         RR(1)=i-1 ; RR(2)=j-1 ; RR(3)=k-1

         do a0=1,MI0
            a=a+1
            asi(:,a)=RR+asi0(:,a0)
            Kion(a)=Kion0(a0)
            mkd(a)=mkd0(a0)
         end do

      end do
      end do
      end do

      if ( a/=MI ) then
         write(*,*) "a,MI=",a,MI
         stop
      end if

      do a=1,MI
         asi(1,a)=asi(1,a)/nc1
         asi(2,a)=asi(2,a)/nc2
         asi(3,a)=asi(3,a)/nc3
      end do

!
! input file
!
      u1=1
      u2=970
      rewind u1
      write(u1,'(a7)') "# RSDFT"
      write(u1,'("NBAND",i8," /")') MB
      write(u1,'("NCG  ",i8," /")') 3
      write(u1,'("ICG  ",i8," /")') 1
      write(u1,'("IPC  ",i8," /")') 1
      write(u1,'("SWGS ",i8," /")') 0
      write(u1,'("NGRID",3i5)') ML1,ML2,ML3
      write(u1,'("Md   ",i3)') 6
      write(u1,'("PP 2 ",a12)') "'Si_psv.dat'"
      write(u1,'("imix ",i4)') 10
      write(u1,'("beta ",f10.5)') 1.0
      write(u1,'("scfconv ",g20.6)') 1.d-15
      write(u1,'("PROCS",4i5,"  /")') 1,1,1,1
      write(u1,'("MBd  ",i3)') 0
      write(u1,'("NBLK ",i3)') 0
      write(u1,'("NBLK1",i3)') 0
      write(u1,'("Diter",i5)') 100
      write(u1,'("Nsweep",i4)') 0
      write(u1,'("IC   ",i5)') 0
      write(u1,'("OC   ",i5)') 3
      write(u1,'("OC2  ",i5)') 100
      write(u1,'("IOCTRL",i5)') 0
      write(u1,'("SWSCF ",i5)') 1
      write(u1,'("SWOPT ",i5)') 0
      write(u1,'("ETLIMIT",g20.8)') 1.d10
      write(u1,'("SCL   ",2i5)') 0,0
      rewind u2
      write(u2,'("AX  ",f20.15)') ax
      write(u2,'("A1  ",3f10.5)') aa(:,1)
      write(u2,'("A2  ",3f10.5)') aa(:,2)
      write(u2,'("A3  ",3f10.5)') aa(:,3)
      write(u2,'("AA")')
      write(u2,'(i4,i6,"  /")') MKI,MI
      do a=1,MI
         write(u2,'(i4,3f15.10,i4,"  /")') Kion(a),asi(:,a),1
      end do

      stop
      END PROGRAM MAIN
