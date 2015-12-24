
      PROGRAM MAIN

      implicit none
      integer,parameter :: MI0=8, NI=200000
      integer,parameter :: MKI=2, nmkd=1
      integer :: nc1,nc2,nc3,nprocs,u,nc
      integer :: a,a0,i,j,k,m,n,m0,n0,MI,MB,ML1,ML2,ML3,ML0,ML
      integer :: NPROW,NPCOL,MBSIZE,NBSIZE,np1,np2,np3,NPCOL0
      integer :: Kion(NI),mkd(NI),Kion0(MI0),mkd0(MI0),np_2d(1:4)
      real(8) :: asi0(3,MI0),asi(3,NI),RR(3),aa(3,3)
      real(8) :: ax,ax0,mem

!
! input
!
      Write(*,*) "input size of the unit cell : nc1,nc2,nc3="
      read(*,*) nc1,nc2,nc3

      write(*,*) "input # of nodes: np_2d(1:4) ="
      read(*,*) np_2d(1:4)

      nprocs=np_2d(1)*np_2d(2)*np_2d(3)*np_2d(4)

!===========================================================================
! # of meshes along each direction
!( If you need more dense grid points, you should make larger the paramter.)
!
      ML0=20
!
! Lattice constant
!
!      ax0 = 8.2392d0
      ax0 = 8.1392d0
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
      Kion0(:)=1 ; Kion0(2)=2 ; Kion0(6:8)=2
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
      u=1
      rewind u
      rewind 970
      write(u,'(a7)') "# RSDFT"
      write(u,'(a9," /")') "'LDAPZ81'"
      write(u,*) ax
      write(u,'(3f10.5)') aa(:,1)
      write(u,'(3f10.5)') aa(:,2)
      write(u,'(3f10.5)') aa(:,3)
      write(970,'(i4,i6,i4)') MKI,MI
      do a=1,MI
         write(970,'(i4,3f15.10,i4)') Kion(a),asi(:,a)
      end do
      write(u,'(2i8," /Nband Nspin")') MB, 1
      write(u,'(f8.3," /Next_electron")') 0.0
      write(u,'(i2," / nk")') 2
      write(u,'(3i2)') 2,2,2
      write(u,'(3i2)') 2,2,2
      write(u,'(2i3," /Ncg, iswitch_gs")') 2,0
      write(u,'(3i5)') ML1,ML2,ML3
      write(u,'(i2," /Md")') 6
      write(u,'(i2,a13," /")') 2," 'Si_psv.dat'"
      write(u,'(i2,a13," /")') 2," 'C_psv.dat'"
      write(u,'(i2,i3,f5.1,g8.1)') 10,4,1.0,1.d-20
      write(u,'(i2)') 2
      write(u,*) 1,0
      write(u,'(3f6.3)') 1.5d0,1.d0,8.d0
      write(u,'(6i4,12x,"/ np_2d(1:6)")') np_2d(1:4),1,1
      write(u,'(i2)') 0
      write(u,'(2i4," /NBLK,NBLK1")') 0,0
      write(u,'(i2," /mloop")') 3
      write(u,'(g10.3," /")') 1.d-6
      write(u,'(i4," / Diter")') 100
      write(u,'(3i5)') 0,100,3
      write(u,'(3i3)') 1,0,0
      write(u,*) 1.d10
      write(u,'(3i5," /")') 0,6,5
      write(u,'(4g10.3)') 0.5d0,1.d-10,1.d-4,1.d-1
      write(u,'(2i5,12x,"/")') 0,0

      stop
      END PROGRAM MAIN
