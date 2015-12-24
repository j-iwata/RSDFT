PROGRAM ufld2gp

  implicit none
  integer,parameter :: u1=1, u2=10, u3=11
  integer,parameter :: maxloop=1000000
  integer :: loop,ik,k,n,i,j,nk,nk0,nb,ne,kmax
  character(5) :: cbuf
  real(8),parameter :: aB=0.529177d0, HT=27.2116d0
  real(8),allocatable :: kxyz_pc(:,:),kxyz_sc(:,:)
  real(8),allocatable :: esp(:,:), weight(:,:)
  real(8) :: emax,emin,e_mergin,eta,de,pi,f,e,x

! ---

  ne = 500

  eta = 0.10d0/27.2116d0

  pi = acos(-1.0d0)

! ---

  write(*,*) "eta(eV), ne="
  read(*,*) eta, ne
  eta = eta/HT

! ---

  open(u1,file="band_ufld",status="old")

  nb=0
  nk=0
  kmax=0
  do loop=1,maxloop
     read(u1,*,END=1) cbuf,k
     kmax=max(k,kmax)
     if ( cbuf == "iktrj" ) then
        nk=nk+1
        nb=0
     else
        nb=nb+1
     end if
  end do
1 continue

  write(*,*) "nk, nb=",nk,nb
  write(*,*) "kmax=",kmax

  nk0=nk
  nk=kmax

! ---

  allocate( kxyz_pc(3,0:nk) ) ; kxyz_pc=0.0d0
  allocate( kxyz_sc(3,nk) ) ; kxyz_sc=0.0d0
  allocate( esp(nb,nk)    ) ; esp=0.0d0
  allocate( weight(nb,nk) ) ; weight=0.0d0

! ---

  rewind u1

  do ik=1,nk0

     read(u1,*) cbuf,k,kxyz_pc(1:3,k),kxyz_sc(1:3,k)

     do n=1,nb

        read(u1,*) i,j,esp(n,k),weight(n,k)

     end do

  end do

! ---

  emin = minval( esp )
  emax = maxval( esp )

  write(*,*) "emin, emax=",emin,emax

  e_mergin = (emax-emin)*0.1d0

  emax = emax + e_mergin
  emin = emin - e_mergin

  write(*,*) "emin, emax, e_mergin=",emin,emax,e_mergin

  de = (emax-emin)/ne

  write(*,*) "ne, de=",ne,de

  write(*,*) "eta=",eta

! ---

  rewind u2
  rewind u3

  kxyz_pc(:,0) = kxyz_pc(:,1)
  x=0.0d0

  do k=1,nk

     x = x + sqrt( sum( (kxyz_pc(:,k)-kxyz_pc(:,k-1))**2 ) )/aB

     do i=1,ne

        e = emin + (i-1)*de

        f=0.0d0
        do n=1,nb
           f = f + eta/( (e-esp(n,k))**2 + eta**2 )*weight(n,k)
        end do
        f=f/pi

        write(u2,'(1x,3f20.15)') x,e*HT,f
        write(u3,'(1x,3f20.15)') x,e*HT,log(f)

     end do ! i

     write(u2,*)
     write(u3,*)

  end do ! k

! ---

  deallocate( weight  )
  deallocate( esp     )
  deallocate( kxyz_sc )
  deallocate( kxyz_pc )

  close(u1)

END PROGRAM ufld2gp
