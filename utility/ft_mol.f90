PROGRAM ft_mol

  implicit none

  integer,parameter :: u0=1, u5=5, u1=11, u2=12, u3=13
  integer,parameter :: max_loop=10000000
  integer :: loop,nt,ne,it,ie
  real(8),parameter :: Tau2fs = 2.418884326505d-2
  real(8),parameter :: HT=27.2116d0
  real(8),allocatable :: dipole(:,:)
  real(8) :: t1,t2,dt,t,Tmax,pi
  real(8) :: emin,emax,de,gamma,e,de_max
  complex(8),parameter :: z0=(0.0d0,0.0d0),zi=(0.0d0,1.0d0)
  complex(8) :: zsum1,zsum2,zsum3,phase,ze

  open(u0,file="tddft_data",status="old")

  read(u0,*) t1
  read(u0,*) t2

  dt = t2-t1
  write(*,*) "dt(au,fs)=",dt,dt*Tau2fs

  rewind u0
  do loop=1,max_loop
     read(u0,*,END=10)
  end do
10 nt=loop-2

  write(*,*) "nt       =",nt

  Tmax = nt*dt
  write(*,*) "Tmax(au,fs)=",Tmax,Tmax*Tau2fs

  pi = acos(-1.0d0)
  de_max = pi/Tmax
  write(*,*) "Energy resolution (au,eV)=",de_max,de_max*HT

! ---

  allocate( dipole(3,0:nt) ) ; dipole=0.0d0

! ---

  rewind u0
  do it=0,nt
     read(u0,*) t1,dipole(1:3,it)
  end do

  do it=1,nt
     dipole(1,it) = dipole(1,it) - dipole(1,0)
     dipole(2,it) = dipole(2,it) - dipole(2,0)
     dipole(3,it) = dipole(3,it) - dipole(3,0)
  end do

! ---

  close(u0)

! ---

  write(*,'(a20)') repeat("-",20)

!  ne = 10000
!  emax = 10.0d0/HT
!  emin = 0.0d0
!  gamma = 0.1d0/HT

  write(*,*) "emin,emax (in eV) = ?,?"
  read(u5,*) emin,emax
  emin=emin/HT
  emax=emax/HT
  write(*,*) "# of energy grid = ?"
  read(u5,*) ne
!  write(*,*) "gamma (in eV) = ?"
!  read(u5,*) gamma
!  gamma=gamma/HT
  gamma=de_max
  write(*,*) "gamma(au,eV)=",gamma,gamma*HT
  de = ( emax - emin )/ne
  write(*,*) "de(au,eV)=",de,de*HT

  write(*,'(a20)') repeat("-",20)

! ---

  rewind u1
  rewind u2
  rewind u3

  do ie=0,ne

     e = emin + ie*de
     ze = e + zi*gamma

     zsum1=z0
     zsum2=z0
     zsum3=z0
     do it=1,nt

        t = it*dt

        phase = exp(zi*ze*t)
        
        zsum1 = zsum1 + phase*dipole(1,it)
        zsum2 = zsum2 + phase*dipole(2,it)
        zsum3 = zsum3 + phase*dipole(3,it)

     end do ! it

     zsum1=zsum1*dt
     zsum2=zsum2*dt
     zsum3=zsum3*dt

     e=e*HT
     write(u1,*) e,real(zsum1),aimag(zsum1)
     write(u2,*) e,real(zsum2),aimag(zsum2)
     write(u3,*) e,real(zsum3),aimag(zsum3)

  end do ! ie

! ---

  deallocate( dipole )

! ---

END PROGRAM ft_mol
