MODULE fermi_module

  use bberf_module

  implicit none

  PRIVATE
  PUBLIC :: read_fermi,ekbt,calc_fermi &
           ,efermi, Eentropy, read_oldformat_fermi

  real(8) :: ekbt, efermi, Eentropy
  integer :: mb1,mb2,kinteg

  integer :: nsetocc,isetocc(2)
  real(8) :: setocc(20)

  logical :: first_time = .true.
  real(8),allocatable :: factor(:)

CONTAINS


  SUBROUTINE read_fermi(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,n1,n2
    character(6) :: cbuf,ckey
    ekbt=1.d-5
    kinteg=5
    nsetocc=0
    isetocc=0
    setocc=0.0d0
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:4) == "EKBT" ) then
             backspace(unit)
             read(unit,*) cbuf,ekbt
          else if ( ckey(1:6) == "KINTEG" ) then
             backspace(unit)
             read(unit,*) cbuf,kinteg
          else if ( ckey(1:6) == "SETOCC" ) then
             backspace(unit)
             read(unit,*) cbuf,n1,n2,setocc(nsetocc+1:nsetocc+n2-n1+1)
             nsetocc=nsetocc+n2-n1+1
             if ( all( isetocc(:) == 0 ) ) then
                isetocc(1)=n1
                isetocc(2)=n2
             else
                if ( isetocc(1) /= n1 .or. isetocc(2) /= n2 ) then
                   write(*,*) "isetocc(:) are different !!!"
                   write(*,*) "isetocc(:)=",isetocc(:)
                   write(*,*) " /= n1,n2 =",n1,n2
                   stop "stop@read_fermi"
                end if
             end if
          end if
       end do
999    continue
       write(*,*) "ekbt   =",ekbt
       write(*,*) "kinteg =",kinteg
       write(*,*) "nsetocc=",nsetocc,setocc(1:nsetocc)
       write(*,*) "isetocc=",isetocc(1:2)
    end if
    call send_fermi(0)
  END SUBROUTINE read_fermi


  SUBROUTINE read_oldformat_fermi(rank,unit)
    implicit none
    integer,intent(IN) :: rank,unit
    if ( rank == 0 ) then
       read(unit,*) ekbt,kinteg
       write(*,*) "ekbt=",ekbt
       write(*,*) "kinteg=",kinteg
    end if
    call send_fermi(0)
  END SUBROUTINE read_oldformat_fermi


  SUBROUTINE send_fermi(rank)
    implicit none
    integer,intent(IN) :: rank
    integer :: ierr
    include 'mpif.h'
    call mpi_bcast(ekbt,1,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kinteg,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nsetocc,1,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(isetocc,2,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
    call mpi_bcast(setocc,nsetocc,MPI_REAL8,rank,MPI_COMM_WORLD,ierr)
  END SUBROUTINE send_fermi


  SUBROUTINE calc_fermi(iter,nfixed,MB,MBZ,MSP,znel,dspn,esp,wbz,occ &
                       ,disp_switch)
    implicit none
    integer,intent(IN)  :: iter,nfixed,MB,MBZ,MSP
    logical,intent(IN)  :: disp_switch
    real(8),intent(IN)  :: esp(MB,MBZ,MSP),wbz(MBZ),znel,dspn
    real(8),intent(OUT) :: occ(MB,MBZ,MSP)
    real(8) :: ef1,ef2,ef,ef0
    integer :: id,n,k,s,efconv,nn
    real(8) :: zne,octmp,ff,xx
    real(8),parameter :: eps=0.d0
    integer,parameter :: mxcycl=1000

    call write_border( 1, " calc_fermi_mp(start)" )

    if ( nsetocc > 0 ) then
       mb1=1
       mb2=MB
       occ(:,:,:)=0.0d0
       nn=nsetocc/MSP
       do s=1,MSP
       do k=1,MBZ
          do n=1,isetocc(1)-1
             occ(n,k,s)=1.0d0
          end do
          occ(isetocc(1):isetocc(2),k,s) = setocc(1+nn*(s-1):nn+nn*(s-1))
       end do
       end do
       ef=maxval( esp(isetocc(2),1:MBZ,1:MSP) )
       goto 100
    end if

    if ( first_time ) then
       first_time = .false.
       if ( kinteg > 0 ) then
          allocate( factor(kinteg) )
          factor(1)=-1.d0/(4.d0*sqrt(acos(-1.d0)))
          do n=2,kinteg
             factor(n)=-factor(n-1)/(4.d0*n)
          end do
       end if
    end if

    mb1 = 1
    mb2 = MB

! Set upper & lower boundarires of Fermi energy

    ef1 = minval( esp(mb1:mb2,1:MBZ,1:MSP) )
    ef2 = maxval( esp(mb1:mb2,1:MBZ,1:MSP) )
    if ( ef1 == ef2 ) then
       ef1 = ef1 - 0.01d0
       ef2 = ef2 + 0.01d0
    end if

!C Safety margin for highly degenerate systems & artificial fault
!C
    ef2 = ef2 + min( ekbt*1.d2, 0.1d0 )

    if ( MSP == 1 .or. iter > Nfixed ) then

       zne = znel - 2.d0*(mb1-1)

       ef0 = 1.d10

       do id=1,mxcycl

          ef = 0.5d0*( ef1 + ef2 )
          if ( ef == ef0 ) goto 100
          octmp = 0.0d0

          do s=1,MSP
          do k=1,MBZ
          do n=mb1,mb2

             xx=(esp(n,k,s)-ef)/ekbt
             ff=ff0(kinteg,xx)

             octmp = octmp + ff*wbz(k)*2.d0/dble(MSP)
             occ(n,k,s)=ff

          end do
          end do
          end do

          if ( octmp-zne > eps ) then
             ef2=ef
          else if ( octmp-zne < -eps ) then
             ef1=ef
          else
             goto 100
          end if

          ef0 = ef

       end do ! id

    else

       efconv=0
       if ( DISP_SWITCH ) then
          write(*,*) "total spin density is fixed!!"
       end if

       do s=1,MSP

          ef1 = minval( esp(mb1:mb2,1:MBZ,s) )
          ef2 = maxval( esp(mb1:mb2,1:MBZ,s) )
          if ( ef1 == ef2 ) then
             ef1 = ef1 - 0.01d0
             ef2 = ef2 + 0.01d0
          end if
          ef2 = ef2 + min( ekbt*1.d2, 0.1d0 )

          zne = 0.5d0*znel + (3-2*s)*0.5d0*dspn

          do id=1,mxcycl

             ef = 0.5d0*(ef1+ef2)
             octmp = 0.0d0

             do n=mb1,mb2

                do k=1,MBZ

                   xx=(esp(n,k,s)-ef)/ekbt
                   ff=ff0(kinteg,xx)

                   octmp=octmp+ff*wbz(k)
                   occ(n,k,s)=ff

                end do

             end do

             if ( octmp-zne > eps ) then
                ef2=ef
             else if ( octmp-zne < -eps ) then
                ef1=ef
             else
                efconv=efconv+1
                exit
             end if

          end do ! id

       end do ! s

       if ( efconv == 2 ) goto 100

    end if

    if ( abs(octmp-zne) > 1.d-10 ) then
       if ( disp_switch ) then
          write(6,*)' EF IS NOT CONVERGED'
          write(6,*)' Check the # of electron, mb1, and mb2'
          write(6,*)' EF1 & EF2=',ef1,ef2
          write(6,*)' octmp,zne=',octmp,zne,octmp-zne
          do s=1,MSP
             do k=1,MBZ
                write(*,*) "s,k =",s,k
                do n=mb1,mb2
                   write(*,*) n,occ(n,k,s),esp(n,k,s)
                end do
             end do
          end do
       end if
       stop'FERMI'
    end if

100 continue

    efermi   = ef
    Eentropy = 0.d0

    do s=1,MSP
    do k=1,MBZ
       do n=1,mb1-1
          occ(n,k,s)=2.d0*wbz(k)/dble(MSP)
       end do
       do n=mb1,mb2
          occ(n,k,s)=2.d0*occ(n,k,s)*wbz(k)/dble(MSP)
       end do
       do n=mb2+1,MB
          occ(n,k,s)=0.d0
       end do
    end do
    end do

    call write_border( 1, " calc_fermi_mp(end)" )

  END SUBROUTINE calc_fermi

  FUNCTION ff0(n,x)
    implicit none
    integer :: n,i
    real(8) :: x,ff,ff0,hp0,hp1,hp2,hp3
    ff0 = 0.5d0*(1.d0-bberf(x))
    if ( n <= 0 ) return
    hp0 = 1.d0
    hp1 = 2.d0*x
    ff  = factor(1)*hp1
    do i=2,n
       hp2 = 2.d0*x*hp1 - 2.d0*(2*i-3)*hp0
       hp3 = 2.d0*x*hp2 - 2.d0*(2*i-2)*hp1
       ff  = ff + factor(i)*hp3
       hp0 = hp2
       hp1 = hp3
    end do
    ff0 = ff0 + ff*exp(-x*x)
    return
  END FUNCTION ff0

END MODULE fermi_module
