MODULE wf_data_module

  use parallel_module
  use watch_module
  use atom_module

  implicit none

  PRIVATE
  PUBLIC :: wf_data, read_parameters, Nspin, Nband, Ngrid, Nbzsm &
       ,lst_wf,nlst,lst_wf_tmp

  integer,parameter :: unit=1
  integer,parameter :: urwf1=3
  integer,parameter :: urwf2=80
  integer,parameter :: uwwf=81
  integer,parameter :: ucwf=20
  integer,parameter :: ugwf=22

  integer :: IO_CTRL=3

  character(30) :: file_wf_r = 'wf.dat1'
  character(30) :: file_wf_w = 'wf.dat2'
  character(30) :: file_cube = 'wfcub_'

  real(8) :: ax,aa(3,3),Vaa,dV
  integer :: Ngrid(0:3),Nband,Nbzsm,Nspin

  integer :: nlst, lst_wf_tmp(4)
  integer,allocatable :: lst_wf(:)
  integer,allocatable :: LL(:,:)
  integer :: zatom(9)
  integer :: node_partition0(6)
  integer :: iswitch_post_wf

  real(8),allocatable :: occ(:,:,:)

#ifdef _DRSDFT_
  integer,parameter :: TYPE_MAIN=MPI_REAL8
  real(8),allocatable :: utmp(:),wf(:,:)
#else
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
  complex(8),allocatable :: utmp(:),wf(:,:)
#endif

CONTAINS


  SUBROUTINE read_parameters
    implicit none
    character(10) :: cbuf,ckey
    integer :: i
    ax=0.0d0
    aa(:,:)=0.0d0
    Ngrid(:)=0
    NBAND=0
    NSPIN=1
    NBZSM=1
    node_partition(:)=1
    node_partition0(:)=1
    lst_wf_tmp(:)=0
    zatom(:)=0
    iswitch_post_wf=0
    if ( myrank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey(1:2) == "AX" ) then
             backspace(unit)
             read(unit,*) cbuf,ax
          else if ( ckey(1:2) == "A1" ) then
             backspace(unit)
             read(unit,*) cbuf,aa(1:3,1)
          else if ( ckey(1:2) == "A2" ) then
             backspace(unit)
             read(unit,*) cbuf,aa(1:3,2)
          else if ( ckey(1:2) == "A3" ) then
             backspace(unit)
             read(unit,*) cbuf,aa(1:3,3)
          else if ( ckey(1:5) == "NBAND" ) then
             backspace(unit)
             read(unit,*) cbuf,Nband
          else if ( ckey(1:5) == "NSPIN" ) then
             backspace(unit)
             read(unit,*) cbuf,Nspin
          else if ( ckey(1:4) == "NPBZ" ) then
             backspace(unit)
             read(unit,*) cbuf,Nbzsm
          else if ( ckey(1:5) == "NGRID" ) then
             backspace(unit)
             read(unit,*) cbuf,Ngrid(1:3)
          else if ( ckey(1:6) == "PROCS" ) then
             backspace(unit)
             read(unit,*) cbuf,node_partition(1:6)
          else if ( ckey(1:6) == "PROCS0" ) then
             backspace(unit)
             read(unit,*) cbuf,node_partition0(1:6)
          else if ( ckey(1:5) == "LSTWF" ) then
             backspace(unit)
             read(unit,*) cbuf,lst_wf_tmp(1:4)
          else if ( ckey(1:5) == "ZATOM" ) then
             backspace(unit)
             read(unit,*) cbuf,zatom(:)
          else if ( ckey(1:8) == "SWPOSTWF" ) then
             backspace(unit)
             read(unit,*) cbuf,iswitch_post_wf
          end if
       end do
999    continue
       Ngrid(0)=Ngrid(1)*Ngrid(2)*Ngrid(3)
       write(*,'(1x,"Ngrid=",i8,2x,3i4)') Ngrid(:)
       write(*,'(1x,"Nband=",i8)') Nband
       write(*,'(1x,"Nzsm,Nspin=",2i4)') Nbzsm,Nspin
       write(*,'(1x,"node_partition =",6i4)') node_partition(:)
       write(*,'(1x,"node_partition0=",6i4)') node_partition0(:)
       write(*,'(1x,"lst_wf_tmp=",6i8)') lst_wf_tmp(:)
       write(*,'(1x,"ax=",f20.15)') ax
       write(*,'(1x,"a1=",3f20.15)') aa(1:3,1)
       write(*,'(1x,"a2=",3f20.15)') aa(1:3,2)
       write(*,'(1x,"a3=",3f20.15)') aa(1:3,3)
       write(*,'(1x,"ZATOM=",5i3)') zatom(1:5)
       write(*,'(1x,"iswitch_post_wf=",i3)') iswitch_post_wf
    end if
    call send_parameters
  END SUBROUTINE read_parameters

  SUBROUTINE send_parameters
    implicit none
    integer :: ierr
    call mpi_bcast(Ngrid,4,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(Nband,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(Nspin,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(Nbzsm,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(node_partition ,6,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(node_partition0,6,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(lst_wf_tmp,4,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(iswitch_post_wf,1,mpi_integer,0,mpi_comm_world,ierr)
  END SUBROUTINE send_parameters


  SUBROUTINE wf_data
    implicit none
    character(5) :: cmyrank
    character(32) :: file_wf_split,fname
    integer :: Nb_1,Nb_2,ML,ML1,ML2,ML3,m1,m2,m3
    integer :: n,k,s,l,ierr,i,j,i1,i2,i3,j1,j2
    integer :: ML_0,ML_1,MB_0,MB_1,MBZ_0,MBZ_1,MSP_0,MSP_1
    integer :: nprocs0,myrank0,iii,n_rank,ML0_0,ML0_1
    integer,allocatable :: ir(:),irank0(:)
    integer,allocatable :: ir_grid0(:),id_grid0(:)
    real(8),allocatable :: utmp_wf(:,:,:),vtmp(:,:,:)
    real(8) :: ct(0:5),et(0:5),rr(3)

    ML    = sum(ir_grid)
    ML_0  = id_grid(myrank_g)+1
    ML_1  = id_grid(myrank_g)+ir_grid(myrank_g)
    MB_0  = id_band(myrank_b)+1
    MB_1  = id_band(myrank_b)+ir_band(myrank_b)
    MBZ_0 = id_bzsm(myrank_k)+1
    MBZ_1 = id_bzsm(myrank_k)+ir_bzsm(myrank_k)
    MSP_0 = id_spin(myrank_s)+1
    MSP_1 = id_spin(myrank_s)+ir_spin(myrank_s)

    rr(:)=0.0d0

    nprocs0=1
    do i=1,6
       nprocs0=nprocs0*node_partition0(i)
    end do

    allocate( ir(0:nprocs-1) ) ; ir=0
    do i=0,nprocs0-1
       j=mod(i,nprocs)
       ir(j)=ir(j)+1
    end do
    n_rank=ir(myrank)
    allocate( irank0(n_rank) ) ; irank0=0
    i1=sum( ir(0:myrank) )-ir(myrank)
    i2=i1+ir(myrank)-1
    deallocate( ir )
    n=0
    do i=i1,i2
       n=n+1
       irank0(n)=i
    end do

!    write(*,'(1x,2i4,2x,10i4)') myrank,n_rank,irank0(1:n_rank)

    allocate( id_grid0(0:nprocs0-1) ) ; id_grid0=0
    allocate( ir_grid0(0:nprocs0-1) ) ; ir_grid0=0

    do i=1,ML
       j=mod(i-1,nprocs0)
       ir_grid0(j)=ir_grid0(j)+1
    end do
    do j=0,nprocs0-1
       id_grid0(j)=sum( ir_grid0(0:j) )-ir_grid0(j)
    end do
    if ( myrank == 0 ) then
       do j=0,nprocs0-1
          write(*,*) j,id_grid0(j)+1,id_grid0(j)+ir_grid0(j),ir_grid0(j)
       end do
    end if

! ---

    aa(:,:)=ax*aa(:,:)

    Vaa = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
         +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
         -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)

    dV = abs(Vaa)/ML

    if ( myrank == 0 ) then
       if ( TYPE_MAIN == MPI_REAL8 ) then
          write(*,*) "size(wf)(MB)r", 8.0d0*ML*nlst/1024.d0**2
       else
          write(*,*) "size(wf)(MB)c=",16.0d0*ML*nlst/1024.d0**2
       end if
    end if

    allocate( wf(ML,nlst) ) ; wf=(0.0d0,0.0d0)

    select case( IO_CTRL )
    case default

    case( 3 )

       call watch(ct(0),et(0))

       open(urwf1,file=file_wf_r,form='unformatted')
       read(urwf1) Ngrid(0:3)
       read(urwf1) Nband,Nb_1,Nb_2

       allocate( LL(3,Ngrid(0)) ) ; LL=0
       allocate( occ(Nband,Nbzsm,Nspin) ) ; occ=0.0d0

       read(urwf1) LL
       read(urwf1) occ

       close(urwf1)

       call watch(ct(1),et(1))
       if ( myrank == 0 ) then
          write(*,*) "time(wf_data(1))",ct(1)-ct(0),et(1)-et(0)
       end if

       select case( iswitch_post_wf )
       case( 0 )

          call watch(ct(1),et(1))

          do iii=1,n_rank

             myrank0 = irank0(iii)
             ML0_0   = id_grid0(myrank0)+1
             ML0_1   = id_grid0(myrank0)+ir_grid0(myrank0)

!write(*,'(1x,10i6)') myrank,iii,myrank0,ML0_0,ML0_1,Nband,Nbzsm,Nspin

             write(cmyrank,'(i5.5)') myrank0
             file_wf_split = trim(file_wf_r)//"."//trim(adjustl(cmyrank))

             allocate( utmp(ML0_0:ML0_1) ) ; utmp=(0.0d0,0.0d0)

             open(urwf2,file=file_wf_split,form='unformatted')

             do s=1,Nspin
             do k=1,Nbzsm
             do n=1,Nband

                read(urwf2) utmp

                do l=1,nlst
                   if ( lst_wf(l) == n ) then
                      wf(ML0_0:ML0_1,l) = utmp(:)
                      exit
                   end if
                end do

             end do ! n
             end do ! k
             end do ! s

             close(urwf2)

             deallocate( utmp )

          end do ! iii

          call watch(ct(2),et(2))
          if ( myrank == 0 ) then
             write(*,*) "time(wf_data(2))",ct(2)-ct(1),et(2)-et(1)
          end if

          do l=1,nlst
             call mpi_allgatherv(wf(ML_0,l),ir_grid(myrank_g),TYPE_MAIN &
                  ,wf(1,l),ir_grid,id_grid,TYPE_MAIN,comm_grid,ierr)
          end do

          call watch(ct(3),et(3))
          if ( myrank == 0 ) then
             write(*,*) "time(wf_data(3))",ct(3)-ct(2),et(3)-et(2)
          end if

          if ( myrank == 0 ) then
             ML1=Ngrid(1)
             ML2=Ngrid(2)
             ML3=Ngrid(3)
             allocate( utmp_wf(ML1,ML2,ML3) ) ; utmp_wf=0.0d0
             do l=1,nlst
                do i=1,ML
                   i1=LL(1,i)+1
                   i2=LL(2,i)+1
                   i3=LL(3,i)+1
                   utmp_wf(i1,i2,i3)=abs(wf(i,l))**2
                end do
goto 10
                m1=ML1/2
                m2=ML2/2
                m3=ML3/2
                allocate( vtmp(-ML1:ML1-1,-ML2:ML2-1,0:ML3-1) )
                vtmp=0.0d0
                do i3=0,ML3-1
                do i2=-ML2,ML2-1
                do i1=-ML1,ML1-1
                   j1=mod(i1+ML1,ML1)
                   j2=mod(i2+ML2,ML2)
                   vtmp(i1,i2,i3)=utmp_wf(j1+1,j2+1,i3+1)
                end do
                end do
                end do
                do i3=0,ML3-1
                do i2=-m2+1,m2
                do i1=-m1+1,m1
                   utmp_wf(i1+m1,i2+m2,i3+1)=vtmp(i1,i2,i3)
                end do
                end do
                end do
                deallocate( vtmp )
                rr(1)=(-m1+1)*aa(1,1)/ML1
                rr(2)=(-m2+1)*aa(2,2)/ML2
                rr(3)=0.0d0
10 continue
                write(cmyrank,'(i5.5)') lst_wf(l)
                fname=trim(file_cube)//trim(adjustl(cmyrank))//".cub"
                open(ucwf,file=fname)
                call cube_gen(ML1,ML2,ML3,utmp_wf,rr)
                close(ucwf)
                fname="wfgp_"//trim(adjustl(cmyrank))
                open(ugwf,file=fname)
                call gpdat_gen(ML1,ML2,ML3,utmp_wf)
                close(ugwf)
             end do
             deallocate( utmp_wf )
          end if

          call watch(ct(4),et(4))
          if ( myrank == 0 ) then
             write(*,*) "time(wf_data(4))",ct(4)-ct(3),et(4)-et(3)
          end if
goto 20
          if ( myrank_g == 0 ) then
             open(uwwf,file=file_wf_w,form='unformatted')
             write(uwwf) TYPE_MAIN,MPI_REAL8,MPI_COMPLEX16
             write(uwwf) Ngrid(0:3)
             write(uwwf) LL
             write(uwwf) nlst
             write(uwwf) lst_wf
             do l=1,nlst
                write(uwwf) wf(:,l)
             end do
             close(uwwf)
          end if
20 continue

          call watch(ct(5),et(5))
          if ( myrank == 0 ) then
             write(*,*) "time(wf_data(5))",ct(5)-ct(4),et(5)-et(4)
          end if

       case( 1 )

          call wf_data_sub1( n_rank, nprocs0, irank0, id_grid0, ir_grid0 )

       end select ! iswitch_post_wf

       deallocate( occ )
       deallocate( LL  )
       deallocate( ir_grid0,id_grid0 )
       deallocate( irank0 )

    end select ! IO_CTRL

    deallocate( wf )

  END SUBROUTINE wf_data


  SUBROUTINE wf_data_sub1( n_rank, nprocs0, irank0, id_grid0, ir_grid0 )
    implicit none
    integer,intent(IN) :: n_rank, nprocs0
    integer,intent(IN) :: irank0(n_rank), id_grid0(0:nprocs0-1)
    integer,intent(IN) :: ir_grid0(0:nprocs0-1)
    integer :: iii, myrank0,ML0_0,ML0_1,s,k,n,l,ierr
    integer :: m1,m2,m3,ML1,Ml2,ML3,i1,i2,i3,i,ML,mL_0
    character(5) :: cmyrank
    character(32) :: file_wf_split,fname
    real(8),allocatable :: rtmp(:), utmp_wf(:,:,:)
    real(8) :: ct(5),et(5),rr(3),c,d

    ML_0 = id_grid(myrank_g)+1
    ML   = Ngrid(0)
    allocate( rtmp(ML) ) ; rtmp=0.0d0

    rr(:) = 0.0d0

    call watch(ct(1),et(1))

    do iii=1,n_rank

       myrank0 = irank0(iii)
       ML0_0   = id_grid0(myrank0)+1
       ML0_1   = id_grid0(myrank0)+ir_grid0(myrank0)

       write(cmyrank,'(i5.5)') myrank0
       file_wf_split = trim(file_wf_r)//"."//trim(adjustl(cmyrank))

       write(*,*) iii,myrank0,ML0_0,ML0_1
       write(*,*) file_wf_split

       allocate( utmp(ML0_0:ML0_1) ) ; utmp=(0.0d0,0.0d0)

       open(urwf2,file=file_wf_split,form='unformatted')

       do s=1,Nspin
       do k=1,Nbzsm
       do n=1,Nband

          read(urwf2) utmp

          if ( lst_wf(1) <= n .and. n <= lst_wf(nlst) ) then
             rtmp(ML0_0:ML0_1) = rtmp(ML0_0:ML0_1) + abs( utmp(:) )**2
             write(*,*) n
          end if

       end do ! n
       end do ! k
       end do ! s

       close(urwf2)

       deallocate( utmp )

    end do ! iii

    call watch(ct(2),et(2))
    if ( myrank == 0 ) then
       write(*,*) "time(wf_data(2))",ct(2)-ct(1),et(2)-et(1)
    end if

    call mpi_allgatherv(rtmp(ML_0),ir_grid(myrank_g),MPI_REAL8 &
         ,rtmp,ir_grid,id_grid,MPI_REAL8,comm_grid,ierr)

    if ( myrank == 0 ) then
       write(*,*) "sum(rtmp)=",sum(rtmp)*dV
    end if

    call watch(ct(3),et(3))
    if ( myrank == 0 ) then
       write(*,*) "time(wf_data(3))",ct(3)-ct(2),et(3)-et(2)
    end if

    if ( myrank == 0 ) then
       ML1=Ngrid(1)
       ML2=Ngrid(2)
       ML3=Ngrid(3)
       allocate( utmp_wf(ML1,ML2,ML3) ) ; utmp_wf=0.0d0
       do i=1,ML
          i1=LL(1,i)+1
          i2=LL(2,i)+1
          i3=LL(3,i)+1
          utmp_wf(i1,i2,i3) = rtmp(i)
       end do
       fname="chrg.cub"
       open(ucwf,file=fname)
       call cube_gen(ML1,ML2,ML3,utmp_wf,rr)
       close(ucwf)
       fname="chrgp"
       open(ugwf,file=fname)
       call gpdat_gen(ML1,ML2,ML3,utmp_wf)
       close(ugwf)
       deallocate( utmp_wf )
    end if

    call watch(ct(4),et(4))
    if ( myrank == 0 ) then
       write(*,*) "time(wf_data(4))",ct(4)-ct(3),et(4)-et(3)
    end if

    deallocate( rtmp )

  END SUBROUTINE wf_data_sub1


  SUBROUTINE cube_gen(ML1,ML2,ML3,utmp_wf,rr)
    implicit none
    integer,intent(IN) :: ML1,ML2,ML3
    real(8),intent(IN) :: utmp_wf(ML1,ML2,ML3),rr(3)
    integer :: j,i1,i2,i3

    write(ucwf,*) "SYS1"
    write(ucwf,*) "SYS1"
    write(ucwf,'(i5,3f14.8)') Natom,rr(1:3)
    write(ucwf,'(i5,3f14.8)') ML1,(aa(j,1)/dble(ML1),j=1,3)
    write(ucwf,'(i5,3f14.8)') ML2,(aa(j,2)/dble(ML2),j=1,3)
    write(ucwf,'(i5,3f14.8)') ML3,(aa(j,3)/dble(ML3),j=1,3)
    do j=1,Natom
       write(ucwf,'(2i5,3f14.8)') zatom(ki_atom(j)),zatom(ki_atom(j)) &
            & ,sum(aa_atom(1:3,j)*aa(1,1:3)) &
            & ,sum(aa_atom(1:3,j)*aa(2,1:3)) &
            & ,sum(aa_atom(1:3,j)*aa(3,1:3))
    end do

    do i1=1,ML1
    do i2=1,ML2
!    do i3=1,ML3,6
  !     if ( ML3-i3 < 6 ) then
  !        select case ( mod(ML3,6) )
  !        case(1)
  !           write(20,'(1f14.8)')  utmp_wf(i1,i2,i3)
  !        case(2)
  !           write(20,'(2f14.8)') (utmp_wf(i1,i2,i3+j),j=0,1)
  !        case(3)
  !           write(20,'(3f14.8)') (utmp_wf(i1,i2,i3+j),j=0,2)
  !        case(4)
  !           write(20,'(4f14.8)') (utmp_wf(i1,i2,i3+j),j=0,3)
  !        case(5)
  !           write(20,'(5f14.8)') (utmp_wf(i1,i2,i3+j),j=0,4)
  !        end select
  !     else
  !        write(20,'(6f14.8)') (utmp_wf(i1,i2,i3+j),j=0,5)
  !     end if
       write(ucwf,'(6f14.8)') utmp_wf(i1,i2,:)
!    end do
    end do
    end do

  END SUBROUTINE cube_gen


  SUBROUTINE gpdat_gen(ML1,ML2,ML3,w)
     implicit none
     integer,intent(IN) :: ML1,ML2,ML3
     real(8),intent(IN) :: w(ML1,ML2,ML3)
     real(8) :: Hgrid(3)
     integer :: i
     Hgrid(1)=sqrt(sum(aa(:,1)**2))/ML1
     Hgrid(2)=sqrt(sum(aa(:,2)**2))/ML2
     Hgrid(3)=sqrt(sum(aa(:,3)**2))/ML3
     do i=1,ML1
        write(ugwf,'(1x,4g20.10)') (i-1)*Hgrid(1),w(i,1,1),sum(w(i,:,:))
     end do
     write(ugwf,*)
     write(ugwf,*)
     do i=1,ML2
        write(ugwf,'(1x,4g20.10)') (i-1)*Hgrid(2),w(1,i,1),sum(w(:,i,:))
     end do
     write(ugwf,*)
     write(ugwf,*)
     do i=1,ML3
        write(ugwf,'(1x,4g20.10)') (i-1)*Hgrid(3),w(1,1,i),sum(w(:,:,i))
     end do
  END SUBROUTINE gpdat_gen


END MODULE wf_data_module
