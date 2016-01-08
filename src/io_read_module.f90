MODULE io_read_module

  use parallel_module
  use wf_module
  use rgrid_module
  use rgrid_mol_module, only: LL

  implicit none

  PRIVATE
  PUBLIC :: simple_wf_io_read

CONTAINS

  SUBROUTINE simple_wf_io_read( file_wf2, SYStype, IO_ctrl, DISP_SWITCH &
       , b, wf_out, occ_out, kbb_out )

    implicit none

    character(*),intent(IN) :: file_wf2
    integer,intent(IN) :: SYStype, IO_ctrl
    logical,intent(IN) :: DISP_SWITCH
    type(wfrange),optional,intent(INOUT) :: b
#ifdef _DRSDFT_
    real(8),allocatable,optional,intent(INOUT) :: wf_out(:,:,:,:)
#else
    complex(8),allocatable,optional,intent(INOUT) :: wf_out(:,:,:,:)
#endif
    real(8),allocatable,optional,intent(INOUT) :: occ_out(:,:,:)
    real(8),allocatable,optional,intent(INOUT) :: kbb_out(:,:)

    integer :: ierr,istatus(MPI_STATUS_SIZE,123),irank,jrank
    integer :: itmp(21),n,k,s,i,i1,i2,i3,j1,j2,j3,mx,my,mz,m1,m2,m3
    integer :: ML_tmp,ML1_tmp,ML2_tmp,ML3_tmp,MB_tmp,MB1_tmp,MB2_tmp
    integer :: ML,ML1,ML2,ML3,MB,n1,n2,ML0,MBZ_tmp,MSP_tmp,MMBZ_tmp
    integer :: MB_0,MB_1,MBZ,MBZ_0,MBZ_1,MSP,MSP_0,MSP_1
    integer :: IO_ctrl0, OC, TYPE_WF
    integer,allocatable :: LL2(:,:),LL_tmp(:,:),ir(:),id(:)
    real(8) :: aa0(3,3),bb0(3,3)
    real(8),allocatable :: kbb0(:,:),occ_tmp(:,:,:)
    logical :: flag_related, flag_newformat
    character(5) :: cmyrank
    character(32) :: file_wf_split
#ifdef _DRSDFT_
    integer,parameter :: TYPE_MAIN=MPI_REAL8
    real(8),parameter :: zero=0.0d0
    real(8),allocatable :: utmp(:), utmp3(:,:,:)
    real(4),allocatable :: utmpSP(:)
#else
    integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: utmp(:), utmp3(:,:,:)
    complex(4),allocatable :: utmpSP(:)
#endif
    real(8),allocatable :: dtmp(:)
    real(4),allocatable :: dtmpSP(:)
    type(para) :: pinfo0, pinfo1

    call write_border( 0, " simple_wf_io_read(start)" )

! ---

    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    MB  = MB_WF
    n1  = ML_0_WF
    n2  = ML_1_WF
    ML0 = n2 - n1 + 1

    MB_0  = MB_0_WF
    MB_1  = MB_1_WF
    MBZ   = MK_WF
    MBZ_0 = MK_0_WF
    MBZ_1 = MK_1_WF
    MSP   = MS_WF
    MSP_0 = MS_0_WF
    MSP_1 = MS_1_WF

! ---

    allocate( LL2(3,ML)    ) ; LL2=0
    allocate( LL_tmp(3,ML) ) ; LL_tmp=0

    if ( SYStype == 0 ) then

       i=n1-1
       do i3=Igrid(1,3),Igrid(2,3)
       do i2=Igrid(1,2),Igrid(2,2)
       do i1=Igrid(1,1),Igrid(2,1)
          i=i+1
          LL2(1,i)=i1
          LL2(2,i)=i2
          LL2(3,i)=i3
       end do
       end do
       end do

    else if ( SYStype == 1 ) then

       LL2(1:3,n1:n2) = LL(1:3,n1:n2)

    end if

    allocate( ir(0:np_grid-1), id(0:np_grid-1) )

    ir(0:np_grid-1)=3*ir_grid(0:np_grid-1)
    id(0:np_grid-1)=3*id_grid(0:np_grid-1)
    call mpi_allgatherv(LL2(1,n1),ir(myrank_g),mpi_integer &
         ,LL2,ir,id,mpi_integer,comm_grid,ierr)

    deallocate( id,ir )

! ---

    if ( myrank == 0 ) then

       itmp(:)=0

       open(3,file=file_wf2,form='unformatted')
       read(3) itmp(1)

       flag_newformat = .false.
       if ( itmp(1) < 0 ) flag_newformat = .true.
       write(*,*) "flag_newformat=",flag_newformat

       if ( flag_newformat ) then
          read(3) itmp(1:4)
          read(3) itmp(5:7)
          read(3) itmp(8:10)
          read(3) itmp(11:13)
!          read(3) itmp(14:21)
       else
          rewind 3
          read(3) itmp(1:4)
          read(3) itmp(5:7)
          itmp(8)=MBZ
          itmp(9)=MSP
       end if

    end if

    call mpi_bcast(flag_newformat,1,mpi_logical,0,mpi_comm_world,ierr)
    call mpi_bcast(itmp,size(itmp),mpi_integer,0,mpi_comm_world,ierr)
    ML_tmp  = itmp(1)
    ML1_tmp = itmp(2)
    ML2_tmp = itmp(3)
    ML3_tmp = itmp(4)
    MB_tmp  = itmp(5)
    MB1_tmp = itmp(6)
    MB2_tmp = itmp(7)
    MBZ_tmp = itmp(8)
    MSP_tmp = itmp(9)
    MMBZ_tmp= itmp(10)
    IO_ctrl0= itmp(11)
    OC      = itmp(12)
    TYPE_WF = itmp(13)

    if ( DISP_SWITCH ) then
       write(*,*) "ML ,ML_tmp                =",ML ,ML_tmp
       write(*,*) "ML1    ,ML2    ,ML3       =",ML1,ML2,ML3
       write(*,*) "ML1_tmp,ML2_tmp,ML3_tmp   =",ML1_tmp,ML2_tmp,ML3_tmp
       write(*,*) "MB ,MB_tmp                =",MB ,MB_tmp
       write(*,*) "MB1_tmp, MB2_tmp          =",MB1_tmp,MB2_tmp
       if ( flag_newformat ) then
          write(*,*) "MBZ_tmp, MSP_tmp,MMBZ_tmp =",MBZ_tmp,MSP_tmp,MMBZ_tmp
          write(*,*) "IO_ctrl0, OC, TYPE_WF     =",IO_ctrl0, OC, TYPE_WF
          write(*,*) "(new format)"
       else
          write(*,*) "(old format)"
       end if
    end if
    if ( ML_tmp /= ML .or. ML1_tmp /= ML1 .or. &
         ML2_tmp /= ML2 .or. ML3_tmp /= ML3 ) stop "stop@simple_wf_io_read"

    if ( myrank == 0 ) then
       read(3) LL_tmp(:,:)
    end if
    call mpi_bcast(LL_tmp,3*ML,mpi_integer,0,mpi_comm_world,ierr)

    if ( DISP_SWITCH ) then
       write(*,*) "minval(LL_tmp),maxval(LL_tmp)"
       write(*,*) minval(LL_tmp(1,1:ML) ),maxval( LL_tmp(1,1:ML))
       write(*,*) minval(LL_tmp(2,1:ML) ),maxval( LL_tmp(2,1:ML))
       write(*,*) minval(LL_tmp(3,1:ML) ),maxval( LL_tmp(3,1:ML))
    end if

    if ( flag_newformat ) then

       if ( present(b) ) then
          if ( b%MB == 0 ) then
             b%MB  = MB
             b%MB0 = 1
             b%MB1 = MB
          end if
          if ( b%MK == 0 ) then
             b%MK  = MBZ_tmp
             b%MK0 = 1
             b%MK1 = MBZ_tmp
             b%MMK = MMBZ_tmp
          end if
          if ( b%MS == 0 ) then
             b%MS  = MSP_tmp
             b%MS0 = 1
             b%MS1 = MSP_tmp
          end if
          if ( present(wf_out) ) call allocate_b_wf( b, wf_out )
          if ( present(occ_out) ) call allocate_b_occ( b, occ_out )
          if ( present(kbb_out) ) then
             allocate( kbb_out(3,b%MK0:b%MK1) ) ; kbb_out=0.0d0
          end if
       end if

    end if

! ---

    allocate( occ_tmp(MB_tmp,MBZ_tmp,MSP_tmp) ) ; occ_tmp=0.0d0
    if ( myrank == 0 ) read(3) occ_tmp
    call mpi_bcast(occ_tmp,size(occ_tmp),mpi_real8,0,mpi_comm_world,ierr)
    
    if ( present(occ_out) ) then
       m1=b%MB0
       m2=min(b%MB1,MB_tmp)
       occ_out(m1:m2,b%MK0:b%MK1,b%MS0:b%MS1) &
            = occ_tmp(m1:m2,b%MK0:b%MK1,b%MS0:b%MS1)
    else
       m1=min(MB ,MB_tmp )
       m2=min(MBZ,MBZ_tmp)
       m3=min(MSP,MSP_tmp)
       occ(1:m1,1:m2,1:m3)=occ_tmp(1:m1,1:m2,1:m3)
       if ( myrank == 0 ) then
          write(*,*) "sum(occ)=",sum(occ)
       end if
    end if

    deallocate( occ_tmp )

! ---

    if ( flag_newformat ) then

       if ( present(kbb_out) ) then
          if ( myrank == 0 ) read(3) aa0,bb0,kbb_out
          call mpi_bcast(kbb_out,size(kbb_out),mpi_real8,0,mpi_comm_world,ierr)
       else
          if ( myrank == 0 ) read(3)
       end if

       if ( myrank == 0 ) then
          read(3) itmp(14:21)
          read(3)
          read(3)
          read(3)
          read(3)
          write(*,'(1x,"np(0:7)=",i5,1x,7i4)') itmp(14:21)
       end if
       call mpi_bcast(itmp(14),8,mpi_integer,0,mpi_comm_world,ierr)

       pinfo0%np(0:7) = itmp(14:21) 
       call construct_para( ML_tmp, MB_tmp, MBZ_tmp, MSP_tmp, pinfo0 )

       call get_np_parallel( pinfo1%np )
       call construct_para( ML, MB, MBZ, MSP, pinfo1 )

    end if

! ---

    if ( SYStype == 0 ) then

       allocate( utmp3(0:ML1-1,0:ML2-1,0:ML3-1) ) ; utmp3=zero

    else if ( SYStype == 1 ) then

       mx=(Ngrid(1)-1)/2
       my=(Ngrid(2)-1)/2
       mz=(Ngrid(3)-1)/2
       allocate( utmp3(-mx:mx,-my:my,-mz:mz) ) ; utmp3=0.d0

    end if

    if ( IO_ctrl == 3 ) then

       if ( myrank == 0 ) close(3)

       write(cmyrank,'(i5.5)') myrank
       file_wf_split = trim(file_wf2)//"."//trim(adjustl(cmyrank))

       open(3,file=file_wf_split,form="unformatted")

    end if

    allocate( utmp(ML) ) ; utmp=zero
    if ( OC == 4 .or. OC == 5 ) then
       allocate( utmpSP(ML) ) ; utmpSP=zero
    end if

    if ( type_wf == 1 ) then

       allocate( dtmp(ML) ) ; dtmp=0.0d0
       if ( OC == 4 .or. OC == 5 ) then
          allocate( dtmpSP(ML) ) ; dtmpSP=0.0d0
       end if

    end if

    do s=1,MSP_tmp
    do k=1,MBZ_tmp
    do n=MB1_tmp,MB2_tmp

       if ( IO_ctrl == 0 ) then

          call chk_rank_relevance( n,k,s, irank, flag_related )

          if ( myrank == 0 ) then

             select case( OC )
             case default
                if ( type_wf == 1 ) then
                   read(3) dtmp
                   utmp=dtmp
                else
                   read(3) utmp
                end if
             case(4,5)
                if ( type_wf == 1 ) then
                   read(3) dtmpSP
                   utmp=dtmpSP
                else
                   read(3) utmpSP
                   utmp=utmpSP
                end if
             end select

             do i=1,ML
                i1=LL_tmp(1,i) ; i2=LL_tmp(2,i) ; i3=LL_tmp(3,i)
                utmp3(i1,i2,i3)=utmp(i)
             end do
             do i=1,ML
                i1=LL2(1,i) ; i2=LL2(2,i) ; i3=LL2(3,i)
                utmp(i)=utmp3(i1,i2,i3)
             end do

          end if ![ myrank == 0 ]

          call mpi_barrier(mpi_comm_world,ierr)

          if ( irank /= 0 .and. myrank_f == 0 ) then

             if ( myrank == 0 ) then
                call mpi_send(utmp,ML,TYPE_MAIN,irank,0, &
                     mpi_comm_world,ierr)
             end if
             if ( myrank == irank ) then
                call mpi_recv(utmp,ML,TYPE_MAIN,0,0, &
                     mpi_comm_world,istatus,ierr)
             end if

          end if

          call mpi_barrier(mpi_comm_world,ierr)

          if ( flag_related ) then
             if ( present(wf_out) ) then
                call mpi_scatterv(utmp,ir_grid,id_grid,TYPE_MAIN &
                     ,wf_out(n1,n,k,s),ML0,TYPE_MAIN,0,comm_grid,ierr)
             else
                call mpi_scatterv(utmp,ir_grid,id_grid,TYPE_MAIN &
                     ,unk(n1,n,k,s),ML0,TYPE_MAIN,0,comm_grid,ierr)
             end if
          end if

       else if ( IO_ctrl == 3 ) then

          call get_corresponding_rank( myrank_g, n,k,s, pinfo1, irank )
          call get_corresponding_rank( myrank_g, n,k,s, pinfo0, jrank )

          if ( irank == myrank .or. jrank == myrank ) then

             select case(OC)
             case default

                if ( type_wf == 1 ) then

                   if ( jrank == myrank ) then
                      read(3) dtmp(n1:n2)
                   end if
                   if ( irank == jrank ) then
                      if ( present(wf_out) ) then
                         wf_out(n1:n2,n,k,s)=dtmp(n1:n2)
                      else
                         unk(n1:n2,n,k,s)=dtmp(n1:n2)
                      end if
                   else if ( irank == myrank ) then
                      call mpi_recv(dtmp(n1),n2-n1+1,MPI_REAL8,jrank,0, &
                           mpi_comm_world,istatus,ierr)
                      if ( present(wf_out) ) then
                         wf_out(n1:n2,n,k,s)=dtmp(n1:n2)
                      else
                         unk(n1:n2,n,k,s)=dtmp(n1:n2)
                      end if
                   else if ( jrank == myrank ) then
                      call mpi_send(dtmp(n1),n2-n1+1,MPI_REAL8,irank,0, &
                           mpi_comm_world,ierr)
                   end if

                else

                   if ( jrank == myrank ) then
                      read(3) utmp(n1:n2)
                   end if
                   if ( irank == jrank ) then
                      if ( present(wf_out) ) then
                         wf_out(n1:n2,n,k,s)=utmp(n1:n2)
                      else
                         unk(n1:n2,n,k,s)=utmp(n1:n2)
                      end if
                   else if ( irank == myrank ) then
                      call mpi_recv(utmp(n1),n2-n1+1,TYPE_MAIN,jrank,0, &
                           mpi_comm_world,istatus,ierr)
                      if ( present(wf_out) ) then
                         wf_out(n1:n2,n,k,s)=utmp(n1:n2)
                      else
                         unk(n1:n2,n,k,s)=utmp(n1:n2)
                      end if
                   else if ( jrank == myrank ) then
                      call mpi_send(utmp(n1),n2-n1+1,TYPE_MAIN,irank,0, &
                           mpi_comm_world,ierr)
                   end if

                end if

             case(4,5)

                if ( type_wf == 1 ) then

                   if ( jrank == myrank ) then
                      read(3) dtmpSP(n1:n2)
                   end if
                   if ( irank == jrank ) then
                      if ( present(wf_out) ) then
                         wf_out(n1:n2,n,k,s)=dtmpSP(n1:n2)
                      else
                         unk(n1:n2,n,k,s)=dtmpSP(n1:n2)
                      end if
                   else if ( irank == myrank ) then
                      call mpi_recv(dtmpSP(n1),n2-n1+1,MPI_REAL,jrank,0, &
                           mpi_comm_world,istatus,ierr)
                      if ( present(wf_out) ) then
                         wf_out(n1:n2,n,k,s)=dtmpSP(n1:n2)
                      else
                         unk(n1:n2,n,k,s)=dtmpSP(n1:n2)
                      end if
                   else if ( jrank == myrank ) then
                      call mpi_send(dtmpSP(n1),n2-n1+1,MPI_REAL,irank,0, &
                           mpi_comm_world,ierr)
                   end if

                else

                   if ( jrank == myrank ) then
                      read(3) utmpSP(n1:n2)
                   end if
                   if ( irank == jrank ) then
                      if ( present(wf_out) ) then
                         wf_out(n1:n2,n,k,s)=utmpSP(n1:n2)
                      else
                         unk(n1:n2,n,k,s)=utmpSP(n1:n2)
                      end if
                   else if ( irank == myrank ) then
                      call mpi_recv(utmpSP(n1),n2-n1+1,TYPE_MAIN,jrank,0, &
                           mpi_comm_world,istatus,ierr)
                      if ( present(wf_out) ) then
                         wf_out(n1:n2,n,k,s)=utmpSP(n1:n2)
                      else
                         unk(n1:n2,n,k,s)=utmpSP(n1:n2)
                      end if
                   else if ( jrank == myrank ) then
                      call mpi_send(utmpSP(n1),n2-n1+1,TYPE_MAIN,irank,0, &
                           mpi_comm_world,ierr)
                   end if

                end if

             end select

          end if

       end if

    end do ! n
    end do ! k
    end do ! s

    if ( IO_ctrl == 0 .and. myrank == 0 .or. IO_ctrl == 3 ) close(3)

    if ( allocated(dtmpSP) ) deallocate( dtmpSP )
    if ( allocated(dtmp) ) deallocate( dtmp )
    if ( allocated(utmpSP) ) deallocate( utmpSP )
    deallocate( utmp )
    deallocate( utmp3 )
    deallocate( LL_tmp )
    deallocate( LL2 )

#ifdef _DRSDFT_
    call mpi_bcast( unk,size(unk),MPI_REAL8,0,comm_fkmb,ierr )
#else
    call mpi_bcast( unk,size(unk),MPI_COMPLEX16,0,comm_fkmb,ierr )
#endif

    if ( DISP_SWITCH ) then
       write(*,*) "read from ",file_wf2
    end if

    call write_border( 0, " simple_wf_io_read(end)" )

    return

  END SUBROUTINE simple_wf_io_read

  SUBROUTINE chk_rank_relevance( n,k,s, irank, flag_related )
    implicit none
    integer,intent(IN) :: n,k,s
    integer,intent(OUT) :: irank
    logical,intent(OUT) :: flag_related
    integer :: i1,i2,i3,j1,j2,j3
    do irank=0,nprocs-1
       i1=id_band(id_class(irank,4))+1
       j1=id_band(id_class(irank,4))+ir_band(id_class(irank,4))
       i2=id_bzsm(id_class(irank,5))+1
       j2=id_bzsm(id_class(irank,5))+ir_bzsm(id_class(irank,5))
       i3=id_spin(id_class(irank,6))+1
       j3=id_spin(id_class(irank,6))+ir_spin(id_class(irank,6))
       if ( id_grid(id_class(irank,0))==0 .and. &
            i1<=n .and. n<=j1 .and. &
            i2<=k .and. k<=j2 .and. &
            i3<=s .and. s<=j3 ) exit
    end do
    if ( irank >= nprocs ) then
       write(*,*) "ERROR!(stpo@simple_wf_io_read): myrank=",myrank
       call stop_program_f( "stop@chk_rank_relevance" )
    end if
    flag_related = .false.
    if ( MK_0_WF <= k .and. k <= MK_1_WF .and. &
         MB_0_WF <= n .and. n <= MB_1_WF  .and. &
         MS_0_WF <= s .and. s <= MS_1_WF ) flag_related = .true.
  END SUBROUTINE chk_rank_relevance

  SUBROUTINE get_corresponding_rank( mrank_g, n,k,s, p, irank )
    implicit none
    integer,intent(IN) :: mrank_g, n,k,s
    type(para),intent(INOUT) :: p
    integer,intent(OUT) :: irank
    integer :: i1,i2,i3,i4,i5,i6,i7
    integer :: n0,n1,k0,k1,s0,s1,irank_g
    irank=-1
    loop: do i7=0,p%np(7)-1
    do i6=0,p%np(6)-1
    do i5=0,p%np(5)-1
    do i4=0,p%np(4)-1
       s0=p%spin%id(i6)+1 ; s1=p%spin%id(i6)+p%spin%ir(i6)
       k0=p%bzsm%id(i5)+1 ; k1=p%bzsm%id(i5)+p%bzsm%ir(i5)
       n0=p%band%id(i4)+1 ; n1=p%band%id(i4)+p%band%ir(i4)
       irank_g=-1
       do i3=0,p%np(3)-1
       do i2=0,p%np(2)-1
       do i1=0,p%np(1)-1
          irank=irank+1
          irank_g=irank_g+1
          if ( s0 <= s .and. s <= s1 .and. &
               k0 <= k .and. k <= k1 .and. &
               n0 <= n .and. n <= n1 .and. irank_g == mrank_g ) exit loop
       end do
       end do
       end do
    end do
    end do
    end do
    end do loop
  END SUBROUTINE get_corresponding_rank

END MODULE io_read_module
