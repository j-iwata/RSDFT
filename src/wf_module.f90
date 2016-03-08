MODULE wf_module

  use parallel_module
  use wf_sub_module

  implicit none

  PRIVATE
  PUBLIC :: unk,esp,esp0,occ,res,init_wf,test_on_wf,gather_wf,gather_b_wf &
           ,ML_WF, ML_0_WF, ML_1_WF, MB_WF, MB_0_WF, MB_1_WF &
           ,MK_WF, MK_0_WF, MK_1_WF, MS_WF, MS_0_WF, MS_1_WF &
           ,Sunk &
           ,write_wf &
           ,hunk, read_wf, iflag_hunk, workwf &
           ,allocate_work_wf, deallocate_work_wf
  PUBLIC :: write_esp_wf
  PUBLIC :: write_info_esp_wf
  PUBLIC :: wfrange
  PUBLIC :: allocate_b_wf, allocate_b_occ

#ifdef _DRSDFT_
  real(8),parameter :: zero=0.d0
  real(8),allocatable :: unk(:,:,:,:)
  real(8),allocatable :: hunk(:,:,:,:)
  real(8),allocatable :: workwf(:,:)
  integer,parameter :: TYPE_MAIN=MPI_REAL8
#else
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),allocatable :: unk(:,:,:,:)
  complex(8),allocatable :: hunk(:,:,:,:)
  complex(8),allocatable :: workwf(:,:)
  integer,parameter :: TYPE_MAIN=MPI_COMPLEX16
#endif

#ifdef _DRSDFT_
    real(8),allocatable :: Sunk(:,:)
#else
    complex(8),allocatable :: Sunk(:,:)
#endif

  real(8),allocatable :: esp(:,:,:), esp0(:,:,:)
  real(8),allocatable :: occ(:,:,:)
  real(8),allocatable :: res(:,:,:)

  integer :: ML_WF, ML_0_WF, ML_1_WF
  integer :: MB_WF, MB_0_WF, MB_1_WF
  integer :: MK_WF, MK_0_WF, MK_1_WF
  integer :: MS_WF, MS_0_WF, MS_1_WF

  integer :: iwork_wf=0
  integer :: iflag_hunk=0

  type wfrange
     integer :: ML,ML0,ML1
     integer :: MB,MB0,MB1
     integer :: MK,MK0,MK1,MMK
     integer :: MS,MS0,MS1
  end type wfrange

CONTAINS


  SUBROUTINE read_wf( rank, unit )
    implicit none
    integer,intent(IN) :: rank,unit
    integer :: i,ierr
    character(6) :: cbuf,ckey
    call write_border( 0, " read_wf(start)" )
    if ( rank == 0 ) then
       rewind unit
       do i=1,10000
          read(unit,*,END=999) cbuf
          call convert_capital(cbuf,ckey)
          if ( ckey == "WORKWF" ) then
             backspace(unit)
             read(unit,*) cbuf,iwork_wf
          end if
       end do
999    continue
       write(*,*) "iwork_wf=",iwork_wf
    end if
    call mpi_bcast(iwork_wf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call write_border( 0, " read_wf(end)" )
  END SUBROUTINE read_wf


  SUBROUTINE init_wf( SYStype_in )
    implicit none
    integer,optional,intent(IN) :: SYStype_in
    integer :: SYStype

    call write_border( 0, " init_wf(start)" )

    SYStype=0
    if ( present(SYStype_in) ) SYStype=SYStype_in

    ML_WF   = sum( ir_grid )
    ML_0_WF = id_grid(myrank_g) + 1
    ML_1_WF = id_grid(myrank_g) + ir_grid(myrank_g)

    MB_WF   = sum( ir_band )
    MB_0_WF = id_band(myrank_b) + 1
    MB_1_WF = id_band(myrank_b) + ir_band(myrank_b)

    MK_WF   = sum( ir_bzsm )
    MK_0_WF = id_bzsm(myrank_k) + 1
    MK_1_WF = id_bzsm(myrank_k) + ir_bzsm(myrank_k)

    MS_WF   = sum( ir_spin )
    MS_0_WF = id_spin(myrank_s) + 1
    MS_1_WF = id_spin(myrank_s) + ir_spin(myrank_s)

    if ( allocated(occ) ) deallocate(occ)
    if ( allocated(res) ) deallocate(res)
    if ( allocated(esp) ) deallocate(esp)
    if ( allocated(unk) ) deallocate(unk)

    allocate( unk(ML_0_WF:ML_1_WF,MB_WF,MK_0_WF:MK_1_WF,MS_0_WF:MS_1_WF) )
    unk=zero
    allocate( esp(MB_WF,MK_WF,MS_WF) )
    esp=0.0d0
    allocate( res(MB_WF,MK_WF,MS_WF) )
    res=0.0d0
    allocate( occ(MB_WF,MK_WF,MS_WF) )
    occ=0.0d0

    call random_initial_wf
!    call fft_initial_wf_sub( ML_WF,MB_WF,MK_WF,MS_WF,ML_0_WF,ML_1_WF &
!        ,MB_0_WF,MB_1_WF,MK_0_WF,MK_1_WF,MS_0_WF,MS_1_WF,unk )
!    call random_initial_wf_sub( ML_WF,MB_WF,MK_WF,MS_WF,ML_0_WF,ML_1_WF &
!         ,MB_0_WF,MB_1_WF,MK_0_WF,MK_1_WF,MS_0_WF,MS_1_WF,unk,SYStype )

    if ( iwork_wf == 1 ) call allocate_work_wf( iwork_wf )

    call write_border( 0, " init_wf(end)" )

  END SUBROUTINE init_wf


  SUBROUTINE random_initial_wf
    implicit none
    integer :: s,k,n,i
    integer,allocatable :: ir(:)
    real(8) :: u(2)

    call random_seed( size=n )
    allocate( ir(n) )
    ir(:)=MB_0_WF+ML_0_WF
    call random_seed( put=ir )
    deallocate( ir )
 
    do s=MS_0_WF,MS_1_WF
       do k=MK_0_WF,MK_1_WF
          do n=MB_0_WF,MB_1_WF
             do i=ML_0_WF,ML_1_WF
                call random_number(u)
                unk(i,n,k,s)=dcmplx(u(1),u(2))
             end do
          end do
       end do
    end do

  END SUBROUTINE random_initial_wf


  SUBROUTINE test_on_wf(dV,disp_switch)
    implicit none
    real(8),intent(IN) :: dV          ! volume element
    logical,intent(IN) :: disp_switch ! diplay switch
    integer :: ierr,s,k,m,n,mm
#ifdef _DRSDFT_
    real(8),allocatable :: uu(:,:)
#else
    complex(8),allocatable :: uu(:,:)
#endif

    allocate( uu(MB_WF,MB_WF) ) ; uu=zero

    mm = size(uu)

    do s=MS_0_WF,MS_1_WF
    do k=MK_0_WF,MK_1_WF
       uu(:,:)=zero
       do n=1,MB_WF
       do m=1,n
#ifdef _DRSDFT_
          uu(m,n)=sum( unk(:,m,k,s)*unk(:,n,k,s) )*dV
#else
          uu(m,n)=sum( conjg(unk(:,m,k,s))*unk(:,n,k,s) )*dV
#endif
       end do ! m
       end do ! n
       call mpi_allreduce(MPI_IN_PLACE,uu,mm,TYPE_MAIN &
            ,MPI_SUM,comm_grid,ierr)
       do n=1,MB_WF
       do m=1,n
          if ( disp_switch ) then
             write(320,'(1x,i2,i5,2i7,2g25.16)') s,k,m,n,uu(m,n)
          end if
       end do ! m
       end do ! n
    end do ! k
    end do ! s

    deallocate( uu )

  END SUBROUTINE test_on_wf


  SUBROUTINE gather_wf
    implicit none
    integer :: k,s,mm,ierr
    call write_border( 1, " gather_wf(start)" )
    mm=ML_1_WF-ML_0_WF+1
    ir_band(:)=ir_band(:)*mm
    id_band(:)=id_band(:)*mm
    do s=MS_0_WF,MS_1_WF
    do k=MK_0_WF,MK_1_WF
       call mpi_allgatherv( unk(ML_0_WF,MB_0_WF,k,s),ir_band(myrank_b) &
            ,TYPE_MAIN,unk(ML_0_WF,1,k,s),ir_band,id_band &
            ,TYPE_MAIN,comm_band,ierr )
       if ( allocated(hunk) ) then
          call mpi_allgatherv( hunk(ML_0_WF,MB_0_WF,k,s),ir_band(myrank_b) &
               ,TYPE_MAIN,hunk(ML_0_WF,1,k,s),ir_band,id_band &
               ,TYPE_MAIN,comm_band,ierr )
       end if
    end do
    end do
    ir_band(:)=ir_band(:)/mm
    id_band(:)=id_band(:)/mm
    call write_border( 1, " gather_wf(end)" )
  END SUBROUTINE gather_wf

  SUBROUTINE write_wf(rankIN)
    implicit none
    integer,optional :: rankIN
    integer :: s,k,n,i
    integer :: rank
    if (present(rankIN)) then
      rank=rankIN+myrank
    else
      rank=myrank
    endif
    write(300+rank,*) 'myrank= ',rank
    write(300+rank,'(A18,2I5)') 'MS_0_WF, MS_1_WF= ',MS_0_WF,MS_1_WF
    write(300+rank,'(A18,2I5)') 'MK_0_WF, MK_1_WF= ',MK_0_WF,MK_1_WF
    write(300+rank,'(A18,2I5)') 'MB_0_WF, MB_1_WF= ',MB_0_WF,MB_1_WF
    write(300+rank,'(A18,2I5)') 'ML_0_WF, ML_1_WF= ',ML_0_WF,ML_1_WF
    do s=MS_0_WF,MS_1_WF
       do k=MK_0_WF,MK_1_WF
          do n=MB_0_WF,MB_1_WF
             do i=ML_0_WF,ML_1_WF
                write(300+rank,'(4I6,2g20.7)') s,k,n,i,unk(i,n,k,s)
             end do
          end do
       end do
    end do
    return

  END SUBROUTINE write_wf

  SUBROUTINE gather_b_wf( k, s )
    implicit none
    integer,intent(IN) :: k,s
    integer :: mm,ierr
    mm=ML_1_WF-ML_0_WF+1
    ir_band(:)=ir_band(:)*mm
    id_band(:)=id_band(:)*mm
    call mpi_allgatherv( unk(ML_0_WF,MB_0_WF,k,s),ir_band(myrank_b) &
            ,TYPE_MAIN,unk(ML_0_WF,1,k,s),ir_band,id_band &
            ,TYPE_MAIN,comm_band,ierr )
    ir_band(:)=ir_band(:)/mm
    id_band(:)=id_band(:)/mm
  END SUBROUTINE gather_b_wf


  SUBROUTINE allocate_work_wf( iflag )
    implicit none
    integer,intent(IN) :: iflag

    call write_border( 0, " allocate_work_wf(start)" )

    iflag_hunk=iflag
    if ( iwork_wf == 0 ) iflag_hunk=0

    if ( myrank == 0 ) then
       write(*,*) "iflag,iwork_wf,iflag_hunk=",iflag,iwork_wf,iflag_hunk
    end if

    if ( iflag == 1 ) then

       if ( allocated(hunk) ) deallocate(hunk)

       allocate( hunk(ML_0_WF:ML_1_WF,MB_WF,MK_0_WF:MK_1_WF,MS_0_WF:MS_1_WF) )

    else if ( iflag == 2 ) then

       if ( allocated(hunk) ) deallocate(hunk)

       allocate( hunk(ML_0_WF:ML_1_WF,MB_WF,MK_WF,MS_0_WF:MS_1_WF) )

    end if

    hunk(:,:,:,:)=zero

    if ( myrank == 0 ) then
       if ( TYPE_MAIN == MPI_COMPLEX16 ) then
          write(*,*) "size(hunk)(MB)=",size(hunk)*16.d0/1024.d0**2
       else if ( TYPE_MAIN == MPI_REAL8 ) then
          write(*,*) "size(hunk)(MB)=",size(hunk)*8.d0/1024.d0**2
       end if
    end if

    call write_border( 0, " allocate_work_wf(end)" )

  END SUBROUTINE allocate_work_wf


  SUBROUTINE deallocate_work_wf
    implicit none
    if ( allocated(hunk) ) deallocate(hunk)
  END SUBROUTINE deallocate_work_wf


  SUBROUTINE write_esp_wf( full_info )
    implicit none
    logical,optional,intent(IN) :: full_info
    integer :: k,n,s,i,n1,n2,nn,fi(2)
    real(8) :: f(6)
    character(57) :: header_string, format_string
    call check_disp_length( i, 0 ) ; if ( i < 1 ) return
    write(header_string,'(a4,a6,a20,2a13,1x)') &
         "k","n","esp(n,k,s)","esp_err  ","occ(n,k,s)  "
    call write_string( "" )
    call write_string( header_string )
    nn=sum(occ)
    n1=max( 1, nn/2-5 )
    n2=min( nn/2+5, size(esp,1) )
    if ( present(full_info) ) then
       if ( full_info ) then
          n1=1
          n2=size(esp,1)
       end if
    end if
    format_string='(i4,i6,2(f20.15,2g13.5,1x))'
    do k=1,size(esp,2)
    do n=n1,n2
       i=0
       do s=1,size(esp,3)
          i=i+1 ; f(i)=esp(n,k,s)
          i=i+1 ; f(i)=esp(n,k,s)-esp0(n,k,s)
          i=i+1 ; f(i)=occ(n,k,s)
       end do
       fi(1:2)=(/ k, n /)
       call write_int_and_real( format_string, 2, fi, i, f ) 
    end do
    end do
  END SUBROUTINE write_esp_wf


  SUBROUTINE write_info_esp_wf( control )
    implicit none
    integer,intent(IN) :: control
    integer :: s,k,n
    integer,parameter :: u=99
    call write_border( 1, " write_info_esp_wf(start)" )
    if ( control > 0 ) then
       if ( control == 2 ) rewind u
       write(u,*) "Eigenvalues"
       write(u,'(a4,a6,a20,2a13,1x)') &
            "k","n","esp(n,k,s)","esp_err  ","occ(n,k,s)  "
       do k=1,MK_WF
       do n=1,MB_WF
          write(u,'(i4,i6,2(f20.15,2g13.5,1x))') k,n &
               ,(esp(n,k,s),esp(n,k,s)-esp0(n,k,s),occ(n,k,s),s=1,MS_WF)
       end do
       end do
    end if
    call write_border( 1, " write_info_esp_wf(end)" )
  END SUBROUTINE write_info_esp_wf


  SUBROUTINE allocate_b_wf( b, wf )
    implicit none
    type(wfrange),intent(INOUT) :: b
#ifdef _DRSDFT_
    real(8),allocatable,intent(INOUT) :: wf(:,:,:,:)
#else
    complex(8),allocatable,intent(INOUT) :: wf(:,:,:,:)
#endif
    allocate( wf(b%ML0:b%ML1,b%MB0:b%MB1,b%MK0:b%MK1,b%MS0:b%MS1) )
    wf=(0.0d0,0.0d0)
  END SUBROUTINE allocate_b_wf

  SUBROUTINE allocate_b_occ( b, occup )
    implicit none
    type(wfrange),intent(INOUT) :: b
    real(8),allocatable,intent(INOUT) :: occup(:,:,:)
    allocate( occup(b%MB0:b%MB1,b%MK0:b%MK1,b%MS0:b%MS1) )
    occup=0.0d0
  END SUBROUTINE allocate_b_occ


END MODULE wf_module
