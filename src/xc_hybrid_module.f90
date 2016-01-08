MODULE xc_hybrid_module

  use io_tools_module
  use wf_module, only: wfrange
  use io_read_module

  implicit none

  PRIVATE
  PUBLIC :: init_xc_hybrid, control_xc_hybrid &
           ,get_flag_xc_hybrid &
           ,omega, R_hf, alpha_hf, q_fock, gamma_hf &
           ,iflag_hf, iflag_pbe0, iflag_hse, iflag_lcwpbe, iflag_hybrid &
           ,FOCK_0, FOCK_1, FKMB_0, FKMB_1, FKBZ_0, FKBZ_1 &
           ,VFunk, unk_hf, occ_hf, occ_factor, npart &
           ,n_kq_fock, i_kq_fock, kq_fock, prep_kq_xc_hybrid
  PUBLIC :: read_xc_hybrid

  integer :: npart
  real(8) :: R_hf

  integer :: iflag_hf     = 0
  integer :: iflag_pbe0   = 0
  integer :: iflag_hse    = 0
  integer :: iflag_lcwpbe = 0
  integer :: iflag_hybrid = 0

#ifdef _DRSDFT_
  real(8),allocatable :: VFunk(:,:,:,:)
  real(8),allocatable :: unk_hf(:,:,:,:)
  real(8),parameter :: byte = 8.0d0
#else
  complex(8),allocatable :: VFunk(:,:,:,:)
  complex(8),allocatable :: unk_hf(:,:,:,:)
  real(8),parameter :: byte = 16.0d0
#endif

  real(8),allocatable :: occ_hf(:,:,:)
  real(8) :: occ_factor

  real(8),allocatable :: kbb_hf(:,:)
  real(8),allocatable :: q_fock(:,:,:)

  integer :: gamma_hf
  integer :: FKBZ,FKBZ_0,FKBZ_1,FKMMBZ
  integer :: FKMB,FKMB_0,FKMB_1
  integer :: FOCK_0, FOCK_1

  integer :: n_kq_fock
  integer,allocatable :: i_kq_fock(:,:,:)
  real(8),allocatable :: kq_fock(:,:)

  logical :: flag_init = .false.

  real(8) :: omega = 0.11d0
  real(8) :: alpha_hf = 0.25d0
  integer :: IC = 0
  integer :: IO_ctrl = 0
  character(30) :: file_wf2 = "wf.dat1"
  character(8) :: XCtype

CONTAINS


  SUBROUTINE read_xc_hybrid
    implicit none
    real(8) :: tmp(2)
    call write_border( 0, " read_xc_hybrid(start)" )
    call IOTools_readStringKeyword( "XCTYPE", XCtype )
    tmp(1:2) = (/ omega, alpha_hf /)
    call IOTools_readReal8Keywords( "HF", tmp )
    omega    = tmp(1)
    alpha_hf = tmp(2)
    call IOTools_readIntegerKeyword( "IC", IC )
    call IOTools_readIntegerKeyword( "IOCTRL", IO_ctrl )
    call write_border( 0, " read_xc_hybrid(end)" )
  END SUBROUTINE read_xc_hybrid


  SUBROUTINE init_xc_hybrid( n1, n2, Ntot, Nspin, MB, MMBZ &
       , MBZ,MBZ_0,MBZ_1, MSP,MSP_0,MSP_1, MB_0,MB_1, kbb, bb, Vcell &
       , SYStype, np_fkmb, disp_switch )
    implicit none
    integer,intent(IN) :: n1, n2, Nspin, MB, MBZ,MMBZ,MBZ_0,MBZ_1, SYStype
    integer,intent(IN) :: MSP, MSP_0, MSP_1, MB_0, MB_1, np_fkmb
    real(8),intent(IN) :: Ntot, kbb(:,:), bb(3,3), Vcell
    logical,intent(IN) :: disp_switch
    integer :: ML0,i,s,k,q,init_num,ierr,m,t
    integer,allocatable :: ir(:),id(:)
    real(8) :: ctime0,ctime1,etime0,etime1,best_time,time
    real(8) :: ctime_hf0,ctime_hf1,etime_hf0,etime_hf1
    real(8),parameter :: eps=1.d-5
    real(8) :: mem(9),qtry(3),c,Pi,k_fock(3)
    real(8),allocatable :: qtmp(:,:)
    type(wfrange) :: b

    call write_border( 0, " init_xc_hybrid(start)" )

    if ( XCtype /= "HF"    .and. XCtype /= "HSE"    .and. &
         XCtype /= "HSE06" .and. XCtype /= "HSE_"   .and. &
         XCtype /= "PBE0"  .and. XCtype /= "LCwPBE" ) return

    if ( flag_init ) return

!
! --- set flags ---
!

    iflag_hf     = 0
    iflag_hse    = 0
    iflag_pbe0   = 0
    iflag_lcwpbe = 0
    iflag_hybrid = 0

    select case( XCtype )
    case( "HF" )
       iflag_hf     = 1
       alpha_hf     = 1.0d0
    case( "HSE","HSE06","HSE_" )
       iflag_hse    = 1
       alpha_hf     = 0.25d0
    case( "PBE0" )
       iflag_pbe0   = 1
       alpha_hf     = 0.25d0
    case( "LCwPBE" )
       iflag_lcwpbe = 1
       alpha_hf     = 1.0d0
    case default
       goto 99
    end select

    if ( disp_switch ) then
       write(*,*) "XCtype       =",XCtype
       write(*,*) "iflag_hf     =",iflag_hf
       write(*,*) "iflag_hse    =",iflag_hse
       write(*,*) "iflag_pbe0   =",iflag_pbe0
       write(*,*) "iflag_lcwpbe =",iflag_lcwpbe
       write(*,*) "iflag_hybrid =",iflag_hybrid
       write(*,*) "alpha_hf     =",alpha_hf
       write(*,*) "omega        =",omega
    end if

    if ( omega <= 0.0d0 ) then
       if ( iflag_hse == 1 .or. iflag_lcwpbe == 1 ) then
          write(*,*) "omega should be greater than zero for HSE or LCwPBE"
          stop "stop@init_xc_hybrid"
       end if
    end if

!
! --- Coefficient of occupation number in exact exchange ---
!

    if ( Ntot < 1.0d0+eps ) then
       occ_factor=0.5d0
    else
       occ_factor=0.25d0*dble(Nspin)
    end if

!
! --- temp ---
!

    FKBZ   = MBZ
    FKBZ_0 = 1
    FKBZ_1 = MBZ
    FKMMBZ = MMBZ
    FKMB   = MB
    FKMB_0 = 1
    FKMB_1 = MB/np_fkmb
    if ( FKMB_1*np_fkmb < MB ) FKMB_1=FKMB_1+1
    FOCK_0 = 1
    FOCK_1 = 2

!
! --- allocate ---
!

    if ( IC == 0 ) then

       mem(1)=byte*(n2-n1+1)*(FKMB_1-FKMB_0+1) &
            *(FKBZ_1-FKBZ_0+1)*(MSP_1-MSP_0+1)

       if ( disp_switch ) then
          write(*,*) "size(unk_hf)(MB)=",mem(1)/1024.d0**2
       end if

       allocate( unk_hf(n1:n2,FKMB_0:FKMB_1,FKBZ_0:FKBZ_1,MSP_0:MSP_1) )
       unk_hf=(0.0d0,0.0d0)
       allocate( occ_hf(FKMB_0:FKMB_1,FKBZ_0:FKBZ_1,MSP) )
       occ_hf=0.0d0

    else if ( IC > 0 ) then

!       call read_xc_hybrid_io( file_wf2, SYStype, IO_ctrl, disp_switch &
!            ,n1,n2, MSP_0,MSP_1, unk_hf, occ_hf, kbb_hf &
!            ,FKMB,FKMB_0,FKMB_1,FKBZ,FKBZ_0,FKBZ_1,FKMMBZ )

       b%ML0 = n1
       b%ML1 = n2
       b%MB  = 0
       b%MB0 = 0
       b%MB1 = 0
       b%MK  = 0
       b%MK0 = 0
       b%MK1 = 0
       b%MMK = 0
       b%MS  = MSP
       b%MS0 = MSP_0
       b%MS1 = MSP_1

       call simple_wf_io_read( file_wf2, SYStype, IO_ctrl, disp_switch &
            , b, unk_hf, occ_hf, kbb_hf )

       FKMB   = b%MB
       FKMB_0 = b%MB0
       FKMB_1 = b%MB1
       FKBZ   = b%MK
       FKBZ_0 = b%MK0
       FKBZ_1 = b%MK1
       FKMMBZ = b%MMK

    end if

! ---

    if ( .not.allocated(kbb_hf) ) then

       allocate( kbb_hf(3,FKBZ_0:FKBZ_1) ) ; kbb_hf=0.0d0
       kbb_hf(1:3,FKBZ_0:FKBZ_1)=kbb(1:3,FKBZ_0:FKBZ_1)

    else

       do k=FKBZ_0,FKBZ_1
          if ( disp_switch ) write(*,'(1x,i4,3f20.15)') k,kbb_hf(1:3,k)
       end do

    end if

!
! --- Coordinates at reciprocal space in hybrid DFT calculation ---
!

    if ( SYStype == 0 ) then

       allocate( q_fock(3,FKBZ_0:FKBZ_1,2) ) ; q_fock=0.0d0

       do q=FKBZ_0,FKBZ_1
          q_fock(:,q,1) = bb(:,1)*kbb_hf(1,q) &
                        + bb(:,2)*kbb_hf(2,q) &
                        + bb(:,3)*kbb_hf(3,q)
       end do

       q_fock(:,:,2) = -q_fock(:,:,1)

       call prep_kq_xc_hybrid( MBZ, MBZ_0, MBZ_1, kbb, bb, 0 )

    end if

!
! --- Switch of hybrid DFT calculation on Gamma-point ---
!

    gamma_hf = 0

    if ( FKMMBZ == 1 .and. all(kbb(:,1)==0.0d0)  ) then

       gamma_hf = 1

       if ( disp_switch ) then      
          write(*,'(1x,"Hybrid DFT calculation on Gamma-point (gamma_hf=" &
                   ,i1,")")') gamma_hf
       end if

    end if

    if ( gamma_hf == 0 .and. np_fkmb > 1 ) then
       call stop_program( " np_fkmb>1 is not available with gamma_hf==0" )
    end if

!
! --- Truncation cutoff of 1/r for HF, PBE0, and LCwPBE functionals ---
!

    Pi = acos(-1.0d0)

    if ( SYStype == 0 ) then

       if ( iflag_hf /= 0 .or. iflag_pbe0 /= 0 .or. iflag_lcwpbe /= 0 ) then

          R_hf = ( 0.75d0*abs(Vcell)*FKMMBZ/Pi )**(1.0d0/3.0d0)

          if ( disp_switch ) then
             write(*,*) "Truncation cutoff of 1/r (Ang.) =",R_hf*0.529177d0
          end if

       end if

    end if

!
! Preparation of Fock potential in HSE and LCwPBE hybrid functionals for RSMOL
!

    if ( SYStype == 1 ) then

!       if ( iflag_hse /= 0 .or. iflag_lcwpbe /= 0 ) call prep_hse_fock
       if ( iflag_hse /= 0 .or. iflag_lcwpbe /= 0 ) stop "stop@init_xc_hybrid"

    end if

!
! --- Divided MPI_Allgatherv ---
!

    npart = 30

!    call gather_wf
!    ML0 = n2 - n1 + 1
!    init_num=min(ML0,nint(size(unk)/1.024d3))
!    best_time=1.d15
!    do i=init_num,init_num-4,-1
!       call bwatch(ctime_hf0,etime_hf0)
!       do s=MSP_0,MSP_1
!          do k=MBZ_0,MBZ_1
!             call rsdft_allgatherv(ML0,MB_0,MB_1,MB,unk(n1,1,k,s),np_band,comm_band,myrank_b,i)
!          end do
!          call rsdft_allgatherv(ML0*MB,MBZ_0,MBZ_1,MBZ,unk(n1,1,1,s),np_bzsm,comm_bzsm,myrank_k,i)
!       end do
!       call bwatch(ctime_hf1,etime_hf1)
!       time=ctime_hf1-ctime_hf0
!       if ( time < best_time ) then
!          best_time=time
!          npart=i
!       end if
!    end do
!    call mpi_bcast(npart,1,mpi_integer,0,mpi_comm_world,ierr)

    if ( disp_switch ) then
       write(*,*) "Division number of mpi_allgatherv =",npart
!       write(*,*) "Time of divided MPI_Allgatherv (s) =",time
!       write(*,*) "Optimal division number of mpi_allgatherv =",npart
!       write(*,*) "Best time of divided MPI_Allgatherv (s) =",best_time
    end if

! ---

99  continue

    flag_init = .true.

    call write_border( 0, " init_xc_hybrid(end)" )

    return 
 
  END SUBROUTINE init_xc_hybrid


  SUBROUTINE control_xc_hybrid( ictrl )
    implicit none
    integer,intent(IN) :: ictrl
    if ( iflag_hybrid == 3 ) return
    call write_border( 1, " control_xc_hybrid(start)" )
    if ( iflag_hf   == 0 .and. iflag_hse    == 0 .and. &
         iflag_pbe0 == 0 .and. iflag_lcwpbe == 0 ) then
       iflag_hybrid = 0
    else
       iflag_hybrid = ictrl
    end if
    call write_border( 1, " control_xc_hybrid(end)" )
  END SUBROUTINE control_xc_hybrid


  SUBROUTINE get_flag_xc_hybrid( iflag )
    implicit none
    integer,intent(OUT) :: iflag
    iflag = iflag_hybrid
  END SUBROUTINE get_flag_xc_hybrid


  SUBROUTINE prep_kq_xc_hybrid( MBZ, MBZ_0, MBZ_1, kbb, bb, iswitch )
    implicit none
    integer,intent(IN) :: MBZ,MBZ_0,MBZ_1
    real(8),intent(IN) :: kbb(3,MBZ), bb(3,3)
    integer,intent(IN) :: iswitch
    integer :: m,i,k,q,t,s
    real(8) :: qtry(3),k_fock(3),c
    real(8),allocatable :: qtmp(:,:)

    if ( allocated(i_kq_fock) ) deallocate( i_kq_fock )
    if ( allocated(kq_fock) ) deallocate( kq_fock )

    if ( iswitch == 0 ) then

       allocate( i_kq_fock(FKBZ_0:FKBZ_1,FKBZ_0:FKBZ_1,2) )
       i_kq_fock=0

       m=2*(FKBZ_1-FKBZ_0+1)**2
       allocate( qtmp(3,m) ) ; qtmp=0.0d0

       i=0
       do k=FKBZ_0,FKBZ_1
       do q=FKBZ_0,FKBZ_1
          do t=1,2
             qtry(1:3)=q_fock(:,k,1)-q_fock(:,q,t)
             if ( i == 0 ) then
                i=i+1
                qtmp(:,i)=qtry(:)
                i_kq_fock(k,q,t)=i
             else
                do s=1,i
                   c=sum((qtry(:)-qtmp(:,s))**2)
                   if ( c < 1.d-12 ) then
                      i_kq_fock(k,q,t)=s
                      exit
                   end if
                end do
                if ( s > i ) then
                   i=i+1
                   qtmp(:,i)=qtry(:)
                   i_kq_fock(k,q,t)=i
                end if
             end if
          end do ! t
       end do ! q
       end do ! k

    else if ( iswitch == 1 ) then

       allocate( i_kq_fock(MBZ_0:MBZ_1,FKBZ_0:FKBZ_1,2) )
       i_kq_fock=0

       m=2*(MBZ_1-MBZ_0+1)*(FKBZ_1-FKBZ_0+1)
       allocate( qtmp(3,m) ) ; qtmp=0.0d0

       i=0
       do k=MBZ_0 ,MBZ_1
          k_fock(:)=bb(:,1)*kbb(1,k)+bb(:,2)*kbb(2,k)+bb(:,3)*kbb(3,k)
       do q=FKBZ_0,FKBZ_1
          do t=1,2
             qtry(1:3)=k_fock(:)-q_fock(:,q,t)
             if ( i == 0 ) then
                i=i+1
                qtmp(:,i)=qtry(:)
                i_kq_fock(k,q,t)=i
             else
                do s=1,i
                   c=sum((qtry(:)-qtmp(:,s))**2)
                   if ( c < 1.d-12 ) then
                      i_kq_fock(k,q,t)=s
                      exit
                   end if
                end do
                if ( s > i ) then
                   i=i+1
                   qtmp(:,i)=qtry(:)
                   i_kq_fock(k,q,t)=i
                end if
             end if
          end do ! t
       end do ! q
       end do ! k

       iflag_hybrid = 3
!       if ( disp_switch ) write(*,*) "iflag_hybrid=",iflag_hybrid

    end if

    n_kq_fock = i
    allocate( kq_fock(3,n_kq_fock) ) ; kq_fock=0.0d0
    kq_fock(:,1:n_kq_fock) = qtmp(:,1:n_kq_fock)

!    if ( disp_switch ) then
!       write(*,*) "n_kq_fock",n_kq_fock
!       write(*,'(1x,4x,5x,a)') "kq_fock(1:3)"
!       do i=1,n_kq_fock
!          write(*,'(1x,i4,3f12.5)') i,kq_fock(:,i)
!       end do
!    end if

    deallocate( qtmp )

  END SUBROUTINE prep_kq_xc_hybrid


END MODULE xc_hybrid_module
