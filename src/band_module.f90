MODULE band_module

  use watch_module
  use parallel_module
  use rgrid_module
  use array_bound_module
  use bz_module
  use aa_module, only: aa
  use bb_module, only: bb
  use kinetic_module, only: init_kinetic
  use momentum_module
  use wf_module
  use cg_module
  use gram_schmidt_t_module
  use subspace_diag_module
  use esp_gather_module
  use electron_module
  use scalapack_module
  use band_variables, only: nfki,nbk,ak,nskip_band &
       ,esp_conv_tol, mb_band, mb2_band, maxiter_band, read_band &
       ,unit_band_eigv,unit_band_dedk,unit_band_ufld
  use sweep_module, only: calc_sweep, init_sweep
  use pseudopot_module, only: pselect
  use ps_nloc2_module, only: prep_uvk_ps_nloc2, prep_rvk_ps_nloc2
  use ps_nloc3_module, only: prep_ps_nloc3, init_ps_nloc3
  use io_module, only: Init_IO
  use xc_hybrid_module, only: iflag_hybrid, prep_kq_xc_hybrid
  use fock_fft_module, only: init_fock_fft
  use band_unfold_module
  use hsort_module
  use fermi_module, only: efermi
  
  implicit none

  PRIVATE
  PUBLIC :: band
  PUBLIC :: read_band

  integer :: unit = 1

  character(64),parameter :: version="version2.0"

CONTAINS


  SUBROUTINE band( MBV_in, disp_switch )

    implicit none
    integer,intent(IN) :: MBV_in
    logical,intent(IN) :: disp_switch
    integer :: MBV,nktrj,i,j,k,s,n,ibz,ierr,iktrj,iter,Diter_band
    integer :: iktrj_0,iktrj_1,iktrj_2,iktrj_00,iktrj_tmp,ireq
    integer,allocatable :: ir_k(:),id_k(:)
    real(8) :: dak(3),sum0,sum1,max_err,max_err0
    real(8),allocatable :: ktrj(:,:),pxyz(:,:,:,:)
    real(8),allocatable :: kbb_tmp(:,:),esp_tmp(:,:,:),esp0_tmp(:,:,:)
    logical :: disp_switch_parallel_bak,flag_end
    character(32) :: loop_info
    character(5) :: cc

#ifdef _DRSDFT_
    if ( disp_switch ) then
       write(*,*) "band calc is not available in REAL8 code"
       write(*,*) "Recompile"
    end if
    call stop_program( "" )
#endif

    call write_border( 0, "" )
    call write_border( 0, " BAND Calc. START -----------" )

    disp_switch_parallel_bak = disp_switch_parallel

    Diter_band = maxiter_band

    flag_end = .false.

    MBV=MBV_in
    if ( MBV < 1 .or. mb_band < MBV ) MBV=1
    if ( disp_switch_parallel ) write(*,*) "MBV=",MBV

! ---

    nktrj = sum( nfki(1:nbk) )
    if ( nktrj > 1 ) nktrj=nktrj+1

    allocate( ktrj(6,nktrj) ) ; ktrj=0.0d0

    if ( iswitch_banduf ) then

       call init_band_unfold( nktrj, ktrj, unit_band_ufld, disp_switch )

    else

       k=0
       do i=1,nbk
          dak(1:3) = ( ak(1:3,i+1) - ak(1:3,i) )/dble( nfki(i) )
          do j=1,nfki(i)
             k=k+1
             ktrj(1:3,k) = ak(1:3,i) + (j-1)*dak(1:3)
             ktrj(4:6,k) = matmul( bb(1:3,1:3),dak(1:3) )
          end do
       end do
       if ( nktrj > 1 ) then
          ktrj(1:3,k+1) = ak(1:3,nbk+1)
          ktrj(4:6,k+1) = 0.d0
       end if

    end if

    if ( disp_switch ) then
       write(*,*) "nktrj=",nktrj
       do k=1,nktrj
          write(*,'(1x,i4,2x,3f15.10,2x,3f15.10)') k,ktrj(1:6,k)
       end do
    end if

! ---

    allocate( ir_k(0:np_bzsm-1),id_k(0:np_bzsm-1) )
    id_k(:)=0
    ir_k(:)=0
    do k=1,nktrj
       i=mod(k-1,np_bzsm)
       ir_k(i)=ir_k(i)+1
    end do
    i=maxval(ir_k)
    ir_k(:)=i
    do i=0,np_bzsm-1
       id_k(i) = sum(ir_k(0:i)) - ir_k(i)
    end do
    if ( myrank == 0 ) then
       write(*,'(1x,3a6)') "rank","id_k","ir_k"
       do i=0,np_bzsm-1
          write(*,'(1x,3i6)') i,id_k(i),ir_k(i)
       end do
       write(*,*) "sum(ir_k),nktrj=",sum(ir_k),nktrj
    end if

! ---

    if ( iflag_hunk == 0 ) call deallocate_work_wf

    if ( MB < mb_band ) then
       call modify_mb
       call modify_arraysize
    end if

! ---

    allocate( esp0_tmp(MB,0:np_bzsm-1,MSP)   ) ; esp0_tmp=0.d0
    allocate( esp_tmp(MB,0:np_bzsm-1,MSP)    ) ; esp_tmp=0.d0
    allocate( kbb_tmp(3,0:np_bzsm-1)         ) ; kbb_tmp=0.d0
    allocate( pxyz(3,MB,0:np_bzsm-1,MSP)     ) ; pxyz=0.d0

! ---

    if ( myrank == 0 ) then
       open(unit_band_eigv,file="band_eigv")
       open(unit_band_dedk,file="band_dedk")
       write(unit_band_eigv,*) version
       write(unit_band_eigv,'(1x,"Fermi_energy(hartree):",f20.15)') efermi
       write(unit_band_eigv,'(1x,"Reciprocal_Lattice_Vectors:")')
       write(unit_band_eigv,'(1x,3f20.15)') bb(1:3,1)
       write(unit_band_eigv,'(1x,3f20.15)') bb(1:3,2)
       write(unit_band_eigv,'(1x,3f20.15)') bb(1:3,3)
    end if

! ---

    iktrj_0 = id_k(myrank_k)+1
    iktrj_1 = id_k(myrank_k)+ir_k(myrank_k)
    iktrj_2 = id_k(myrank_k)+maxval(ir_k)

    call init_sweep( 2, mb2_band, esp_conv_tol )

    call Init_IO( "band" )

    MBZ_1 = MBZ_0

    loop_iktrj : do iktrj = iktrj_0, iktrj_2

       iktrj_00 = id_k(0) + 1 + iktrj - iktrj_0

       if ( iktrj_00 <= nskip_band ) then
          if ( DISP_SWITCH_PARALLEL ) then
             do i=0,np_bzsm-1
                write(*,*) "Band ",iktrj_00-id_k(0)+id_k(i)," is skipped"
             end do
          end if
          cycle
       end if

       if ( iktrj <= nktrj ) then
          kbb(1:3,MBZ_0) = ktrj(1:3,iktrj)
       else
          kbb(1:3,MBZ_0) = 0.0d0
       end if
       call mpi_allgather &
            (kbb(1,MBZ_0),3,mpi_real8,kbb_tmp,3,mpi_real8,comm_bzsm,ierr)
       if ( myrank == 0 ) then
          write(*,*) "kbb_tmp"
          do ibz=0,np_bzsm-1
             iktrj_tmp = id_k(ibz)+iktrj-iktrj_0+1
             write(*,'(1x,i8,3f10.5)') iktrj_tmp,kbb_tmp(1:3,ibz)
          end do
       endif

       call init_kinetic(aa,bb,Nbzsm,kbb,DISP_SWITCH=disp_switch)

       select case( pselect )
       case( 2 )
          call prep_uvk_ps_nloc2(MBZ_0,MBZ_0,kbb(1,MBZ_0))
       case( 3 )
          call init_ps_nloc3
          call prep_ps_nloc3
       end select

       if ( iflag_hybrid > 0 ) then
          if ( disp_switch ) write(*,*) "iflag_hybrid=",iflag_hybrid
          call prep_kq_xc_hybrid(Nbzsm,MBZ_0,MBZ_0,kbb,bb,1)
          call init_fock_fft
       end if

! --- sweep ---

       write(loop_info,'("( iktrj=",i4," )")') iktrj
       call calc_sweep( disp_switch, ierr, Diter_band, loop_info )

! ---

       if ( ierr == -1 ) then
          if ( myrank == 0 ) write(*,*) "etime limit !!!"
          exit loop_iktrj
       end if
       if ( ierr == -2 ) then
          write(*,*) "band is not converged"
          return
       end if

! --- band unfolding ---

       call band_unfold( iktrj, disp_switch )

! ---

       select case( pselect )
       case( 2 )
          call prep_rvk_ps_nloc2(MBZ_0,MBZ_0,kbb(1,MBZ_0))
       case( 3 )
!          call prep_rvk_ps_nloc3(MBZ_0,MBZ_0,kbb(1,MBZ_0))
       end select

       pxyz(:,:,:,:)=0.d0
#ifndef _DRSDFT_
       do n=1,mb2_band
       do s=MSP_0,MSP_1
          call calc_expectval_momentum &
               (MBZ_0,ML_0,ML_1,1,1,unk(ML_0,n,MBZ_0,s),pxyz(1,n,myrank_k,s))
       end do ! s
       end do ! n
#endif

       do s=MSP_0,MSP_1
          call mpi_allgather(pxyz(1,1,myrank_k,s),3*MB,MPI_REAL8 &
               ,pxyz(1,1,0,s),3*MB,MPI_REAL8,comm_bzsm,ierr)
       end do ! s
       call mpi_allgather(pxyz(1,1,0,MSP_0),3*MB*np_bzsm*(MSP_1-MSP_0+1),MPI_REAL8 &
            ,pxyz,3*MB*np_bzsm*(MSP_1-MSP_0+1),MPI_REAL8,comm_spin,ierr)

       do s=1,MSP
          call mpi_allgather(esp(1,MBZ_0,s),MB,mpi_real8,esp_tmp(1,0,s),MB,mpi_real8,comm_bzsm,ierr)
       end do
       do s=1,MSP
          esp0_tmp(:,myrank_k,s)=esp(:,MBZ_0,s)
          call mpi_allgather(esp0_tmp(1,myrank_k,s),MB,mpi_real8,esp0_tmp(1,0,s),MB,mpi_real8,comm_bzsm,ierr)
       end do

       do ibz=0,np_bzsm-1
          iktrj_tmp = iktrj_00 + id_k(ibz)
          if ( iktrj_tmp > nktrj ) exit
          if ( myrank == 0 ) then
             write(unit_band_eigv,'(1x,2i6,3f20.12,i8)') &
                  iktrj_tmp,mb2_band,kbb_tmp(1:3,ibz),MBV
             write(unit_band_dedk,'(1x,2i6,3f20.12)') &
                  iktrj_tmp,mb2_band,kbb_tmp(1:3,ibz)
          end if
          do n=1,mb2_band
             if ( disp_switch ) then
                write(*,'(1x,i5,2x,2(1x,g22.12,1x,g15.5))') n &
           ,( esp_tmp(n,ibz,s),abs(esp_tmp(n,ibz,s)-esp0_tmp(n,ibz,s)),s=1,MSP )
             end if
             if ( myrank == 0 ) then
                write(unit_band_eigv,'(1x,i5,2x,2(1x,g22.12,1x,g15.5))') n &
                     ,( esp_tmp(n,ibz,s),abs(esp_tmp(n,ibz,s)-esp0_tmp(n,ibz,s)),s=1,MSP )
                write(unit_band_dedk,'(1x,i5,2x,2(g22.12,1x,3g15.5))') n &
                     ,( esp_tmp(n,ibz,s),pxyz(1:3,n,ibz,s),s=1,MSP )
             end if
          end do ! n
       end do ! ibz

    end do loop_iktrj

    if ( myrank == 0 ) then
       close(unit_band_dedk)
       close(unit_band_eigv)
    end if

    deallocate( pxyz )
    deallocate( kbb_tmp )
    deallocate( esp_tmp )
    deallocate( esp0_tmp )
    deallocate( id_k )
    deallocate( ir_k )
    deallocate( ktrj )

    disp_switch_parallel = disp_switch_parallel_bak

    call finalize_band_unfold

    call write_border( 0, " BAND Calc. END -----------" )
    call write_border( 0, "" )

  END SUBROUTINE band

  SUBROUTINE send_wf_band(mm,nn,f)
    implicit none
    integer,intent(IN) :: mm,nn
    complex(8),intent(IN) :: f(mm,nn)
    integer :: ierr,ireq,istatus(MPI_STATUS_SIZE)
    if ( myrank_k > 0 ) then
       call MPI_ISEND(f,mm*nn,MPI_COMPLEX16,myrank_k-1,1,comm_bzsm,ireq,ierr)
       call MPI_WAIT(ireq,istatus,ierr)
    end if
  END SUBROUTINE send_wf_band
  SUBROUTINE recv_wf_band(mm,nn,f)
    implicit none
    integer,intent(IN) :: mm,nn
    complex(8),intent(OUT) :: f(mm,nn)
    integer :: ierr,ireq,istatus(MPI_STATUS_SIZE)
    if ( myrank_k < np_bzsm-1 ) then
       call MPI_IRECV(f,mm*nn,MPI_COMPLEX16,myrank_k+1,1,comm_bzsm,ireq,ierr)
       call MPI_WAIT(ireq,istatus,ierr)
    end if
  END SUBROUTINE recv_wf_band


  SUBROUTINE modify_mb

    implicit none

    integer :: i,j

    if ( mb_band <= Nband ) return

    Nband = mb_band

    call init_scalapack( Nband )

    ir_band(:)=0
    do i=0,Nband-1
       j=mod(i,np_band)
       ir_band(j) = ir_band(j) + 1
    end do
    do i=0,np_band-1
       id_band(i) = sum( ir_band(0:i) ) - ir_band(i)
    end do
    MB   = Nband
    MB_0 = id_band(myrank_b) + 1
    MB_1 = id_band(myrank_b) + ir_band(myrank_b)

    call init_subspace_diag( Nband )

  END SUBROUTINE modify_mb


  SUBROUTINE modify_arraysize

    implicit none

    complex(8),allocatable :: utmp(:,:,:,:)
    integer :: ml_old,mb_old,mbz_old,msp_old

    ml_old  = size( unk,1 )
    mb_old  = size( unk,2 )
    mbz_old = size( unk,3 )
    msp_old = size( unk,4 )

    allocate( utmp(ml_old,mb_old,mbz_old,msp_old) )
    utmp=unk

    call init_wf

    unk(ML_0:ML_0+ml_old-1,1:mb_old,MBZ_0:MBZ_0+mbz_old-1,MSP_0:MSP_0+msp_old-1) &
         = utmp(1:ml_old,1:mb_old,1:mbz_old,1:msp_old)

    deallocate( utmp )

  END SUBROUTINE modify_arraysize


END MODULE band_module
