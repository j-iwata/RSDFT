MODULE test_rtsol_module

  use wf_module
  use hamiltonian_module
  use watch_module

  implicit none

  PRIVATE
  PUBLIC :: test_rtsol

  real(8) :: dt    = 0.001d0
  integer :: nt    = 10
  integer :: nalg  = 4
  integer :: MB_RT = 16
  complex(8),parameter :: zi=(0.0d0,1.0d0)
  complex(8),parameter :: z0=(0.0d0,0.0d0)

CONTAINS


  SUBROUTINE test_rtsol

#ifdef _DRSDFT_
    write(*,*) "test_rtsol is available only for COMPLEX16 version !!!"
    return
#else

    implicit none
    logical :: disp_sw
    integer :: i,it,itaylor,n,k,s
    real(8) :: c,ct(2),et(2)
    complex(8),allocatable :: zcoef(:)
    complex(8),allocatable :: tpsi(:),hpsi(:)

    call write_border(60," test_rtsol(start)")
    call check_disp_switch( disp_sw, 0 )

    if ( disp_sw ) then
       write(*,*) "dt   =",dt
       write(*,*) "nt   =",nt
       write(*,*) "nalg =",nalg
       write(*,*) "MB_RT=",MB_RT
       write(*,*) "MK_WF=",MK_WF
       write(*,*) "ML_WF=",ML_WF
       write(*,*) "MS_WF=",MS_WF
    end if

    allocate( zcoef(nalg) ) ; zcoef=z0
    allocate( tpsi(ML_0_WF:ML_1_WF) ) ; tpsi=z0
    allocate( hpsi(ML_0_WF:ML_1_WF) ) ; hpsi=z0

    do itaylor=1,nalg
       c=1.0d0
       do i=1,itaylor
          c=c*i
       end do
       zcoef(itaylor) = (-zi*dt)**itaylor/c
    end do

    time_hmlt(:,:)=0.0d0
    call watch(ct(1),et(1))

    do it=1,nt

       do s=MS_0_WF,MS_1_WF
       do k=MK_0_WF,MK_1_WF
       do n=1,MB_RT
          do itaylor=1,nalg
             call hamiltonian(k,s,tpsi,hpsi,ML_0_WF,ML_1_WF,1,1)
!$OMP parallel do
             do i=ML_0_WF,ML_1_WF
                unk(i,n,k,s) = unk(i,n,k,s) + zcoef(itaylor)*hpsi(i)
                tpsi(i) = hpsi(i)
             end do
!$OMP end parallel do
          end do
       end do ! n
       end do ! k
       end do ! s

    end do ! it

    call watch(ct(2),et(2))
    if ( disp_sw ) then
       call write_watchb( time_hmlt, 4, time_hmlt_indx )
       write(*,*) "time(total)",ct(2)-ct(1),et(2)-et(1)
    end if

    deallocate( hpsi  )
    deallocate( tpsi  )
    deallocate( zcoef )

    call write_border(60," test_rtsol(end)")

#endif

  END SUBROUTINE test_rtsol


END MODULE test_rtsol_module
