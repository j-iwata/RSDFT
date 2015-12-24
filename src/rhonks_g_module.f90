MODULE rhonks_g_module

  use rgrid_module, only: dV
  use wf_module, only: unk
  use parallel_module, only: myrank
  use ps_nloc2_variables, only: Mlma,sbufnl,rbufnl,sendmap,recvmap,lma_nsend &
                               ,nzlma,nrlma_xyz,num_2_rank,uVk,JJP,MJJ
  use pseudopot_module, only: pselect
  use var_ps_member_g, only: uVunk, QRij, nzqr_pair, N_nzqr
  use var_para_ps_nloc_g, only: MJJ_Q,JJP_Q
  use para_rgrid_comm, only: do3StepComm

  implicit none

  include 'mpif.h'

  PRIVATE
  PUBLIC :: get_rhonks

CONTAINS

  SUBROUTINE get_rhonks( rhonks,nn1,nn2,n,k,s )
    ! IN:	nn1,nn2,n,k,s,
    !		unk(nn1:nn2,n,k,s),dV
    !		uVnk(MJJ(lma),lma,k),uVnk0(MJJ(lma),lma,k),
    !		Mlma,nzlma,MJJ(lma),JJP(MJJ(lma),lma),
    !		nrlma_xyz(1:6),num_2_rank(nrlma_xyz(i),i),
    !           lma_nsend(num_2_rank(m,i)),sbufnl(i2,irank),
    !           rbufnl(i2,jrank),TYPE_MAIN,
    !		sendmap(i1,irank),recvmap(i1,jrank),uVk(:,:,:)
    !		N_nzqr,nzqr_pair(N_nzqr,2),
    !		MJJ_Q(N_nzqr),JJP_Q(MJJ_Q(kk1),kk1),QRij(MJJ_Q(kk1),kk1),
    !		myrank
    ! OUT:	rhonks(:)

    implicit none
    integer,intent(IN) :: nn1,nn2
    integer,intent(IN) :: n,k,s
    real(8),intent(INOUT) :: rhonks(nn1:nn2)
#ifdef _DRSDFT_
    real(8) :: Qrhonks(nn1:nn2)
    real(8) :: p_uVunk1,p_uVunk2
    real(8) :: uVunk1,uVunk2
    real(8),parameter :: zero=0.d0
    real(8) :: uVunk2_z(Mlma)
    real(8) :: tmp
#else
    complex(8) :: Qrhonks(nn1:nn2)
    complex(8) :: p_uVunk1,p_uVunk2
    complex(8) :: uVunk1,uVunk2
    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8) :: uVunk2_z(Mlma)
    complex(8) :: tmp
#endif
    integer :: kk1,i,j
    integer :: ierr
    integer :: lma,irank,jrank,i1,i2,nreq
    integer :: ib1,ib2,nb,ib,lma1,lma2,m
    integer :: istatus(MPI_STATUS_SIZE,512),ireq(512)

    select case ( pselect )
    case ( 102 )

!----- term1 -----
       rhonks(nn1:nn2) = abs( unk(nn1:nn2,n,k,s) )**2
!===== term1 =====

!----- term2 -----

!----- get_uVunk -----

       ib1=n
       ib2=n
       nb = ib2 - ib1 + 1
       allocate( uVunk(nzlma,ib1:ib2)  ) ; uVunk(:,:) =zero
       do ib=ib1,ib2
          do lma=1,nzlma
             do j=1,MJJ(lma)
                i=JJP(j,lma)
#ifdef _DRSDFT_
                uVunk(lma,ib)=uVunk(lma,ib)+uVk(j,lma,k)*unk(i,ib,k,s)
#else
                uVunk(lma,ib)=uVunk(lma,ib)+conjg(uVk(j,lma,k))*unk(i,ib,k,s)
#endif
             end do
             uVunk(lma,ib) = dV*uVunk(lma,ib)
          end do	! lma
       end do	! ib


! 3WayComm
       call do3StepComm(nrlma_xyz,num_2_rank,sendmap,recvmap &
                       ,lma_nsend,sbufnl,rbufnl,nzlma,ib1,ib2,uVunk)

!===== term2 =====

!----- get Qrhonks -----

       Qrhonks=zero
       do ib=ib1,ib2
          do kk1=1,N_nzqr
             lma1=nzqr_pair(kk1,1)
             lma2=nzqr_pair(kk1,2)
             if ( lma1 < lma2 ) stop 'NZQR_PAIR is stange'

             if ( lma1 == lma2 ) then

                do j=1,MJJ_Q(kk1)
                   i=JJP_Q(j,kk1)
                   if ( i<nn1 .or. nn2<i ) then
                      write(*,'(A,6I6)') &
                           'STOP *****  myrank,i,nn1,nn2,kk1,j=' &
                           ,myrank,i,nn1,nn2,kk1,j
                      stop
                   end if
#ifdef _DRSDFT_
                   Qrhonks(i) = Qrhonks(i) &
                        + QRij(j,kk1)*uVunk(lma1,ib)*uVunk(lma2,ib)
#else
                   Qrhonks(i) = Qrhonks(i) &
                        + QRij(j,kk1)*conjg(uVunk(lma1,ib))*uVunk(lma2,ib)
#endif
                end do ! j

             else

                do j=1,MJJ_Q(kk1)
                   i=JJP_Q(j,kk1)
#ifdef _DRSDFT_
                   Qrhonks(i) = Qrhonks(i) &
                        + QRij(j,kk1)*uVunk(lma1,ib)*uVunk(lma2,ib)
                   Qrhonks(i) = Qrhonks(i) &
                        + QRij(j,kk1)*uVunk(lma2,ib)*uVunk(lma1,ib)
#else
                   Qrhonks(i) = Qrhonks(i) &
                        + QRij(j,kk1)*conjg(uVunk(lma1,ib))*uVunk(lma2,ib)
                   Qrhonks(i) = Qrhonks(i) &
                        + QRij(j,kk1)*conjg(uVunk(lma2,ib))*uVunk(lma1,ib)
#endif
                end do ! j

             end if

          end do ! kk1

       end do ! ib

!===== get Qrhonks =====

!----- total = term1 + term2 -----
       rhonks(nn1:nn2) = rhonks(nn1:nn2) + dble( Qrhonks(nn1:nn2) )
!===== total = term1 + term2 =====

       deallocate( uVunk  )

    end select

  END SUBROUTINE get_rhonks

END MODULE rhonks_g_module
