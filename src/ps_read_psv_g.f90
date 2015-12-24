MODULE ps_read_PSV_g

  use var_ps_member, only: ps1d, psg_allocate_ps1d

  implicit none

  PRIVATE
  PUBLIC :: readPSVG

CONTAINS

  SUBROUTINE readPSVG( unit_ps,ddi_,qqr_,psi_,phi_,bet_,psp )
    implicit none
    integer,intent(IN) :: unit_ps
    type(ps1d) :: psp
    character(80) :: cbuf9
    integer :: i,j,l,l1,l2,i1,i2,ic
    integer :: k2,ll1,ll2,ll3,ii1,ii2
    real(8) :: dif,r2
    integer,allocatable :: nl3v_(:),l3v_(:,:),ncf(:,:),nrin(:,:)
    real(8),allocatable :: rin(:,:),coe(:,:,:)
    real(8),allocatable :: qrL_(:,:,:)
    real(8),allocatable,intent(IN) :: psi_(:,:,:),phi_(:,:,:),bet_(:,:,:)
    real(8),allocatable,intent(IN) :: ddi_(:,:,:),qqr_(:,:,:)
    ! psi(1:nsmpl-1,Rrefmax,Lrefmax)
    integer :: Rrefmax,Lrefmax,lpsmax,npqmax,ncfmax,nsmpl
    integer :: nl3vmax,npq_,iunit
    integer,parameter :: max_loop=1000000

    Lrefmax = psp%nlf
    Rrefmax = maxval(psp%nrf)
    lpsmax  = Lrefmax*Rrefmax
    npqmax  = lpsmax*(lpsmax+1)/2
    ncfmax  = 10
    nsmpl   = maxval(psp%NRps)-1

    allocate( nl3v_(npqmax)                ) ; nl3v_=0
    allocate( l3v_(lpsmax,npqmax)          ) ; l3v_ =0
    allocate( ncf(lpsmax,npqmax)           ) ; ncf  =0
    allocate( nrin(lpsmax,npqmax)          ) ; nrin =0
    allocate( rin(lpsmax,npqmax)           ) ; rin  =0.d0
    allocate( coe(ncfmax,lpsmax,npqmax)    ) ; coe  =0.d0
    allocate( qrL_(nsmpl+1,lpsmax,npqmax)  ) ; qrL_ =0.d0

    rewind unit_ps

    do j=1,max_loop
      read(unit_ps,'(A)') cbuf9
      if ( index(cbuf9,'#### DATA') > 0 ) exit
    end do
    if ( j > max_loop ) stop "stop@readPSVG(1)" 
    read(unit_ps,*) npq_

    k2=0
    do l1=1,psp%nlf
    do i1=1,psp%nrf(l1)
    do l2=1,l1
    do i2=1,psp%nrf(l2)
       if ( l1 == l2 .and. i2 > i1 ) cycle
       k2=k2+1
       read(unit_ps,*) ll1,ii1,ll2,ii2
       if ( (ll1/=l1) .or. (ll2/=l2) .or. (ii1/=i1) .or. (ii2/=i2) ) then
          stop 'ERROR in pseudization data'
       end if
       read(unit_ps,*) nl3v_(k2)
       do ll3=1,nl3v_(k2)
          read(unit_ps,*) l3v_(ll3,k2),ncf(ll3,k2),nrin(ll3,k2),rin(ll3,k2)
          if ( ncf(ll3,k2) > 0 ) then
             read(unit_ps,*) coe(1:ncf(ll3,k2),ll3,k2)
             dif=abs(psp%rad(nrin(ll3,k2)+1)-rin(ll3,k2))
             if ( dif > 1.d-8 ) then
                stop 'ERROR in PSV read space cutoff'
             end if
          end if
          qrL_(:,ll3,k2)=0.d0
          do ic=ncf(ll3,k2),1,-1
             do i=1,nrin(ll3,k2)
                r2=(psp%rad(i+1))**2
                qrL_(i,ll3,k2)=qrL_(i,ll3,k2)*r2+coe(ic,ll3,k2)
             end do
          end do
          do i=1,nrin(ll3,k2)
             qrL_(i,ll3,k2)=qrL_(i,ll3,k2)*psp%rad(i+1)**(l3v_(ll3,k2)+1)
          end do
          do i=nrin(ll3,k2)+1,nsmpl
             qrL_(i,ll3,k2)=psi_(i,i1,l1)*psi_(i,i2,l2) &
                           -phi_(i,i1,l1)*phi_(i,i2,l2)
          end do
       end do ! ll3
    end do ! i2
    end do ! l2
    end do ! i1
    end do ! l1

    if ( npq_ /= k2 ) stop 'ERROR npq/=k2'

    psp%npq = npqmax  ! npqmax includes too large mergin ( npq_ is enough ) 
    psp%nrf_max = maxval( psp%nrf )
    call psg_allocate_ps1d( psp )

    do l=1,psp%nlf
       do j=1,psp%nrf(l)
       do i=1,j
          psp%ddi(i,j,l) = ddi_(i,j,l)
          psp%ddi(j,i,l) = ddi_(i,j,l)
          psp%qqr(i,j,l) = qqr_(i,j,l)
          psp%qqr(j,i,l) = qqr_(i,j,l)
       end do
       end do
    end do

    k2=0
    do l1=1,psp%nlf
    do i1=1,psp%nrf(l1)
    do l2=1,l1
    do i2=1,psp%nrf(l2)
       if ( l1 == l2 .and. i2 > i1 ) cycle
       k2=k2+1
       psp%nl3v(k2)=nl3v_(k2)
       do ll3=1,nl3v_(k2)
          psp%l3v(ll3,k2)=l3v_(ll3,k2)
          do i=nsmpl,1,-1
             psp%qrL(i+1,ll3,k2)=qrL_(i,ll3,k2)
          end do
       end do ! ll3
    end do ! i2
    end do ! l2
    end do ! i1
    end do ! l1

    deallocate( nl3v_,l3v_ )
    deallocate( ncf,nrin,rin,coe )
    deallocate( qrL_ )

    return

  END SUBROUTINE readPSVG

END MODULE ps_read_PSV_g
