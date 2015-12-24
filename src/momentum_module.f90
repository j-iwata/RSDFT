MODULE momentum_module

  use rgrid_module
  use parallel_module
  use bz_module
  use aa_module
  use bb_module
  use kinetic_variables, only: Md, coef_nab
  use bc_module
  use ps_nloc2_variables
  use atom_module
  use pseudopot_module, only: pselect

  implicit none

  PRIVATE
  PUBLIC :: calc_expectval_momentum

CONTAINS


  SUBROUTINE calc_expectval_momentum(k,n1,n2,b1,b2,tpsi,pxyz)
    implicit none
    integer,intent(IN)    :: k,n1,n2,b1,b2
    complex(8),intent(IN) :: tpsi(n1:n2,b1:b2)
    real(8),intent(OUT)   :: pxyz(3,b1:b2)
    integer :: ib,ierr
    real(8) :: kxyz(3),pxyz0(3,b1:b2)
    pxyz0=0.d0
    pxyz =0.d0
    call momentum_kine(n1,n2,b1,b2,tpsi,pxyz)
    pxyz0=pxyz0+pxyz
    pxyz=0.d0
    select case( pselect )
    case( 2 )
       call momentum_nloc(k,n1,n2,b1,b2,tpsi,pxyz)
    case( 3 )
!       call momentum_ps_nloc3(k,n1,n2,b1,b2,tpsi,pxyz)
    end select
    pxyz0=pxyz0+pxyz
    call mpi_allreduce(pxyz0,pxyz,3*(b2-b1+1),mpi_real8,mpi_sum,comm_grid,ierr)

    kxyz(1) = bb(1,1)*kbb(1,k) + bb(1,2)*kbb(2,k) + bb(1,3)*kbb(3,k)
    kxyz(2) = bb(2,1)*kbb(1,k) + bb(2,2)*kbb(2,k) + bb(2,3)*kbb(3,k)
    kxyz(3) = bb(3,1)*kbb(1,k) + bb(3,2)*kbb(2,k) + bb(3,3)*kbb(3,k)

    do ib=b1,b2
       pxyz(1:3,ib) = pxyz(1:3,ib) + kxyz(1:3)
!       if ( myrank == 0 ) write(*,'(1x,5x,2x,3g22.12)') kxyz(1:3)
    end do

  END SUBROUTINE calc_expectval_momentum


  SUBROUTINE momentum_kine(n1,n2,b1,b2,tpsi,pxyz)
    implicit none
    integer,intent(IN)    :: n1,n2,b1,b2
    complex(8),intent(IN) :: tpsi(n1:n2,b1:b2)
    real(8),intent(INOUT) :: pxyz(3,b1:b2)
    integer :: n,nb,ib,i,i1,i2,i3,j,a1b,b1b,a2b,b2b,a3b,b3b
    real(8) :: c,a1,a2,a3,a2x_nab(3),a2y_nab(3),a2z_nab(3),pi2
    complex(8),allocatable :: uuu(:,:)
    complex(8),parameter :: zi=(0.d0,1.d0)

    a1b = Igrid(1,1)
    b1b = Igrid(2,1)
    a2b = Igrid(1,2)
    b2b = Igrid(2,2)
    a3b = Igrid(1,3)
    b3b = Igrid(2,3)

    pi2 = 2.d0*acos(-1.d0)
    a1  = sqrt(sum(aa(1:3,1)**2))/pi2
    a2  = sqrt(sum(aa(1:3,2)**2))/pi2
    a3  = sqrt(sum(aa(1:3,3)**2))/pi2
    a2x_nab(1) = bb(1,1)*a1
    a2x_nab(2) = bb(1,2)*a2
    a2x_nab(3) = bb(1,3)*a3
    a2y_nab(1) = bb(2,1)*a1
    a2y_nab(2) = bb(2,2)*a2
    a2y_nab(3) = bb(2,3)*a3
    a2z_nab(1) = bb(3,1)*a1
    a2z_nab(2) = bb(3,2)*a2
    a2z_nab(3) = bb(3,3)*a3

    nb=b2-b1+1

    do n=1,nb
       ib=n-1+b1
       i=n1-1
       do i3=a3b,b3b
       do i2=a2b,b2b
       do i1=a1b,b1b
          i=i+1
          www(i1,i2,i3,n) = tpsi(i,ib)
       end do
       end do
       end do
    end do

    call bcset(1,1,Md,0)

    allocate( uuu(3,nb) )
    uuu=0.d0

    do n=1,nb
       ib=n-1+b1
       do j=1,Md
          c=-coef_nab(1,j)
          i=n1-1
          do i3=a3b,b3b
          do i2=a2b,b2b
          do i1=a1b,b1b
             i=i+1
             uuu(1,n)=uuu(1,n)+conjg( tpsi(i,ib) )*c*( www(i1-j,i2,i3,n)-www(i1+j,i2,i3,n) )
          end do
          end do
          end do
          c=-coef_nab(2,j)
          i=n1-1
          do i3=a3b,b3b
          do i2=a2b,b2b
          do i1=a1b,b1b
             i=i+1
             uuu(2,n)=uuu(2,n)+conjg( tpsi(i,ib) )*c*( www(i1,i2-j,i3,n)-www(i1,i2+j,i3,n) )
          end do
          end do
          end do
          c=-coef_nab(3,j)
          i=n1-1
          do i3=a3b,b3b
          do i2=a2b,b2b
          do i1=a1b,b1b
             i=i+1
             uuu(3,n)=uuu(3,n)+conjg( tpsi(i,ib) )*c*( www(i1,i2,i3-j,n)-www(i1,i2,i3+j,n) )
          end do
          end do
          end do
       end do ! j
    end do ! n

    uuu(:,:) = -zi*uuu(:,:)*dV

    do n=1,nb
       ib=n-1+b1
       pxyz(1,ib)=pxyz(1,ib)+a2x_nab(1)*uuu(1,n)+a2x_nab(2)*uuu(2,n)+a2x_nab(3)*uuu(3,n)
       pxyz(2,ib)=pxyz(2,ib)+a2y_nab(1)*uuu(1,n)+a2y_nab(2)*uuu(2,n)+a2y_nab(3)*uuu(3,n)
       pxyz(3,ib)=pxyz(3,ib)+a2z_nab(1)*uuu(1,n)+a2z_nab(2)*uuu(2,n)+a2z_nab(3)*uuu(3,n)
!       if ( myrank == 0 ) write(*,'(1x,5x,2x,3g22.12)') pxyz(:,ib)
    end do

    deallocate( uuu )

  END SUBROUTINE momentum_kine


  SUBROUTINE momentum_nloc(k,n1,n2,b1,b2,tpsi,pxyz)
    implicit none
    integer,intent(IN)    :: k,n1,n2,b1,b2
    complex(8),intent(IN) :: tpsi(n1:n2,b1:b2)
    real(8),intent(INOUT) :: pxyz(3,b1:b2)
    integer :: lma,i,j,ib,i1,i2,i3,nb,a,k1,k2,k3,m,nreq,jrank,irank,ierr
    integer,allocatable :: istatus(:,:),ireq(:)
    complex(8),allocatable :: uuu(:,:,:),uuu0(:,:,:)
    logical,allocatable :: a_rank(:)

    nb=b2-b1+1

    allocate( uuu(0:3,nzlma,b1:b2)  ) ; uuu=0.0d0
    allocate( uuu0(0:3,nzlma,b1:b2) ) ; uuu0=0.0d0

    do ib=b1,b2
       do lma=1,nzlma
          do j=1,MJJ(lma)
             i=JJP(j,lma)
#ifndef _DRSDFT_
             uuu(0,lma,ib)=uuu(0,lma,ib)+conjg( uVk(j,lma,k) )*tpsi(i,ib)
             uuu(1,lma,ib)=uuu(1,lma,ib)+conjg( xVk(j,lma,k) )*tpsi(i,ib)
             uuu(2,lma,ib)=uuu(2,lma,ib)+conjg( yVk(j,lma,k) )*tpsi(i,ib)
             uuu(3,lma,ib)=uuu(3,lma,ib)+conjg( zVk(j,lma,k) )*tpsi(i,ib)
#endif
          end do ! j
          uuu(0,lma,ib)=iuV(lma)*uuu(0,lma,ib)*dV*dV
       end do ! lma
    end do ! ib

    nreq = 2*maxval( nrlma_xyz )
    allocate( ireq(nreq) )
    allocate( istatus(MPI_STATUS_SIZE,nreq) )

    nreq=0
    do i=1,6
       select case(i)
       case(1,3,5)
          j=i+1
          uuu0(:,:,:)=uuu(:,:,:)
       case(2,4,6)
          j=i-1
       end select
       do m=1,nrlma_xyz(i)
          nreq=0
          irank=num_2_rank(m,i)
          jrank=num_2_rank(m,j)
          if( irank>=0 )then
             i2=0
             do ib=b1,b2
                do i1=1,lma_nsend(irank)
                   do i3=0,3
                      i2=i2+1
                      sbufnl(i2,irank)=uuu0(i3,sendmap(i1,irank),ib)
                   end do
                end do
             end do
             nreq=nreq+1
             call mpi_isend(sbufnl(1,irank),lma_nsend(irank)*nb*4 &
                  ,TYPE_MAIN,irank,1,comm_grid,ireq(nreq),ierr)
          end if
          if ( jrank >= 0 ) then
             nreq=nreq+1
             call mpi_irecv(rbufnl(1,jrank),lma_nsend(jrank)*nb*4 &
                  ,TYPE_MAIN,jrank,1,comm_grid,ireq(nreq),ierr)
          end if
          call mpi_waitall(nreq,ireq,istatus,ierr)
          if ( jrank >= 0 ) then
             i2=0
             do ib=b1,b2
                do i1=1,lma_nsend(jrank)
                   do i3=0,3
                      i2=i2+1
                      uuu(i3,recvmap(i1,jrank),ib) &
                         = uuu(i3,recvmap(i1,jrank),ib) + rbufnl(i2,jrank)
                   end do ! i3
                end do ! i1
             end do ! ib
          end if
       end do ! m
    end do ! i

    deallocate( istatus )
    deallocate( ireq )

    allocate( a_rank(Natom) )
    a_rank(:)=.false.
    do a=1,Natom
       i1 = nint( aa_atom(1,a)*Ngrid(1) )
       i2 = nint( aa_atom(2,a)*Ngrid(2) )
       i3 = nint( aa_atom(3,a)*Ngrid(3) )
       k1 = i1/Ngrid(1) ; if ( i1<0 ) k1=(i1+1)/Ngrid(1)-1
       k2 = i2/Ngrid(2) ; if ( i2<0 ) k2=(i2+1)/Ngrid(2)-1
       k3 = i3/Ngrid(3) ; if ( i3<0 ) k3=(i3+1)/Ngrid(3)-1
       i1 = i1 - k1*Ngrid(1)
       i2 = i2 - k2*Ngrid(2)
       i3 = i3 - k3*Ngrid(3)
       if ( Igrid(1,1) <= i1 .and. i1 <= Igrid(2,1) .and. &
            Igrid(1,2) <= i2 .and. i2 <= Igrid(2,2) .and. &
            Igrid(1,3) <= i3 .and. i3 <= Igrid(2,3) ) then
          a_rank(a)=.true.
       end if
    end do

    do ib=b1,b2
    do lma=1,nzlma
       a=amap(lma)
       if ( a <= 0 ) cycle
       if ( a_rank(a) ) then
          pxyz(1,ib)=pxyz(1,ib)+2.d0*aimag( conjg(uuu(1,lma,ib))*uuu(0,lma,ib) )
          pxyz(2,ib)=pxyz(2,ib)+2.d0*aimag( conjg(uuu(2,lma,ib))*uuu(0,lma,ib) )
          pxyz(3,ib)=pxyz(3,ib)+2.d0*aimag( conjg(uuu(3,lma,ib))*uuu(0,lma,ib) )
       end if
    end do
!    if ( myrank == 0 ) write(*,'(1x,5x,2x,3g22.12)') pxyz(:,ib)
    end do

    deallocate( a_rank )

    deallocate( uuu0 )
    deallocate( uuu  )

  END SUBROUTINE momentum_nloc


END MODULE momentum_module
