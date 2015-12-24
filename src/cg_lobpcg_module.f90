MODULE cg_lobpcg_module

  use hamiltonian_module
  use cgpc_module

  implicit none

  PRIVATE
  PUBLIC :: init_lobpcg, lobpcg

  integer :: ML_0,ML_1
  integer :: MB_0,MB_1
  integer :: MB_d
  real(8) :: dV
  integer :: comm_grid

  include 'mpif.h'

CONTAINS


  SUBROUTINE init_lobpcg( n1,n2, m1,m2, dV_in, MB_d_in, comm_in )
    implicit none
    integer,intent(IN) :: n1,n2,m1,m2,MB_d_in,comm_in
    real(8),intent(IN) :: dV_in
    ML_0 = n1
    ML_1 = n2
    MB_0 = m1
    MB_1 = m2
    dV   = dV_in
    MB_d = MB_d_in
    comm_grid = comm_in
  END SUBROUTINE init_lobpcg


#ifdef _DRSDFT_
  SUBROUTINE lobpcg( k, s, Mcg, igs, unk, esp, res )

    implicit none

    integer,intent(IN) :: k,s,Mcg,igs
    real(8),intent(INOUT) :: unk(:,:)
    real(8),intent(INOUT) :: esp(:),res(:)

    integer :: ns,ne,nn,n,m,icg,n1,n2,ierr
    integer :: mm,i,ML0,ld,j
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: work(999),W(999),c,d,r
    real(8),allocatable :: sb(:),rb(:)
    real(8),allocatable :: E0(:),E1(:)
    real(8) :: ztmp,zdV

    real(8),allocatable :: vv(:,:,:),hv(:,:,:)
    real(8),allocatable :: vt(:,:),ut(:,:),wt(:,:)
    real(8),allocatable :: vv0(:,:,:),vv1(:,:,:)
    real(8),parameter :: zero=0.d0,one=1.d0

    n1  = ML_0
    n2  = ML_1
    ML0 = ML_1-ML_0+1
    ld  = 3*MB_d

    allocate( sb(MB_d) )
    allocate( rb(MB_d) )
    allocate( E0(MB_d) )
    allocate( E1(MB_d) )
    allocate( vv(n1:n2,MB_d,4) )
    allocate( hv(n1:n2,MB_d,3) )
    allocate( vt(n1:n2,MB_d*3) )
    allocate( ut(n1:n2,MB_d*3) )
    allocate( wt(n1:n2,MB_d) )
    allocate( vv0(ld,ld,2) )
    allocate( vv1(ld,ld,2) )
    vv=zero
    hv=zero

    res(:) = 0.0d0
    esp(:) = 0.0d0

    do ns=MB_0,MB_1,MB_d

       ne = min( ns+MB_d-1, MB_1 )
       nn = ne - ns + 1

       E0(1:nn) = 1.d10

       vv(:,1:nn,1) = unk(:,ns:ne) ! xk

       call hamiltonian( k, s, vv, hv, n1, n2, ns, ne )

       do n=1,nn
          sb(n)=sum( vv(:,n,1)*hv(:,n,1) )*dV
       end do

       call mpi_allreduce(sb,E1,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       do n=1,nn
          do i=n1,n2
             vv(i,n,4) = -( hv(i,n,1) - E1(n)*vv(i,n,1) ) ! gk
          end do
          sb(n)=sum( abs(vv(:,n,4))**2 )*dV
       end do

       call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

!
! --- CG iteration Start
!
       do icg=1,Mcg+1

! --- Convergence check ---

          res(ns:ne)=rb(1:nn)

          if ( all(rb(1:nn) < ep0) ) exit
          if ( all(abs(E1(1:nn)-E0(1:nn)) < ep1) ) exit
          if ( icg == Mcg+1 ) exit

! ---

          do n=1,nn
             do i=n1,n2
                vv(i,n,3)=vv(i,n,4) ! Pgk=gk
             end do
          end do

! --- Preconditioning ---

          call preconditioning(E1,k,s,nn,ML0,vv(n1,1,1),vv(n1,1,4),vv(n1,1,3))

! ---

          call hamiltonian(k,s,vv(n1,1,3),hv(n1,1,3),n1,n2,ns,ne)

          if ( icg == 1 ) then
             mm = 2*nn
             vv(:,1:nn,2) = vv(:,1:nn,3)
             hv(:,1:nn,2) = hv(:,1:nn,3)
          else
             mm = 3*nn
          end if

          m = 0
          do j=1,3
             do n=1,nn
                m = m + 1
                vt(:,m) = vv(:,n,j)
                ut(:,m) = hv(:,n,j)
             end do
          end do
          vv0(:,:,:)=zero
          zdV = dV
          ztmp=0.5d0*dV
          call dsyrk( 'U','C',mm,ML0,zdV,vt,ML0,zero,vv0(1,1,1),ld)
          call dsyr2k('U','C',mm,ML0,ztmp,vt,ML0,ut,ML0,zero,vv0(1,1,2),ld)
          call mpi_allreduce &
               (vv0,vv1,ld*ld*2,MPI_REAL8,mpi_sum,comm_grid,ierr)

          call dsygv &
          (1,'V','U',mm,vv1(1,1,2),ld,vv1(1,1,1),ld,W,work,999,ierr)

          E0(1:nn) = E1(1:nn)
          E1(1:nn) = W(1:nn)

          call dgemm('N','N',ML0,nn,mm,one,vt,ML0,vv1(1,1,2),ld,zero,wt,ML0)
          do n=1,nn
             vv(:,n,1) = wt(:,n) ! xk
          end do

          call dgemm('N','N',ML0,nn,mm,one,ut,ML0,vv1(1,1,2),ld,zero,wt,ML0)
          do n=1,nn
             hv(:,n,1) = wt(:,n) ! hxk
          end do

          wt(n1:n2,1:nn) = matmul( vt(n1:n2,nn+1:mm),vv1(nn+1:mm,nn+1:2*nn,2) )
          do n=1,nn
             vv(:,n,2) = wt(:,n) ! pk
          end do

          wt(n1:n2,1:nn) = matmul( ut(n1:n2,nn+1:mm),vv1(nn+1:mm,nn+1:2*nn,2) )
          do n=1,nn
             hv(:,n,2) = wt(:,n) ! hpk
          end do

          do n=1,nn
             vv(:,n,4) = -( hv(:,n,1) - E1(n)*vv(:,n,1) )
             sb(n) = sum( abs(vv(:,n,4))**2 )*dV
          end do
          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       end do ! icg

       esp(ns:ne)   = E1(1:nn)
       unk(:,ns:ne) = vv(:,1:nn,1)

    end do ! ns [band-loop]

    deallocate( vv1,vv0 )
    deallocate( wt )
    deallocate( ut,vt )
    deallocate( vv,hv )
    deallocate( E0,E1 )
    deallocate( sb,rb )

    return

  END SUBROUTINE lobpcg

#else

  SUBROUTINE lobpcg( k, s, Mcg, igs, unk, esp, res )

    implicit none

    integer,intent(IN) :: k,s,Mcg,igs
    complex(8),intent(INOUT) :: unk(:,:)
    real(8),intent(INOUT) :: esp(:),res(:)

    integer :: ns,ne,nn,n,m,icg,n1,n2,ierr
    integer :: mm,i,ML0,ld,j
    real(8),parameter :: ep0=0.d0
    real(8),parameter :: ep1=1.d-15
    real(8) :: rwork(999),W(999),c,d,r
    real(8),allocatable :: sb(:),rb(:)
    real(8),allocatable :: E0(:),E1(:)
    complex(8) :: work(999),ztmp,zdV

    complex(8),allocatable :: vv(:,:,:),hv(:,:,:)
    complex(8),allocatable :: vt(:,:),ut(:,:),wt(:,:)
    complex(8),allocatable :: vv0(:,:,:),vv1(:,:,:)
    complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)

    n1  = ML_0
    n2  = ML_1
    ML0 = ML_1-ML_0+1
    ld  = 3*MB_d

    allocate( sb(MB_d) )
    allocate( rb(MB_d) )
    allocate( E0(MB_d) )
    allocate( E1(MB_d) )
    allocate( vv(n1:n2,MB_d,4) )
    allocate( hv(n1:n2,MB_d,3) )
    allocate( vt(n1:n2,MB_d*3) )
    allocate( ut(n1:n2,MB_d*3) )
    allocate( wt(n1:n2,MB_d) )
    allocate( vv0(ld,ld,2) )
    allocate( vv1(ld,ld,2) )
    vv=zero
    hv=zero

    res(:) = 0.0d0
    esp(:) = 0.0d0

    do ns=MB_0,MB_1,MB_d

       ne = min( ns+MB_d-1, MB_1 )
       nn = ne - ns + 1

       E0(1:nn) = 1.d10

       vv(:,1:nn,1) = unk(:,ns:ne) ! xk

       call hamiltonian( k, s, vv, hv, n1, n2, ns, ne )

       do n=1,nn
          sb(n)=sum( conjg(vv(:,n,1))*hv(:,n,1) )*dV
       end do

       call mpi_allreduce(sb,E1,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       do n=1,nn
          do i=n1,n2
             vv(i,n,4) = -( hv(i,n,1) - E1(n)*vv(i,n,1) ) ! gk
          end do
          sb(n)=sum( abs(vv(:,n,4))**2 )*dV
       end do

       call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

!
! --- CG iteration Start
!
       do icg=1,Mcg+1

! --- Convergence check ---

          res(ns:ne)=rb(1:nn)

          if ( all(rb(1:nn) < ep0) ) exit
          if ( all(abs(E1(1:nn)-E0(1:nn)) < ep1) ) exit
          if ( icg == Mcg+1 ) exit

! ---

          do n=1,nn
             do i=n1,n2
                vv(i,n,3)=vv(i,n,4) ! Pgk=gk
             end do
          end do

! --- Preconditioning ---

          call preconditioning(E1,k,s,nn,ML0,vv(n1,1,1),vv(n1,1,4),vv(n1,1,3))

! ---

          call hamiltonian(k,s,vv(n1,1,3),hv(n1,1,3),n1,n2,ns,ne)

          if ( icg == 1 ) then
             mm = 2*nn
             vv(:,1:nn,2) = vv(:,1:nn,3)
             hv(:,1:nn,2) = hv(:,1:nn,3)
          else
             mm = 3*nn
          end if

          m = 0
          do j=1,3
             do n=1,nn
                m = m + 1
                vt(:,m) = vv(:,n,j)
                ut(:,m) = hv(:,n,j)
             end do
          end do
          vv0(:,:,:)=zero
          zdV = dV
          ztmp=0.5d0*dV
          call zherk( 'U','C',mm,ML0,zdV,vt,ML0,zero,vv0(1,1,1),ld)
          call zher2k('U','C',mm,ML0,ztmp,vt,ML0,ut,ML0,zero,vv0(1,1,2),ld)
          call mpi_allreduce &
               (vv0,vv1,ld*ld*2,MPI_COMPLEX16,mpi_sum,comm_grid,ierr)

          call zhegv &
          (1,'V','U',mm,vv1(1,1,2),ld,vv1(1,1,1),ld,W,work,999,rwork,ierr)

          E0(1:nn) = E1(1:nn)
          E1(1:nn) = W(1:nn)

          call zgemm('N','N',ML0,nn,mm,one,vt,ML0,vv1(1,1,2),ld,zero,wt,ML0)
          do n=1,nn
             vv(:,n,1) = wt(:,n) ! xk
          end do

          call zgemm('N','N',ML0,nn,mm,one,ut,ML0,vv1(1,1,2),ld,zero,wt,ML0)
          do n=1,nn
             hv(:,n,1) = wt(:,n) ! hxk
          end do

          wt(n1:n2,1:nn) = matmul( vt(n1:n2,nn+1:mm),vv1(nn+1:mm,nn+1:2*nn,2) )
          do n=1,nn
             vv(:,n,2) = wt(:,n) ! pk
          end do

          wt(n1:n2,1:nn) = matmul( ut(n1:n2,nn+1:mm),vv1(nn+1:mm,nn+1:2*nn,2) )
          do n=1,nn
             hv(:,n,2) = wt(:,n) ! hpk
          end do

          do n=1,nn
             vv(:,n,4) = -( hv(:,n,1) - E1(n)*vv(:,n,1) )
             sb(n) = sum( abs(vv(:,n,4))**2 )*dV
          end do
          call mpi_allreduce(sb,rb,nn,mpi_real8,mpi_sum,comm_grid,ierr)

       end do ! icg

       esp(ns:ne)   = E1(1:nn)
       unk(:,ns:ne) = vv(:,1:nn,1)

    end do ! ns [band-loop]

    deallocate( vv1,vv0 )
    deallocate( wt )
    deallocate( ut,vt )
    deallocate( vv,hv )
    deallocate( E0,E1 )
    deallocate( sb,rb )

    return

  END SUBROUTINE lobpcg
#endif


END MODULE cg_lobpcg_module
