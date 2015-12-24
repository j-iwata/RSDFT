MODULE ps_pcc_force_module

  use atom_module
  use bb_module
  use ggrid_module
  use rgrid_module
  use parallel_module
  use xc_module, only: Vxc
  use ps_pcc_module, only: cdcg
  use fft_module

  implicit none

  PRIVATE
  PUBLIC :: calc_ps_pcc_force

CONTAINS


  SUBROUTINE calc_ps_pcc_force( MI, force )
    implicit none
    integer,intent(IN) :: MI
    real(8),intent(OUT) :: force(3,MI)
    integer :: a,ik,i,j,ierr,i1,i2,i3,j1,j2,j3,irank,n
    integer :: MG,ML,ML1,ML2,ML3,N_MI,MI_0,MI_1,MSP_0,MSP_1
    integer,allocatable :: icnt(:),idis(:)
    real(8) :: pi2,a1,a2,a3,Gx,Gy,Gz,Gr
    complex(8),allocatable :: z0(:),z1(:,:,:),z2(:,:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    real(8),allocatable :: w0(:),w1(:),w2(:,:,:)
    include 'mpif.h'

    force(:,:) = 0.0d0

    pi2 = 2.0d0*acos(-1.0d0)
    MG  = NGgrid(0)
    ML  = Ngrid(0)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)

! ---

    allocate( icnt(0:nprocs-1), idis(0:nprocs-1) )
    N_MI = Natom/nprocs
    icnt(0:nprocs-1) = N_MI
    n = Natom - N_MI*nprocs
    if ( n>0 ) then
       do irank=0,n-1
          icnt(irank)=icnt(irank)+1
       end do
    end if
    do irank=0,nprocs-1
       idis(irank) = sum( icnt(0:irank) ) - icnt(irank)
    end do
    MI_0 = idis(myrank)+1
    MI_1 = idis(myrank)+icnt(myrank)
    deallocate( idis, icnt )

! ---

    allocate( w2(0:ML1-1,0:ML2-1,0:ML3-1) ) ; w2=0.0d0
    allocate( w1(ML)                      ) ; w1=0.0d0

    MSP_0 = id_spin(myrank_s) + 1
    MSP_1 = id_spin(myrank_s) + ir_spin(myrank_s)
    n=size(Vxc,1)
    allocate( w0(n) )
    w0(:)=0.0d0
    do j=MSP_0,MSP_1
       w0(:) = w0(:) + Vxc(:,j)
    end do
    call mpi_allreduce( mpi_in_place, w0,n,MPI_REAL8,MPI_SUM,comm_spin,ierr)
    call mpi_allgatherv(w0,ir_grid(myrank_g),MPI_REAL8 &
         ,w1,ir_grid,id_grid,MPI_REAL8,comm_grid,ierr)
    deallocate( w0 )

    i=0
    irank=-1
    do i3=1,node_partition(3)
    do i2=1,node_partition(2)
    do i1=1,node_partition(1)
       irank=irank+1
       do j3=pinfo_grid(5,irank),pinfo_grid(5,irank)+pinfo_grid(6,irank)-1
       do j2=pinfo_grid(3,irank),pinfo_grid(3,irank)+pinfo_grid(4,irank)-1
       do j1=pinfo_grid(1,irank),pinfo_grid(1,irank)+pinfo_grid(2,irank)-1
          i=i+1
          w2(j1,j2,j3)=w1(i)
       end do ! i1
       end do ! i2
       end do ! i3
    end do ! j1
    end do ! j2
    end do ! j3

    deallocate( w1 )

! ---

    call construct_Ggrid(0)

    call init_fft

    allocate( z0(MG) ) ; z0=zero
    allocate( z1(0:ML1-1,0:ML2-1,0:ML3-1) ) ; z1=zero
    allocate( z2(0:ML1-1,0:ML2-1,0:ML3-1) ) ; z2=zero

    do a=MI_0,MI_1

       ik=ki_atom(a)
       a1=pi2*aa_atom(1,a)
       a2=pi2*aa_atom(2,a)
       a3=pi2*aa_atom(3,a)

       do i=1,MG
          j=MGL(i)
          Gr=a1*LLG(1,i)+a2*LLG(2,i)+a3*LLG(3,i)
          z0(i)=cdcg(j,ik)*dcmplx(sin(Gr),cos(Gr))
       end do

       z1(:,:,:)=zero
       do i=1,MG
          i1=mod(ML1+LLG(1,i),ML1)
          i2=mod(ML2+LLG(2,i),ML2)
          i3=mod(ML3+LLG(3,i),ML3)
          Gx=bb(1,1)*LLG(1,i)+bb(1,2)*LLG(2,i)+bb(1,3)*LLG(3,i)
          z1(i1,i2,i3) = Gx*z0(i)
       end do
       call backward_fft( z1, z2 )
       force(1,a) = sum( w2(:,:,:)*z1(:,:,:) )*dV

       z1(:,:,:)=zero
       do i=1,MG
          i1=mod(ML1+LLG(1,i),ML1)
          i2=mod(ML2+LLG(2,i),ML2)
          i3=mod(ML3+LLG(3,i),ML3)
          Gy=bb(2,1)*LLG(1,i)+bb(2,2)*LLG(2,i)+bb(2,3)*LLG(3,i)
          z1(i1,i2,i3) = Gy*z0(i)
       end do
       call backward_fft( z1, z2 )
       force(2,a) = sum( w2(:,:,:)*z1(:,:,:) )*dV

       z1(:,:,:)=zero
       do i=1,MG
          i1=mod(ML1+LLG(1,i),ML1)
          i2=mod(ML2+LLG(2,i),ML2)
          i3=mod(ML3+LLG(3,i),ML3)
          Gz=bb(3,1)*LLG(1,i)+bb(3,2)*LLG(2,i)+bb(3,3)*LLG(3,i)
          z1(i1,i2,i3) = Gz*z0(i)
       end do
       call backward_fft( z1, z2 )
       force(3,a) = sum( w2(:,:,:)*z1(:,:,:) )*dV

    end do ! a

    deallocate( z2 )
    deallocate( z1 )
    deallocate( z0 )

    call finalize_fft

    call destruct_Ggrid

    deallocate( w2 )

    allocate( w2(3,Natom,1) ) ; w2=0.0d0
    w2(1:3,MI_0:MI_1,1) = force(1:3,MI_0:MI_1)
    call MPI_ALLREDUCE(w2,force,3*Natom,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    deallocate( w2 )

  END SUBROUTINE calc_ps_pcc_force


END MODULE ps_pcc_force_module
