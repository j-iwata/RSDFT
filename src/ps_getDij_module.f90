MODULE ps_getDij_module

  use parallel_module, only: myrank,nprocs_g ! no need of nprocs_g
  use array_bound_module, only: MSP_0,MSP_1
  use var_ps_member_g, only: Dij00,Dij,N_nzqr,QRij,nzqr_pair
  use var_para_ps_nloc_g, only: MJJ_Q,JJP_Q
  use parallel_module, only: comm_grid
  use localpot_module, only: Vloc
  use rgrid_module, only: dV
  use pseudopot_module, only: pselect
  use ps_nloc2_variables, only: amap,iorbmap,mmap
  use var_ps_member, only: norb, lo
  use atom_module, only: Natom

  implicit none

  PRIVATE
  PUBLIC :: getDij

  include 'mpif.h'

CONTAINS

  SUBROUTINE getDij

    implicit none
    integer :: s,kk1,i,j,i1,i2,m1,m2,lma1,lma2,a1,ierr
    real(8),allocatable :: IntQV(:,:,:,:,:)

    select case ( pselect )
    case default

       return

    case ( 102 )

       call write_border( 1, " getDij(start)" )

       i=maxval( norb )
       j=maxval( lo )
       allocate( IntQV(Natom,i,i,-j:j,-j:j) ) ; IntQV=0.0d0

       do s=MSP_0,MSP_1

          IntQV(:,:,:,:,:)=0.0d0

          do kk1=1,N_nzqr

             lma1 = nzqr_pair(kk1,1)
             lma2 = nzqr_pair(kk1,2)

             a1 = amap(lma1)

             i1 = iorbmap(lma1)
             i2 = iorbmap(lma2)

             m1 = mmap(lma1)
             m2 = mmap(lma2)

             do j=1,MJJ_Q(kk1)
                i=JJP_Q(j,kk1)
                IntQV(a1,i1,i2,m1,m2) = IntQV(a1,i1,i2,m1,m2) &
                     + QRij(j,kk1)*Vloc(i,s)
             end do

          end do ! kk1

          call mpi_allreduce( MPI_IN_PLACE, IntQV, size(IntQV) &
               ,MPI_REAL8, MPI_SUM, comm_grid, ierr )

          do kk1=1,N_nzqr

             lma1 = nzqr_pair(kk1,1)
             lma2 = nzqr_pair(kk1,2)

             a1 = amap(lma1)

             i1 = iorbmap(lma1)
             i2 = iorbmap(lma2)

             m1 = mmap(lma1)
             m2 = mmap(lma2)

             Dij(kk1,s) = Dij00(kk1) + IntQV(a1,i1,i2,m1,m2)*dV

          end do ! kk1

       end do ! s

       deallocate( IntQV )

       call write_border( 1, " getDij(end)" )

    end select

  END SUBROUTINE getDij

END MODULE ps_getDij_module
