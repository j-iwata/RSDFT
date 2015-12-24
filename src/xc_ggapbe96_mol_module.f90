MODULE xc_ggapbe96_mol_module

  use parallel_module
  use rgrid_mol_module
  use bc_module, only: www, bcset

  implicit none

  PRIVATE
  PUBLIC :: get_LLL_mol, construct_gradient_mol, calc_ve_mol

CONTAINS


  SUBROUTINE get_LLL_mol( mx,my,mz,LLL )
    implicit none
    integer,intent(IN)  :: mx,my,mz
    integer,intent(OUT) :: LLL(-mx:mx,-my:my,-mz:mz)
    integer :: i, ML_0, ML_1, ierr

    LLL=0

    ML_0 = id_grid(myrank_g) + 1
    ML_1 = id_grid(myrank_g) + ir_grid(myrank_g)

    do i=ML_0,ML_1
       LLL( LL(1,i),LL(2,i),LL(3,i) ) = i
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, LLL, size(LLL), MPI_INTEGER &
         , MPI_SUM, comm_grid, ierr )

  END SUBROUTINE get_LLL_mol


  SUBROUTINE construct_gradient_mol( Md,nab,rho,gx,gy,gz )
    implicit none
    integer,intent(IN)  :: Md
    real(8),intent(IN)  :: nab(-Md:Md), rho(:,:)
    real(8),intent(OUT) :: gx(:),gy(:),gz(:)
    integer :: s,i,m,MSP,ML_0,ML_1,i0,i1,i2,i3
    real(8) :: c

    ML_0 = id_grid(myrank_g) + 1
    ML_1 = id_grid(myrank_g) + ir_grid(myrank_g)
    MSP  = sum(ir_spin)

    www(:,:,:,:)=0.0d0
    do s=1,MSP
       do i=ML_0,ML_1
          i0 = i-ML_0+1
          i1 = LL(1,i)
          i2 = LL(2,i)
          i3 = LL(3,i)
          www(i1,i2,i3,1) = www(i1,i2,i3,1) + rho(i0,s)
       end do
    end do

    call bcset(1,1,Md,0)

    c = 1.0d0/Hsize

    gx(:)=0.0d0
    gy(:)=0.0d0
    gz(:)=0.0d0

    i=0
    do i=ML_0,ML_1
       i0=i-ML_0+1
       i1 = LL(1,i)
       i2 = LL(2,i)
       i3 = LL(3,i)
       do m=1,Md
          gx(i0) = gx(i0) - c*nab(m)*( www(i1-m,i2,i3,1) - www(i1+m,i2,i3,1) )
          gy(i0) = gy(i0) - c*nab(m)*( www(i1,i2-m,i3,1) - www(i1,i2+m,i3,1) )
          gz(i0) = gz(i0) - c*nab(m)*( www(i1,i2,i3-m,1) - www(i1,i2,i3+m,1) )
       end do
    end do

  END SUBROUTINE construct_gradient_mol


  SUBROUTINE calc_ve_mol( ML_0, ML_1, m1,m2,m3, Md, nab, ve, LLL, rrrr )
    implicit none
    integer,intent(IN) :: ML_0,ML_1,m1,m2,m3,Md,LLL(-m1:m1,-m2:m2,-m3:m3)
    real(8),intent(IN) :: nab(-Md:Md),rrrr(:,:)
    real(8),intent(INOUT) :: ve(ML_0:)
    real(8) :: c,cm
    integer :: i,i1,i2,i3,j,m

    c = 1.0d0/Hsize

    do i=ML_0,ML_1

       i1 = LL(1,i)
       i2 = LL(2,i)
       i3 = LL(3,i)

       do m=-Md,Md

          cm=-nab(m)*sign(1,m)*c ! The minus sign comes from anti symmetry
                                 !  of nabla matrix (Dji=-Dij). See XC.doc.

          j = LLL(i1+m,i2,i3)
          if ( j /= 0 ) ve(i) = ve(i) + cm*rrrr(j,1)

          j = LLL(i1,i2+m,i3)
          if ( j /= 0 ) ve(i) = ve(i) + cm*rrrr(j,2)

          j = LLL(i1,i2,i3+m)
          if ( j /= 0 ) ve(i) = ve(i) + cm*rrrr(j,3)

       end do ! m

    end do ! i

  END SUBROUTINE calc_ve_mol


END MODULE xc_ggapbe96_mol_module
