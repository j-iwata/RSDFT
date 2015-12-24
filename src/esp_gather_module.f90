MODULE esp_gather_module

  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: esp_gather

CONTAINS

  SUBROUTINE esp_gather(MB,MBZ,MSP,esp)
    integer,intent(IN) :: MB,MBZ,MSP
    real(8),intent(IN) :: esp(MB,MBZ,MSP)
    integer :: k,s,n,ierr,np
    integer :: MB_0,MBZ_0,MBZ_1,MSP_0,MSP_1
    integer,allocatable :: ir(:),id(:)

    call write_border( 1, " esp_gather(start)" )

    np = max( np_band,np_bzsm,np_spin )

    MB_0  = id_band(myrank_b)+1
    MBZ_0 = id_bzsm(myrank_k)+1
    MBZ_1 = id_bzsm(myrank_k)+ir_bzsm(myrank_k)
    MSP_0 = id_spin(myrank_s)+1
    MSP_1 = id_spin(myrank_s)+ir_spin(myrank_s)

    allocate( ir(0:np-1),id(0:np-1) )

    id(0:np_band-1)=id_band(0:np_band-1)
    ir(0:np_band-1)=ir_band(0:np_band-1)

    do s=MSP_0,MSP_1
    do k=MBZ_0,MBZ_1
       call mpi_allgatherv(esp(MB_0,k,s),ir(myrank_b),mpi_real8, &
                           esp(1,k,s),ir,id,mpi_real8,comm_band,ierr)
    end do
    end do

    id(0:np_bzsm-1)=id_bzsm(0:np_bzsm-1)*MB
    ir(0:np_bzsm-1)=ir_bzsm(0:np_bzsm-1)*MB

    do s=MSP_0,MSP_1
       call mpi_allgatherv(esp(1,MBZ_0,s),ir(myrank_k),mpi_real8, &
                           esp(1,1,s),ir,id,mpi_real8,comm_bzsm,ierr)
    end do

    id(0:np_spin-1)=id_spin(0:np_spin-1)*MB*MBZ
    ir(0:np_spin-1)=ir_spin(0:np_spin-1)*MB*MBZ

    call mpi_allgatherv(esp(1,1,MSP_0),ir(myrank_s),mpi_real8, &
                        esp(1,1,1),ir,id,mpi_real8,comm_spin,ierr)

    deallocate( id,ir )

    call write_border( 1, " esp_gather(end)" )

  END SUBROUTINE esp_gather

END MODULE esp_gather_module
