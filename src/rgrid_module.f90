MODULE rgrid_module

  use aa_module
  use rgrid_variables
  use rgrid_sol_module
  use rgrid_mol_module
 !use esm_rgrid_module
  use parallel_module

  implicit none

  PRIVATE
  PUBLIC :: Init_Rgrid, InitParallel_Rgrid &
       ,Igrid,Ngrid,Hgrid,dV,zdV

  integer :: SYStype=0

CONTAINS


  SUBROUTINE Init_Rgrid( SYStype_in, Md, unit )
    implicit none
    integer,intent(IN) :: SYStype_in, Md, unit

    call write_border( 80," init_rgrid(start)")

    SYStype = SYStype_in

    select case( SYStype )

    case(0) ! ----- SOL Sol sol

       call Init_RgridSol(aa)

    case(1) ! ----- MOL Mol mol

       call Read_RgridMol(myrank,unit)

       call GetNumGrids_RgridMol(Ngrid)

       call GetSimBox_RgridMol(aa)

       call GetGridSize_RgridMol(Hgrid)

       ax  = aa(1,1)
       Va  = aa(1,1)*aa(2,2)*aa(3,3)
       dV  = Hgrid(1)*Hgrid(2)*Hgrid(3)
       zdV = dV

!    case(2) ! ----- ESM Esm esm
!
!       call Init_RgridSol(aa)
!
!       call Read_RgridESM(myrank,unit)
!
!       call Init_RgridESM(aa,Ngrid,Md)

    end select

    if ( disp_switch_parallel ) then
       write(*,*) "SYStype=",SYStype
       write(*,'(1x,"aa",2x,3f15.10)') aa(1:3,1)
       write(*,'(1x,"  ",2x,3f15.10)') aa(1:3,2)
       write(*,'(1x,"  ",2x,3f15.10)') aa(1:3,3)
       write(*,'(1x,"Ngrid(1:3)=",3i5)') Ngrid(1:3)
       write(*,'(1x,"Ngrid(0)  =",i8 )') Ngrid(0)
       write(*,'(1x,"Hgrid(1:3)=",3f15.8)') Hgrid(1:3)
    end if

    call write_border( 80," init_rgrid(end)")

  END SUBROUTINE Init_Rgrid


  SUBROUTINE InitParallel_Rgrid
    implicit none
    integer :: ierr, Nshift(3)

    call write_border( 80, " InitParallel_Rgrid(start)" )

    Nshift(:) = 0

    select case(SYStype)

    case(0) ! ----- SOL sol

       call InitParallel_RgridSol( node_partition, np_grid, pinfo_grid )

    case(1) ! ----- MOL mol

       call InitParallel_RgridMol(node_partition,np_grid,pinfo_grid,myrank==0)

       Nshift(1:3) = -1

!    case(2)
!
!       call InitParallel_RgridSol( node_partition, np_grid, pinfo_grid )
!       call InitParallel_RgridESM( np_grid, pinfo_grid, Nshift )

    end select

    Igrid(1,0) = pinfo_grid(7,myrank_g) + 1
    Igrid(2,0) = pinfo_grid(7,myrank_g) + pinfo_grid(8,myrank_g)
    Igrid(1,1) = pinfo_grid(1,myrank_g)
    Igrid(2,1) = pinfo_grid(1,myrank_g) + pinfo_grid(2,myrank_g) - 1
    Igrid(1,2) = pinfo_grid(3,myrank_g)
    Igrid(2,2) = pinfo_grid(3,myrank_g) + pinfo_grid(4,myrank_g) - 1
    Igrid(1,3) = pinfo_grid(5,myrank_g)
    Igrid(2,3) = pinfo_grid(5,myrank_g) + pinfo_grid(6,myrank_g) - 1

    Igrid(1,1:3) = Igrid(1,1:3) - Nshift(1:3)
    Igrid(2,1:3) = Igrid(2,1:3) - Nshift(1:3)

    id_grid(:) = pinfo_grid(7,:)
    ir_grid(:) = pinfo_grid(8,:)

    idisp(myrank) = id_grid(myrank_g)
    ircnt(myrank) = ir_grid(myrank_g)
    call mpi_allgather(idisp(myrank),1,mpi_integer,idisp,1,mpi_integer,mpi_comm_world,ierr)
    call mpi_allgather(ircnt(myrank),1,mpi_integer,ircnt,1,mpi_integer,mpi_comm_world,ierr)

    call write_border( 80, " InitParallel_Rgrid(end)" )

  END SUBROUTINE InitParallel_Rgrid


END MODULE rgrid_module
