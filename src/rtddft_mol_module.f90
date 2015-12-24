MODULE rtddft_mol_module

  use wf_module
  use rgrid_mol_module, only: Hsize, LL
  use grid_module, only: grid, get_range_rgrid
  use hamiltonian_module
  use parallel_module
  use density_module
  use total_energy_module
  use hartree_module
  use xc_module
  use watch_module
  use io_module

  implicit none

  PRIVATE
  PUBLIC :: rtddft_mol, init_rtddft_mol

  type td
     real(8) :: dt
     real(8) :: tmax
     integer :: nt
     integer :: ialg
     integer :: nalg
     real(8) :: field(3)
     real(8) :: strength
  END type td

  type(td) :: tddft
  real(8),allocatable :: dipole(:,:)
  integer,parameter :: unit=91
  character(10),parameter :: file_tddft = 'tddft_data'

  real(8),parameter :: Tau2fs = 2.418884326505d-2

CONTAINS


  SUBROUTINE init_rtddft_mol( unit, rank )
     implicit none
     integer,intent(IN) :: unit, rank
     integer :: i
     real(8) :: sbuf(7)
     character(5) :: cbuf,ckey
     call write_border(40," init_rtddft_mol")
     tddft%dt  =0.0d0
     tddft%nt  =0
     tddft%ialg=1
     tddft%nalg=4
     tddft%field=0.0d0
     tddft%strength=0.0d0
     if ( rank == 0 ) then
        rewind unit
        do i=1,100000
           read(unit,*,END=900) cbuf
           call convert_capital(cbuf,ckey)
           if ( ckey == "TDDFT" ) then
              backspace(unit)
              read(unit,*) cbuf, tddft%dt, tddft%nt, tddft%ialg, tddft%nalg
           else if ( ckey == "FIELD" ) then
              backspace(unit)
              read(unit,*) cbuf, tddft%field(1:3)
           end if
        end do ! i
900     continue
     end if
     sbuf(1)=tddft%dt
     sbuf(2)=tddft%nt
     sbuf(3)=tddft%ialg
     sbuf(4)=tddft%nalg
     sbuf(5:7)=tddft%field(1:3)
     call MPI_BCAST(sbuf,7,MPI_REAL8,0,MPI_COMM_WORLD,i)
     tddft%dt         = sbuf(1)
     tddft%nt         = nint( sbuf(2) )
     tddft%ialg       = nint( sbuf(3) )
     tddft%nalg       = nint( sbuf(4) )
     tddft%field(1:3) = sbuf(5:7)

     tddft%tmax = tddft%dt * tddft%nt
     tddft%strength = sqrt(sum(tddft%field(:)**2))

     if ( rank == 0 ) then
        write(*,*) "dt=",tddft%dt
        write(*,*) "nt=",tddft%nt
        write(*,*) "tmax=",tddft%tmax
        write(*,*) "ialgorithm=",tddft%ialg
        write(*,*) "nalgorithm=",tddft%nalg
        write(*,'("field=",3f15.8)') tddft%field(1:3)
        write(*,'("strength=",3f15.8)') tddft%strength
     end if
  END SUBROUTINE init_rtddft_mol


  SUBROUTINE rtddft_mol( iswitch_tddft )
    implicit none
    integer,intent(IN) :: iswitch_tddft
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),parameter :: zi=(0.0d0,1.0d0)
    complex(8),allocatable :: tpsi(:),hpsi(:),zcoef(:)
    integer :: itaylor,i,n,k,s,it,ierr
    real(8) :: c,t,ct(0:9),et(0:9)
    type(grid) :: rgrid
    logical :: disp_sw,flag_end
    real(8) :: Etot

#ifdef _DRSDFT_

    write(*,*) "rtddft_mol is available only for COMPLEX16 calculations"
    return

#elif defined _DRSDFT_

    call write_border(40," rtddft_mol(start)")
    call check_disp_switch( disp_sw, 0 )
    call check_disp_switch( .false., 1 )

    ct(:)=0.0d0
    et(:)=0.0d0

    call get_range_rgrid( rgrid )

    allocate( dipole(0:3,0:tddft%nt) ) ; dipole=0.0d0
    allocate( tpsi(ML_0_WF:ML_1_WF) ) ; tpsi=zero
    allocate( hpsi(ML_0_WF:ML_1_WF) ) ; hpsi=zero
    allocate( zcoef(tddft%nalg) ) ; zcoef=zero

    do itaylor=1,tddft%nalg
       c=1.0d0
       do i=1,itaylor
          c=c*i
       end do
       zcoef(itaylor) = (-zi*tddft%dt)**itaylor/c
    end do

! ---

    call Init_IO( "tddft" )

! ---

    if ( myrank == 0 ) open(unit,file=file_tddft,position="append")

! ---

    if ( iswitch_tddft == 1 ) then

       call calc_dipole( dipole(0,0), rgrid%VolumeElement )

       if ( disp_sw ) then
          write(*,'(6f16.10,1x,2f10.5)') &
               0.0,dipole(1:3,0),dipole(0,0),Etot
       end if
       if ( myrank == 0 ) then
          write(unit,'(f10.5,1x,3f20.15,1x,2f20.15)') &
               0.0,dipole(1:3,0),dipole(0,0),Etot
       end if

       call initial_condition

    end if

! ---

    do it=1,tddft%nt

       call watch(ct(0),et(0))

       t = it*tddft%dt

       do s=MS_0_WF,MS_1_WF
       do k=MK_0_WF,MK_1_WF
       do n=MB_0_WF,MB_1_WF

          tpsi(:) = unk(:,n,k,s)
          do itaylor=1,tddft%nalg
             call hamiltonian(k,s,tpsi,hpsi,ML_0_WF,ML_1_WF,1,1)
             unk(:,n,k,s) = unk(:,n,k,s) + zcoef(itaylor)*hpsi(:)
             tpsi(:) = hpsi(:)
          end do

       end do ! n
       end do ! k
       end do ! s

       call calc_density
       call calc_hartree( ML_0_WF, ML_1_WF, MS_WF, rho )
       call calc_xc
       call calc_total_energy( .true., Etot )

       call calc_dipole( dipole(0,it), rgrid%VolumeElement )

       call watch(ct(1),et(1))
       call global_watch( .false., flag_end )

       if ( disp_sw ) then
          write(*,'(6f16.10,1x,2f10.5)') &
               t,dipole(1:3,it),dipole(0,it),Etot,ct(1)-ct(0),et(1)-et(0)
       end if
       if ( myrank == 0 ) then
          write(unit,'(f10.5,1x,3f20.15,1x,2f20.15)') &
               t,dipole(1:3,it),dipole(0,it),Etot
       end if

       if ( flag_end ) exit

    end do ! it

    deallocate( zcoef )
    deallocate( hpsi )
    deallocate( tpsi )
    deallocate( dipole )

! ---

    if ( myrank == 0 ) close(unit)

! ---

    call write_data( disp_sw, .true. )

! ---

#endif

  END SUBROUTINE rtddft_mol


  SUBROUTINE initial_condition
    implicit none
    integer :: i,n,k,s
    real(8) :: x,y,z,kr,kx,ky,kz

    kx = tddft%field(1)
    ky = tddft%field(2)
    kz = tddft%field(3)

    do s=MS_0_WF,MS_1_WF
    do k=MK_0_WF,MK_1_WF
    do n=MB_0_WF,MB_1_WF

       do i=ML_0_WF,ML_1_WF
          x=LL(1,i)*Hsize
          y=LL(2,i)*Hsize
          z=LL(3,i)*Hsize
          kr=x*kx+y*ky+z*kz
          unk(i,n,k,s) = dcmplx(cos(kr),sin(kr))*unk(i,n,k,s)
       end do ! i

    end do ! n
    end do ! k
    end do ! s

  END SUBROUTINE initial_condition


  SUBROUTINE calc_dipole( d, dV )
    implicit none
    real(8),intent(OUT) :: d(0:3)
    real(8),intent(IN) :: dV
    integer :: i
    real(8) :: x,y,z,c,d0(0:3),trho

    c=0.5d0 ; if ( MS_WF == 2 ) c=1.0d0

    d0(:)=0.0d0
    do i=ML_0_WF,ML_1_WF
       x=LL(1,i)*Hsize
       y=LL(2,i)*Hsize
       z=LL(3,i)*Hsize
       trho=c*( rho(i,1) + rho(i,MS_WF) )
       d0(0) = d0(0) + trho
       d0(1) = d0(1) + x*trho
       d0(2) = d0(2) + y*trho
       d0(3) = d0(3) + z*trho
    end do ! i
    d0(:)=d0(:)*dV

    call MPI_ALLREDUCE(d0,d,4,MPI_REAL8,MPI_SUM,comm_grid,i)

  END SUBROUTINE calc_dipole


END MODULE rtddft_mol_module
