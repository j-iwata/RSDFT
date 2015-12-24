MODULE pseudopot_module

  use var_ps_member
  use ps_read_PSV
  use var_ps_member_g, only: allocatePSG, sendPSG, npq, ddi, qqr, nl3v, l3v, qrL
  use ps_read_TM_module
  use ps_read_YB_module
  use ps_read_UPF_module
  use ps_gth_module
  use io_tools_module

  implicit none

  PRIVATE
  PUBLIC :: ippform,file_ps,inorm,NRps,norb,Mr,lo,no,vql,cdd,cdc,rad &
           ,anorm,viod,Rps,Zps,parloc,rab,cdd_coef,ps_type,Rcloc &
           ,hnml,knml,hnl,knl

  PUBLIC :: read_pseudopot

  integer,PUBLIC :: pselect = 2

  integer :: Nelement
  integer :: unit_ps,ielm

CONTAINS


  SUBROUTINE read_ppname_pseudopot
    implicit none
    integer :: i
    call IOTools_readIntegerString( "PP", ippform(1), file_ps(1) )
    do i=2,Nelement
       call IOTools_readIntegerString( "PP", ippform(i), file_ps(i), norewind=.true. )
    end do
  END SUBROUTINE read_ppname_pseudopot


  SUBROUTINE read_param_pseudopot
    implicit none
    call IOTools_readIntegerKeyword( "PSELECT", pselect )
  END SUBROUTINE read_param_pseudopot


  SUBROUTINE read_pseudopot( Nelement_in, rank )

    implicit none
    integer,intent(IN) :: Nelement_in, rank
    real(8),allocatable :: psi_(:,:,:),phi_(:,:,:),bet_(:,:,:)
    real(8),allocatable :: ddi_(:,:,:),qqr_(:,:,:)
    integer :: i,j,io,jo,li,lj
    integer :: Lrefmax,Rrefmax,npqmax,nsmpl

    call write_border( 0, " read_pseudopot(start)" )

    Nelement = Nelement_in

    Nelement_PP = Nelement
    Nelement_   = Nelement

    allocate( ippform(Nelement) ) ; ippform=0
    allocate( file_ps(Nelement) ) ; file_ps=""

    call read_ppname_pseudopot

    call read_param_pseudopot

    if ( any(ippform>100) ) pselect=102

    if ( .not.( pselect==2 .or. pselect==3 .or. pselect==102 ) ) then
       stop "invalid pselect(stop@read_param_pseudopot)"
    end if

    allocate( ps(Nelement) )

    if ( rank == 0 ) then

       max_psgrd=0
       max_psorb=0

       do ielm=1,Nelement

          unit_ps=33+ielm
          open(unit_ps,FILE=file_ps(ielm),STATUS='old')

          select case( ippform(ielm) )
          case( 1 )

             close(unit_ps)
             open(unit_ps,FILE=file_ps(ielm),form='unformatted',STATUS='old')

             call ps_read_TM( unit_ps, ps(ielm) )

             call ps_allocate( ps(ielm)%Mr, ps(ielm)%norb )
             Mr(ielm)                 = ps(ielm)%Mr
             norb(ielm)               = ps(ielm)%norb
             Zps(ielm)                = ps(ielm)%Zps
             anorm(1:norb(ielm),ielm) = ps(ielm)%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps(ielm)%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps(ielm)%Rps(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps(ielm)%NRps(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps(ielm)%lo(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps(ielm)%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps(ielm)%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps(ielm)%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps(ielm)%rad(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps(ielm)%rab(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                  = ps(ielm)%viod(1:Mr(ielm),1:norb(ielm))

          case( 2, 102 )

             call read_PSV( unit_ps, ielm, ps(ielm) )

             call ps_allocate( ps(ielm)%Mr, ps(ielm)%norb )
             Mr(ielm)                 = ps(ielm)%Mr
             norb(ielm)               = ps(ielm)%norb
             Zps(ielm)                = ps(ielm)%Zps
             anorm(1:norb(ielm),ielm) = ps(ielm)%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps(ielm)%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps(ielm)%Rps(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps(ielm)%NRps(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps(ielm)%lo(1:norb(ielm))
             no(1:norb(ielm),ielm)    = ps(ielm)%no(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps(ielm)%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps(ielm)%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps(ielm)%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps(ielm)%rad(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps(ielm)%rab(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                  = ps(ielm)%viod(1:Mr(ielm),1:norb(ielm))
             parloc(1:4,ielm)         = ps(ielm)%parloc(1:4)
             nlf(ielm)                = ps(ielm)%nlf
             nrf(1:norb(ielm),ielm)   = ps(ielm)%nrf(1:norb(ielm))

             if ( pselect == 102 ) then !-----> uspp

                Lrefmax = ps(ielm)%nlf
                Rrefmax = maxval( ps(ielm)%nrf )
                npqmax  = (Lrefmax*Rrefmax*(Lrefmax*Rrefmax+1))/2
                nsmpl   = maxval( ps(ielm)%NRps(:) )
                call allocatePSG( Lrefmax,Rrefmax,npqmax,nsmpl,Nelement )

                npq(ielm) = ps(ielm)%npq

                ddi(1:Rrefmax,1:Rrefmax,1:Lrefmax,ielm) &
                     = ps(ielm)%ddi(1:Rrefmax,1:Rrefmax,1:Lrefmax)
                qqr(1:Rrefmax,1:Rrefmax,1:Lrefmax,ielm) &
                     = ps(ielm)%qqr(1:Rrefmax,1:Rrefmax,1:Lrefmax)
                nl3v(1:npqmax,ielm) = ps(ielm)%nl3v(1:npqmax)
                l3v(1:Lrefmax,1:npqmax,ielm) &
                     = ps(ielm)%l3v(1:Lrefmax,1:npqmax)
                qrL(1:nsmpl,1:Lrefmax,1:npqmax,ielm) = &
                     ps(ielm)%qrL(1:nsmpl,1:Lrefmax,1:npqmax)

             else if (  count( ps(ielm)%anorm /= 0.0d0 ) &
                      < count( ps(ielm)%Dij /= 0.0d0 )  ) then !--> MultiRef

                ps_type=1
                do j=1,ps(ielm)%norb
                   jo=ps(ielm)%no(j)
                   lj=ps(ielm)%lo(j)
                do i=1,ps(ielm)%norb
                   io=ps(ielm)%no(i)
                   li=ps(ielm)%lo(i)
                   if ( li /= lj ) cycle
                   hnml(io,jo,li,ielm) = ps(ielm)%Dij(i,j)
                end do
                end do

             end if

          case( 3 )

             call ps_read_YB( unit_ps, ps(ielm) )

             call ps_allocate( ps(ielm)%Mr, ps(ielm)%norb )
             Mr(ielm)                 = ps(ielm)%Mr
             norb(ielm)               = ps(ielm)%norb
             Zps(ielm)                = ps(ielm)%Zps
             anorm(1:norb(ielm),ielm) = ps(ielm)%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps(ielm)%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps(ielm)%Rps(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps(ielm)%NRps(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps(ielm)%lo(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps(ielm)%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps(ielm)%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps(ielm)%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps(ielm)%rad(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps(ielm)%rab(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                  = ps(ielm)%viod(1:Mr(ielm),1:norb(ielm))

          case( 4 )

             call read_ps_gth( unit_ps, ps(ielm) )

             call ps_allocate( 1, ps(ielm)%norb )
             norb(ielm)               = ps(ielm)%norb
             Zps(ielm)                = ps(ielm)%Zps
             Rps(1:norb(ielm),ielm)   = ps(ielm)%Rps(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps(ielm)%lo(1:norb(ielm))
             no(1:norb(ielm),ielm)    = ps(ielm)%no(1:norb(ielm))
             parloc(1:4,ielm)         = ps(ielm)%parloc(1:4)
             Rcloc(ielm)              = ps(ielm)%Rcloc
             hnl(:,:,ielm)            = ps(ielm)%hnl(:,:)
             knl(:,:,ielm)            = ps(ielm)%knl(:,:)
             hnml(:,:,:,ielm)         = ps(ielm)%hnml(:,:,:)
             knml(:,:,:,ielm)         = ps(ielm)%knml(:,:,:)
             inorm(1:norb(ielm),ielm) = ps(ielm)%inorm(1:norb(ielm))

             if ( any( hnml /= 0.0d0 ) ) ps_type=1

          case( 5 )

             call ps_read_UPF( unit_ps, ps(ielm) )

             call ps_allocate( ps(ielm)%Mr, ps(ielm)%norb )
             Mr(ielm)                 = ps(ielm)%Mr
             norb(ielm)               = ps(ielm)%norb
             Zps(ielm)                = ps(ielm)%Zps
             anorm(1:norb(ielm),ielm) = ps(ielm)%anorm(1:norb(ielm))
             inorm(1:norb(ielm),ielm) = ps(ielm)%inorm(1:norb(ielm))
             Rps(1:norb(ielm),ielm)   = ps(ielm)%Rps(1:norb(ielm))
             NRps(1:norb(ielm),ielm)  = ps(ielm)%NRps(1:norb(ielm))
             lo(1:norb(ielm),ielm)    = ps(ielm)%lo(1:norb(ielm))
             no(1:norb(ielm),ielm)    = ps(ielm)%no(1:norb(ielm))
             vql(1:Mr(ielm),ielm)     = ps(ielm)%vql(1:Mr(ielm))
             cdd(1:Mr(ielm),ielm)     = ps(ielm)%cdd(1:Mr(ielm))
             cdc(1:Mr(ielm),ielm)     = ps(ielm)%cdc(1:Mr(ielm))
             rad(1:Mr(ielm),ielm)     = ps(ielm)%rad(1:Mr(ielm))
             rab(1:Mr(ielm),ielm)     = ps(ielm)%rab(1:Mr(ielm))
             viod(1:Mr(ielm),1:norb(ielm),ielm) &
                                      = ps(ielm)%viod(1:Mr(ielm),1:norb(ielm))

             if ( any( ps(ielm)%Dij /= 0.0d0 ) ) then ! Multireference
                ps_type = 1
                do j=1,norb(ielm)
                   jo=no(j,ielm)
                   lj=lo(j,ielm)
                do i=1,norb(ielm)
                   io=no(i,ielm)
                   li=lo(i,ielm)
                   if ( li /= lj ) cycle
                   hnml(io,jo,li,ielm) = ps(ielm)%Dij(i,j)
                end do
                end do
             end if

          case default

             stop "ippform error"

          end select ! ippform

          close(unit_ps)

       end do ! ielm

       write(*,*) "ps_type = ",ps_type

    end if ! [ rank == 0 ]

! --- bcast pseudopotential data

    call send_pseudopot(rank)
    if ( pselect > 100 ) call sendPSG( rank, Nelement )

    do ielm=1,Nelement
       call ps_send_ps1d( ps(ielm) )
       if ( pselect > 100 ) call psg_send_ps1d( ps(ielm) )
    end do

! ---

!    call chk_pot(1,rank)

    call write_border( 0, " read_pseudopot(end)" )

  END SUBROUTINE read_pseudopot


  SUBROUTINE chk_pot(iflag,rank)
    implicit none
    integer,intent(IN) :: iflag,rank
    integer :: u,ielm,i,j
    if ( rank == 0 ) then
       do ielm=1,Nelement
          u=9+ielm
          rewind u
          do i=1,Mr(ielm)
             if ( iflag == 2 ) then
                write(u,'(1x,5f20.10)') rad(i,ielm),cdd(i,ielm),cdc(i,ielm)
             else
                write(u,'(1x,5f20.10)') &
                     rad(i,ielm),vql(i,ielm),(viod(i,j,ielm),j=1,norb(ielm))
             end if
          end do
          write(*,'(1x,"chk_pot(",i1,"): fort.",i2)') iflag,u
       end do
    end if
    stop "stop@chk_pot"
  END SUBROUTINE chk_pot


  SUBROUTINE send_pseudopot_1(myrank)
    implicit none
    integer,intent(IN) :: myrank
    integer :: m,n,ierr
    include 'mpif.h'
    m=max_psgrd
    n=max_psorb
    call mpi_bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Nelement,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( myrank /= 0 ) then
       call ps_allocate(m,n)
    end if
    call mpi_bcast(Mr    ,Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(norb  ,Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Zps   ,Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(parloc,4*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(anorm ,n*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(inorm ,n*Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Rps   ,n*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(NRps  ,n*Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lo    ,n*Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(no    ,n*Nelement,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vql   ,m*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cdd   ,m*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cdc   ,m*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rad   ,m*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rab   ,m*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(viod  ,m*n*Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Rcloc ,Nelement,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!
    call mpi_bcast(max_ngauss,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( max_ngauss /= 0 ) then
       if ( myrank /= 0 ) then
          allocate( cdd_coef(3,max_ngauss,Nelement) ) ; cdd_coef=0.0d0
       end if
       call mpi_bcast(cdd_coef,size(cdd_coef),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    end if
!
    call mpi_bcast(hnl ,size(hnl) ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(knl ,size(knl) ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hnml,size(hnml),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(knml,size(knml),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ps_type,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  END SUBROUTINE send_pseudopot_1


  SUBROUTINE ps_allocate_1(n_grd,n_orb)
    implicit none
    integer,intent(IN) :: n_grd,n_orb
    integer :: mg,mo
    real(8),allocatable :: vql_tmp(:,:),cdd_tmp(:,:),rad_tmp(:,:)
    real(8),allocatable :: cdc_tmp(:,:),viod_tmp(:,:,:)
    real(8),allocatable :: anorm_tmp(:,:),Rps_tmp(:,:),rab_tmp(:,:)
    integer,allocatable :: inorm_tmp(:,:),lo_tmp(:,:),NRps_tmp(:,:)
    integer,allocatable :: no_tmp(:,:)
    if ( .not.allocated(hnl) ) then
       allocate( hnl(3,0:2,Nelement) ) ; hnl=0.0d0
       allocate( knl(3,1:2,Nelement) ) ; knl=0.0d0
    end if
    if ( .not.allocated(hnml) ) then
       allocate( hnml(3,3,0:2,Nelement) ) ; hnml=0.0d0
       allocate( knml(3,3,1:2,Nelement) ) ; knml=0.0d0
    end if
    if ( max_psgrd==0 .or. max_psorb==0 ) then
       allocate( Mr(Nelement)   ) ; Mr=0
       allocate( norb(Nelement) ) ; norb=0
       allocate( Zps(Nelement)  ) ; Zps=0.d0
       allocate( parloc(4,Nelement)    ) ; parloc=0.d0
       allocate( anorm(n_orb,Nelement) ) ; anorm=0.d0
       allocate( inorm(n_orb,Nelement) ) ; inorm=0
       allocate( Rps(n_orb,Nelement)   ) ; Rps=0.d0
       allocate( NRps(n_orb,Nelement)  ) ; NRps=0
       allocate( lo(n_orb,Nelement)    ) ; lo=0
       allocate( no(n_orb,Nelement)    ) ; no=0
       allocate( vql(n_grd,Nelement)   ) ; vql=0.d0
       allocate( cdd(n_grd,Nelement)   ) ; cdd=0.d0
       allocate( cdc(n_grd,Nelement)   ) ; cdc=0.d0
       allocate( rad(n_grd,Nelement)   ) ; rad=0.d0
       allocate( rab(n_grd,Nelement)   ) ; rab=0.d0
       allocate( viod(n_grd,n_orb,Nelement) ) ; viod=0.d0
       allocate( Rcloc(Nelement) ) ; Rcloc=0.0d0
       max_psgrd=n_grd
       max_psorb=n_orb
       return
    end if
    mg = max( max_psgrd, n_grd )
    mo = max( max_psorb, n_orb )
    if ( max_psgrd < mg ) then
       allocate( vql_tmp(mg,Nelement) ) ; vql_tmp=0.0d0
       allocate( cdd_tmp(mg,Nelement) ) ; cdd_tmp=0.0d0
       allocate( rad_tmp(mg,Nelement) ) ; rad_tmp=0.0d0
       allocate( rab_tmp(mg,Nelement) ) ; rab_tmp=0.0d0
       allocate( cdc_tmp(mg,Nelement) ) ; cdc_tmp=0.0d0
       vql_tmp(1:max_psgrd,1:Nelement) = vql(1:max_psgrd,1:Nelement)
       cdd_tmp(1:max_psgrd,1:Nelement) = cdd(1:max_psgrd,1:Nelement)
       rad_tmp(1:max_psgrd,1:Nelement) = rad(1:max_psgrd,1:Nelement)
       rab_tmp(1:max_psgrd,1:Nelement) = rab(1:max_psgrd,1:Nelement)
       cdc_tmp(1:max_psgrd,1:Nelement) = cdc(1:max_psgrd,1:Nelement)
       deallocate( cdc )
       deallocate( rab )
       deallocate( rad )
       deallocate( cdd )
       deallocate( vql )
       allocate( vql(mg,Nelement) ) ; vql=0.d0
       allocate( cdd(mg,Nelement) ) ; cdd=0.d0
       allocate( rad(mg,Nelement) ) ; rad=0.d0
       allocate( rab(mg,Nelement) ) ; rab=0.d0
       allocate( cdc(mg,Nelement) ) ; cdc=0.d0
       vql(:,:)=vql_tmp(:,:)
       cdd(:,:)=cdd_tmp(:,:)
       rad(:,:)=rad_tmp(:,:)
       rab(:,:)=rab_tmp(:,:)
       cdc(:,:)=cdc_tmp(:,:)
       deallocate( cdc_tmp )
       deallocate( rab_tmp )
       deallocate( rad_tmp )
       deallocate( cdd_tmp )
       deallocate( vql_tmp )
       allocate( viod_tmp(mg,mo,Nelement) )
       viod_tmp(1:max_psgrd,1:max_psorb,1:Nelement) &
            = viod(1:max_psgrd,1:max_psorb,1:Nelement)
       deallocate( viod )
       allocate( viod(mg,mo,Nelement) ) ; viod=0.d0
       viod(:,:,:)=viod_tmp(:,:,:)
       deallocate( viod_tmp )
    end if
    if ( max_psorb < mo ) then
       allocate( anorm_tmp(mo,Nelement) ) ; anorm_tmp=0.0d0
       allocate( inorm_tmp(mo,Nelement) ) ; inorm_tmp=0
       allocate( lo_tmp(mo,Nelement) ) ; lo_tmp=0
       allocate( no_tmp(mo,Nelement) ) ; no_tmp=0
       allocate( Rps_tmp(mo,Nelement) ) ; Rps_tmp=0.0d0
       allocate( NRps_tmp(mo,Nelement) ) ; NRps_tmp=0
       anorm_tmp(1:max_psorb,1:Nelement) = anorm(1:max_psorb,1:Nelement)
       inorm_tmp(1:max_psorb,1:Nelement) = inorm(1:max_psorb,1:Nelement)
       lo_tmp(1:max_psorb,1:Nelement) = lo(1:max_psorb,1:Nelement)
       no_tmp(1:max_psorb,1:Nelement) = no(1:max_psorb,1:Nelement)
       Rps_tmp(1:max_psorb,1:Nelement) = Rps(1:max_psorb,1:Nelement)
       NRps_tmp(1:max_psorb,1:Nelement) = NRps(1:max_psorb,1:Nelement)
       deallocate( NRps )
       deallocate( Rps )
       deallocate( no )
       deallocate( lo )
       deallocate( inorm )
       deallocate( anorm )
       allocate( anorm(mo,Nelement) ) ; anorm=0.d0
       allocate( inorm(mo,Nelement) ) ; inorm=0
       allocate( lo(mo,Nelement)    ) ; lo=0
       allocate( no(mo,Nelement)    ) ; no=0
       allocate( Rps(mo,Nelement)   ) ; Rps=0.d0
       allocate( NRps(mo,Nelement)  ) ; NRps=0
       anorm(:,:) = anorm_tmp(:,:)
       inorm(:,:) = inorm_tmp(:,:)
       lo(:,:) = lo_tmp(:,:)
       no(:,:) = no_tmp(:,:)
       Rps(:,:) = Rps_tmp(:,:)
       NRps(:,:) = NRps_tmp(:,:)
       deallocate( NRps_tmp )
       deallocate( Rps_tmp )
       deallocate( lo_tmp )
       deallocate( inorm_tmp )
       deallocate( anorm_tmp )
       if ( max_psgrd >= mg ) then
          allocate( viod_tmp(mg,mo,Nelement) )
          viod_tmp(1:max_psgrd,1:max_psorb,1:Nelement) &
               = viod(1:max_psgrd,1:max_psorb,1:Nelement)
          deallocate( viod )
          allocate( viod(mg,mo,Nelement) ) ; viod=0.d0
          viod(:,:,:)=viod_tmp(:,:,:)
          deallocate( viod_tmp )
       end if
    end if
    max_psgrd = mg
    max_psorb = mo
  END SUBROUTINE ps_allocate_1


END MODULE pseudopot_module
