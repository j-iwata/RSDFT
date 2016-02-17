!--------------------------------------------
! Unified Pseudopotential Format (UPF)
! This potential is adopted in QUANTUM ESPRESSO 
! The unit of energy is in Rydberg, and converted to Hartree in this routine
!--------------------------------------------
MODULE ps_read_UPF_module

  use var_ps_member, only: ps1d, ps_allocate_ps1d

  implicit none

  PRIVATE
  PUBLIC :: ps_read_UPF

CONTAINS

  SUBROUTINE ps_read_UPF( g, psp )
    implicit none
    integer,intent(IN) :: g
    type(ps1d),intent(INOUT) :: psp
    integer,parameter :: max_loop = 100000
    integer :: loop
    character(30) :: cbuf

    do loop=1,max_loop

       read(g,'(a)',END=10) cbuf
       write(*,*) cbuf

       if ( cbuf(1:21) == '<UPF version="2.0.1">' ) then
          rewind g
          call ps_read_upf_ver201( g, psp )
          return
       else if ( cbuf(1:9) == "<PP_INFO>" ) then
          rewind g
          call ps_read_upf_verorg( g, psp )
          return
       end if

    end do ! loop

10  stop "Format is invalid (stop@ps_read_upf)"

  END SUBROUTINE ps_read_UPF


  SUBROUTINE ps_read_upf_verorg( g, psp )
    implicit none
    integer,intent(IN) :: g
    type(ps1d),intent(INOUT) :: psp
    integer,parameter :: max_loop = 100000, max_array_size = 8
    integer :: i,j,n,i0,i1,loop,nrr,norb,nrc
    integer,allocatable :: lo(:),NRps(:),inorm(:)
    real(8) :: tmp,Zps
    real(8),allocatable :: rr(:),rx(:),vql(:),cdc(:),cdd(:)
    real(8),allocatable :: viod(:,:),anorm(:)
    character(30) :: cbuf

    write(*,'(a40," ps_read_upf_verorg")') repeat("-",40)

! Read

    nrr=0

    do loop=1,max_loop

       read(g,*,END=10) cbuf

       if ( cbuf(1:11) == "<PP_HEADER>" ) then

          do i=1,5
             read(g,*)
          end do
          read(g,*) Zps

          do i=1,3
             read(g,*)
          end do
          read(g,*) nrr

       else if ( cbuf(1:11) == "<PP_HEADER " ) then

          do i=1,max_loop
             read(g,*) cbuf
             if ( cbuf(1:9) == "z_valence" ) then
                stop "test"
             end if
          end do

       end if

       if ( nrr > 0 ) then
          if ( .not.allocated(rr) ) then
             allocate( rr(1:nrr)  ) ; rr =0.0d0
             allocate( rx(1:nrr)  ) ; rx =0.0d0
             allocate( vql(1:nrr) ) ; vql=0.0d0
             allocate( cdc(1:nrr) ) ; cdc=0.0d0
             allocate( cdd(1:nrr) ) ; cdd=0.0d0
          end if
       end if

       if ( cbuf(1:9) == "<PP_MESH>" ) then

          do i=1,max_loop
             read(g,*) cbuf
             if ( cbuf(1:6) == "<PP_R>" ) exit
          end do

          read(g,*) rr(1:nrr)

          do i=1,max_loop
             read(g,*) cbuf
             if ( cbuf(1:8) == "<PP_RAB>" ) exit
          end do ! i

          read(g,*) rx(1:nrr)

       end if ! </PP_MESH>

       if ( cbuf(1:9) == "<PP_NLCC>" ) then

          read(g,*) cdc(1:nrr)

       end if

       if ( cbuf(1:10) == "<PP_LOCAL>" ) then

          read(g,*) vql(1:nrr)

       end if ! </PP_LOCAL>

       if ( nrr > 0 ) then
          if ( .not.allocated(lo) ) then
             allocate( lo(max_array_size)       ) ; lo=0
             allocate( NRps(max_array_size)     ) ; NRps=0
             allocate( inorm(max_array_size)    ) ; inorm=0
             allocate( anorm(max_array_size)    ) ; anorm=0.0d0
             allocate( viod(nrr,max_array_size) ) ; viod=0.0d0
          end if
       end if

       if ( cbuf(1:13) == "<PP_NONLOCAL>" ) then

          norb=0
          do i=1,max_loop
             read(g,*) cbuf
             if ( cbuf(1:9) == "<PP_BETA>" ) then
                read(g,*) j, lo(j)
                read(g,*) nrc ; if ( nrc > nrr) stop"stop@ps_read_upf(1)"
                read(g,*) viod(1:nrc,j)
                read(g,*) cbuf
                NRps(j)=nrc
                norb=max( j, norb )
             else if ( cbuf(1:8) == "<PP_DIJ>" ) then
                read(g,*) n
                do j=1,n
                   read(g,*) i0,i1,anorm(j)
                end do
                do j=1,n
                   inorm(j)=nint( sign(1.0d0,anorm(j)) )
                   anorm(j)=abs( anorm(j) )
                end do
                read(g,*) cbuf
             else if ( cbuf(1:8) == "<" ) then
                write(*,*) cbuf,"exit"
                exit
             end if
          end do

       end if

       if ( cbuf(1:10) == "<PP_PSWFC>" ) then
       end if

       if ( cbuf(1:12) == "<PP_RHOATOM>" ) then

          read(g,*) cdd(1:nrr)

       end if

    end do ! loop
10 continue

    psp%Mr = nrr + 1
    psp%norb = norb
    call ps_allocate_ps1d( psp )

    psp%Zps = Zps

    do i=nrr+1,2,-1
       psp%rad(i) =  rr(i-1)
       psp%rab(i) =  rx(i-1)
       psp%cdd(i) = cdd(i-1)
       psp%cdc(i) = cdc(i-1)
       psp%vql(i) = vql(i-1)
    end do
    psp%rad(1) = 0.0d0
    psp%rab(1) = 0.0d0
    psp%cdd(1) = 0.0d0
    psp%cdc(1) = 0.0d0
    psp%vql(1) = psp%vql(2)

    psp%vql(:) = 0.5d0*psp%vql(:)

    do j=1,norb
       psp%lo(j)    = lo(j)
       psp%inorm(j) = inorm(j)
       psp%anorm(j) = anorm(j)
    end do

    do j=1,norb
       do i=NRps(j)+1,2,-1
          psp%viod(i,j) = viod(i-1,j)
       end do
       psp%viod(1,j) = 0.0d0
    end do

    do j=1,norb
       do i=nrr,1,-1
          if ( abs(psp%viod(i,j)) >= 1.d-13 ) then
             psp%Rps(j)  = psp%rad(i-1)
             psp%NRps(j) = i-1
             exit
          end if
       end do
    end do

    do j=1,norb
       tmp = sqrt( abs( psp%anorm(j) ) )
       psp%viod(:,j) = sqrt(0.5d0)*psp%viod(:,j)*tmp
    end do

    deallocate( viod, anorm, inorm, NRps, lo )
    deallocate( cdd, cdc, vql, rx, rr )

    write(*,*) "*** Unified Pseudopotenetial Format (UPF) ***"
    write(*,*) "Znuc=",psp%Zps
    write(*,*) "# of radial mesh points =",psp%Mr
    write(*,*) "# of orbitals =",psp%norb
    write(*,*) "angular momentum =",psp%lo(1:psp%norb)
    write(*,*) "cut off radius =",psp%Rps(1:psp%norb)
    write(*,*) "# of grid points within cut off radius",psp%NRps(1:psp%norb)
    write(*,*) "uVu integral (anorm) ="
    write(*,'(1x,8f10.5)') ( psp%inorm(i)*psp%anorm(i),i=1,psp%norb )
    write(*,*) "sum(rhov)=",sum(psp%cdd*psp%rab)
    write(*,*) "sum(rhoc)=",sum(psp%cdc*psp%rab*(psp%rad)**2)*4*acos(-1.d0)

    write(*,'(a40," ps_read_upf_verorg(end)")') repeat("-",40)

    return
  END SUBROUTINE ps_read_upf_verorg


  SUBROUTINE ps_read_upf_ver201( g, psp )
    implicit none
    integer,intent(IN) :: g
    type(ps1d),intent(INOUT) :: psp
    integer,parameter :: max_loop=1000000, max_array_size = 8
    integer :: loop,i,j,k,l,ir,ic,nr
    integer,allocatable :: lo(:),no(:)
    character(100) :: cbuf, ckey
    integer :: norb,nrr,nsize,ltmp,columns
    real(8) :: tmp,Zps
    real(8),allocatable :: rr(:),rx(:),vql(:),cdc(:),cdd(:)
    real(8),allocatable :: viod(:,:),anorm(:),Dij(:,:),work(:)

    write(*,'(a40," ps_read_upf_ver201")') repeat("-",40)

    nrr=0
    norb=0

    do loop=1,max_loop

       read(g,'(a)',END=10) cbuf
       ckey = adjustl( cbuf )

       if ( ckey(1:6) == "</UPF>" ) then
          write(*,*) ckey(1:6)
          exit
       end if

       if ( ckey(1:11) == "<PP_HEADER " ) then

          write(*,*) ckey(1:11)

          backspace(g)

          do i=1,max_loop

             read(g,'(a)') cbuf

             j = index( cbuf, "z_valence=" )
             if ( j /= 0 ) then
                j = j + len_trim( "z_valence=" )
                k = index( cbuf(j+1:), '"' )
                ckey=cbuf(j+1:j+k-1)
                read(ckey,*) Zps
             end if

             j = index( cbuf, "mesh_size=" )
             if ( j /= 0 ) then
                j = j + len_trim( "mesh_size=" )
                k = index( cbuf(j+1:), '"' )
                ckey=cbuf(j+1:j+k-1)
                read(ckey,*) nrr
             end if

             j = index( cbuf, "number_of_proj=" )
             if ( j /= 0 ) then
                j = j + len_trim( "number_of_proj=" )
                k = index( cbuf(j+1:), '"' )
                ckey=cbuf(j+1:j+k-1)
                read(ckey,*) norb
             end if

             j = index( cbuf, "/>" )
             if ( j /= 0 ) exit

          end do ! i

       end if ! </PP_HEADER>

       if ( nrr > 0 ) then
          if ( .not.allocated(rr) ) then
             allocate( rr(nrr)  ) ; rr =0.0d0
             allocate( rx(nrr)  ) ; rx =0.0d0
             allocate( cdc(nrr) ) ; cdc=0.0d0
             allocate( vql(nrr) ) ; vql=0.0d0
             allocate( cdd(nrr) ) ; cdd=0.0d0
          end if
       end if

       if ( ckey(1:9) == "<PP_MESH " ) then

          write(*,*) ckey(1:9)

          do i=1,max_loop
             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )
             if ( ckey(1:6) == "<PP_R " ) exit
          end do
          read(g,*) rr(1:nrr)

          do i=1,max_loop
             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )
             if ( ckey(1:8) == "<PP_RAB " ) exit
          end do
          read(g,*) rx(1:nrr)

       end if ! </PP_MESH>

       if ( ckey(1:9) == "<PP_NLCC " ) then

          write(*,*) ckey(1:9)

          read(g,*) cdc(1:nrr)

       end if

       if ( ckey(1:10) == "<PP_LOCAL " ) then

          write(*,*) ckey(1:10)

          read(g,*) vql(1:nrr)

       end if ! </PP_LOCAL>

       if ( nrr > 0 ) then
          if ( .not.allocated(lo) ) then
             allocate( lo(max_array_size)       ) ; lo=0
             allocate( no(0:max_array_size)     ) ; no=0
             allocate( viod(nrr,max_array_size) ) ; viod=0.0d0
          end if
       end if

       if ( ckey(1:13) == "<PP_NONLOCAL>" ) then

          write(*,*) ckey(1:13)

          j=0
          do i=1,max_loop

             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )

             if ( ckey(1:9) == "<PP_BETA." ) then

                do k=1,max_loop

                   call get_num_from_string( cbuf, "angular_momentum=", ltmp )

                   l=index( cbuf, ">" )
                   if ( l /= 0 ) then
                      j=j+1
                      read(g,*) viod(1:nrr,j)
                      lo(j)=ltmp
                      no(lo(j)) = no(lo(j)) + 1
                      exit
                   end if

                   read(g,'(a)') cbuf

                end do ! k

                if ( j == norb ) exit

             end if

          end do ! i

          if ( norb > 0 ) then
             if ( .not.allocated(Dij) ) then
                allocate( Dij(norb,norb) ) ; Dij=0.0d0
                allocate( anorm(norb)    ) ; anorm=0.0d0
             end if
          end if

          do i=1,max_loop

             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )

             if ( ckey(1:8) == "<PP_DIJ " ) then
                write(*,*) ckey(1:8)
                call get_num_from_string( cbuf, "size=", nsize )
                call get_num_from_string( cbuf, "columns=", columns )
                write(*,*) "nsize,columns=",nsize,columns
! ---
                allocate( work(nsize) ) ; work=0.0d0
                nr=nsize/columns
                do ir=1,nr
                   read(g,*) work((ir-1)*columns+1:ir*columns)
                end do
                if ( nr*columns < nsize ) read(g,*) work(nr*columns+1:nsize)
                do ic=1,norb
                   do ir=1,norb
                      Dij(ir,ic) = work(ir+(ic-1)*norb)
                   end do
                end do
                deallocate( work )
! ---
                j=count( Dij /= 0.0d0 )
                k=0 !------> check single or multi reference
                do l=1,norb
                   if ( Dij(l,l) /= 0.0d0 ) k=k+1
                end do
                if ( j == k ) then
                   do l=1,norb
                      anorm(l) = Dij(l,l)
                   end do
                   Dij(:,:)=0.0d0
                end if !------> check single or multi reference (end)
                exit
             end if

          end do ! i

       end if ! </PP_NONLOCAL>

       if ( ckey(1:12) == "<PP_RHOATOM " ) then

          write(*,*) ckey(1:12)

          read(g,*) cdd(1:nrr)
          exit

       end if

    end do ! loop
10 continue

    psp%Mr = nrr + 1
    psp%norb = norb
    call ps_allocate_ps1d( psp )

    psp%Zps = Zps

    do i=nrr+1,2,-1
       psp%rad(i) =  rr(i-1)
       psp%rab(i) =  rx(i-1)
       psp%cdd(i) = cdd(i-1)
       psp%cdc(i) = cdc(i-1)
       psp%vql(i) = vql(i-1)
    end do
    psp%rad(1) = 0.0d0
    psp%rab(1) = 0.0d0
    psp%cdd(1) = 0.0d0
    psp%cdc(1) = 0.0d0
    psp%vql(1) = psp%vql(2)

    psp%vql(:) = 0.5d0*psp%vql(:)

    no(:)=0
    do j=1,norb
       psp%lo(j) = lo(j)
       no( lo(j) ) = no( lo(j) ) + 1
       psp%no(j) = no( lo(j) )
    end do

    do j=1,norb
       do i=nrr+1,2,-1
          psp%viod(i,j) = viod(i-1,j)
       end do
       psp%viod(1,j) = 0.0d0
    end do

    do j=1,norb
       do i=nrr,1,-1
          if ( abs(psp%viod(i,j)) >=1.d-13 ) then
             psp%Rps(j) = psp%rad(i-1)
             psp%NRps(j)= i-1
             exit
          end if
       end do
    end do

    if ( any( anorm /= 0.0d0 ) ) then !------> single reference
       do j=1,norb
          psp%anorm(j)=anorm(j)
          psp%inorm(j)=1
          if ( psp%anorm(j) < 0.0d0 ) psp%inorm(j)=-1
          tmp = sqrt( abs( psp%anorm(j) ) )
          psp%viod(:,j) = sqrt(0.5d0)*psp%viod(:,j)*tmp
       end do
    else !------> multi reference
       psp%Dij(1:norb,1:norb) = Dij(1:norb,1:norb)
       do j=1,norb
          psp%viod(:,j) = sqrt(0.5d0)*psp%viod(:,j)
       end do
    end if

    deallocate( Dij, anorm )
    deallocate( viod, no, lo )
    deallocate( cdd, vql, cdc, rx, rr )

    write(*,*) "*** Unified Pseudopotenetial Format (UPF) ***"
    write(*,*) "Znuc=",psp%Zps
    write(*,*) "# of radial mesh points =",psp%Mr
    write(*,*) "# of orbitals =",psp%norb
    write(*,*) "angular momentum =",psp%lo(1:psp%norb)
    write(*,*) "cut off radius =",psp%Rps(1:psp%norb)
    write(*,*) "# of grid points within cut off radius",psp%NRps(1:psp%norb)
    write(*,*) "uVu integral (anorm) ="
    write(*,'(1x,8f10.5)') ( psp%inorm(i)*psp%anorm(i),i=1,psp%norb )
    write(*,*) "Dij ="
    write(*,'(1x,9f10.5)') (( psp%Dij(i,j),i=1,psp%norb ),j=1,psp%norb)
    write(*,*) "sum(rhov)=",sum(psp%cdd*psp%rab)
    write(*,*) "sum(rhoc)=",sum(psp%cdc*psp%rab*(psp%rad)**2)*4*acos(-1.d0)

    write(*,'(a40," ps_read_upf_ver201(end)")') repeat("-",40)

    return
  END SUBROUTINE ps_read_upf_ver201


  SUBROUTINE get_num_from_string( buf, str, n )
    implicit none
    character(*),intent(IN) :: buf, str
    integer,intent(OUT) :: n
    integer :: i, l, len_str
    character(8) :: str_out
    l = index( buf, str )
    if ( l /= 0 ) then
       len_str = len_trim( adjustl(str) )
       i=l+len_str+1
       str_out = buf(i:i)
       read(str_out,*) n
    end if
  END SUBROUTINE get_num_from_string


END MODULE ps_read_UPF_module
