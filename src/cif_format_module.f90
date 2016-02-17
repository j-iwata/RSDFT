MODULE cif_format_module

  use lattice_module

  implicit none

  PRIVATE
  PUBLIC :: check_cif_format
  PUBLIC :: read_atom_cif

  character(30),parameter :: keyword(11)=(/ &
       "_cell_length_a"   , "_cell_length_b"  , "_cell_length_c", &
       "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma", &
       "_symmetry_equiv_pos_as_xyz", &
       "_atom_site_label", &
       "_atom_site_fract_x","_atom_site_fract_y","_atom_site_fract_z" /)

CONTAINS


  SUBROUTINE check_cif_format( unit, ierr )
    implicit none
    integer,intent(IN)  :: unit
    integer,intent(OUT) :: ierr
    character(30) :: cbuf
    call write_border( 1, " check_cif_format(start)" )
    ierr = -1
    rewind unit
10  continue
    read(unit,*,END=90) cbuf
    if ( any(cbuf==keyword) ) then
       ierr = 0
    else
       goto 10
    end if
90 continue
    call write_border( 1, " check_cif_format(end)" )
  END SUBROUTINE check_cif_format


  SUBROUTINE read_atom_cif &
       ( rank, unit, aa_obj, aa_atom, ki_atom, md_atom, zn_atom)
    implicit none
    integer,intent(IN) :: rank, unit
    type(lattice),intent(INOUT) :: aa_obj
    real(8),intent(OUT),allocatable :: aa_atom(:,:)
    integer,intent(OUT),allocatable :: ki_atom(:),md_atom(:),zn_atom(:)
    logical :: flag
    integer,parameter :: u1=10,u2=20
    integer :: i,j,m,n,i0,i1,i2,i3,z,nsym,itmp(0:3)
    integer :: nbas,isym,natm
    integer,allocatable :: iatm(:),katm(:)
    character(30) :: cbuf, cbuf2
    character(30),allocatable :: cdummy(:)
    real(8),parameter :: bohr=0.529177d0
    real(8) :: alatl(3),angle(3),deg2rad,rr
    real(8) :: R(3,4),asi(3),rsi(3),Rasi(3)
    real(8),allocatable :: rot(:,:,:),atm(:,:)
    include 'mpif.h'

    call write_border( 1, " read_atom_cif(start)" )

    if ( rank == 0 ) then

       do i=1,3
          call find_key( keyword(i), unit, flag )
          if ( flag ) read(unit,*) cbuf, cbuf2
          n=index(cbuf2,"(")-1
          if ( n == -1 ) n=len_trim(cbuf2)
          if ( n > 0 ) read(cbuf2(1:n),*) alatl(i)
          write(*,'(1x,3f15.10)') alatl(i)
       end do
       do i=1,3
          call find_key( keyword(i+3), unit, flag )
          if ( flag ) read(unit,*) cbuf, angle(i)
          write(*,'(1x,3f15.10)') angle(i)
       end do

    end if

#ifndef _NO_MPI_
    call MPI_BCAST( alatl,3,MPI_REAL8,0,MPI_COMM_WORLD,i )
    call MPI_BCAST( angle,3,MPI_REAL8,0,MPI_COMM_WORLD,i )
#endif

    deg2rad = acos(-1.0d0)/180.0d0

    aa_obj%LatticeVector(:,3) = (/ 0.0d0, 0.0d0, 1.0d0 /)

    aa_obj%LatticeVector(:,1) = (/ sin(angle(2)*deg2rad)**2, 0.0d0, 0.0d0 /)

    aa_obj%LatticeVector(:,2) = 0.0d0
    aa_obj%LatticeVector(1,2) = cos(angle(3)*deg2rad)/aa_obj%LatticeVector(1,1)
    aa_obj%LatticeVector(2,2) = sqrt( sin(angle(1)*deg2rad)**2 &
                                    - aa_obj%LatticeVector(1,2)**2 )

    aa_obj%LatticeVector(:,1) = aa_obj%LatticeVector(:,1)*alatl(1)/bohr
    aa_obj%LatticeVector(:,2) = aa_obj%LatticeVector(:,2)*alatl(2)/bohr
    aa_obj%LatticeVector(:,3) = aa_obj%LatticeVector(:,3)*alatl(3)/bohr

    aa_obj%LatticeConstant = alatl(1)/bohr

    if ( rank == 0 ) then
       do i=1,3
          write(*,'(1x,3f15.10)') aa_obj%LatticeVector(:,i)
       end do
    end if

    if ( rank == 0 ) then

       call find_key( keyword(7), unit, flag )

       nsym=0

       if ( flag ) then

          read(unit,*) cbuf

          open(u1,file="cif_sym.dat")

          i=0
20        read(unit,*) cbuf
          if ( cbuf /= "loop_" ) then
             i=i+1
             call chr_to_matrix( cbuf, R )
             write(u1,'(1x,3f20.15,2x,f20.15)') (R(1,j),j=1,4)
             write(u1,'(1x,3f20.15,2x,f20.15)') (R(2,j),j=1,4)
             write(u1,'(1x,3f20.15,2x,f20.15)') (R(3,j),j=1,4)
             goto 20
          endif
          nsym=i

          close(u1)

       end if

       if ( nsym == 0 ) then
          nsym=1
          allocate( rot(3,4,nsym) ) ; rot=0.0d0
          rot(1,1,1)=1.0d0
          rot(2,2,1)=1.0d0
          rot(3,3,1)=1.0d0
       else if ( nsym > 0 ) then
          allocate( rot(3,4,nsym) ) ; rot=0.0d0
          open(u1,file="cif_sym.dat",status="old")
          do i=1,nsym
             read(u1,*) (rot(1,j,i),j=1,4)
             read(u1,*) (rot(2,j,i),j=1,4)
             read(u1,*) (rot(3,j,i),j=1,4)
          end do
          close(u1)
       end if

       itmp=0
       i=0
30     read(unit,*,END=92) cbuf
       i=i+1
       if ( cbuf == keyword( 8) ) itmp(0)=i
       if ( cbuf == keyword( 9) ) itmp(1)=i
       if ( cbuf == keyword(10) ) itmp(2)=i
       if ( cbuf == keyword(11) ) itmp(3)=i
       if ( cbuf(1:5) /= "_atom" ) then
          backspace(unit)
       else
          goto 30
       end if
92     continue

       n=maxval(itmp)
       allocate( cdummy(n) ) ; cdummy=""

       open(u2,file="cif_bas.dat")

       nbas=0

40     read(unit,*,END=94) cdummy(1:n)

       call get_atomic_number( cdummy(itmp(0)), z )
       if ( z > 0 ) then
          do i=1,3
             cbuf=cdummy(itmp(i))
             m=index(cbuf,"(")-1
             if ( m == -1 ) m=len_trim(cbuf)
             read(cbuf(1:m),*) asi(i)
          end do
          nbas=nbas+1
          write(*,'(1x,i3,2x,a2,2x,i3,2x,3f15.10)') &
               nbas,cdummy(itmp(0)),z,(asi(i),i=1,3)
          do isym=1,nsym
             Rasi(:) = matmul( rot(:,1:3,isym), asi(:) )
             Rasi(:) = Rasi(:) + rot(:,4,isym)
             write(u2,'(1x,a3,2x,i3,2x,3f20.15)') cdummy(itmp(0)),z,Rasi(:)
          end do
          goto 40
       end if  
94     continue

       close(u2)

       deallocate( cdummy )
       deallocate( rot )

       natm = nbas*nsym

       allocate( atm(3,natm) ) ; atm=0.0d0
       allocate( katm(natm)  ) ; katm=0
       allocate( iatm(118)   ) ; iatm=0

       open(u2,file="cif_bas.dat",status="old")

       n=0
       loop_i : do i=1,natm
          read(u2,*) cbuf,z,asi(:)
          call shift_aa_coordinates_atom( asi )
          do j=1,n
             rr=sum( (asi-atm(:,j))**2 )
             if ( rr < 1.d-8 ) cycle loop_i
          end do
          if ( iatm(z) == 0 ) iatm(z)=maxval(iatm)+1
          n=n+1
          atm(:,n)=asi(:)
          katm(n)=iatm(z)
          write(*,'(1x,2(i3,2x),3f10.5)') n,katm(n),(asi(j),j=1,3)
       end do loop_i

       close(u2)


       allocate( aa_atom(3,n) ) ; aa_atom=0.0d0
       allocate( ki_atom(n)   ) ; ki_atom=0
       allocate( md_atom(n)   ) ; md_atom=1

       aa_atom(:,:) = atm(:,1:n)
       ki_atom(:)   = katm(1:n)

       m=maxval( iatm )
       allocate( zn_atom(m) ) ; zn_atom=1

       do z=1,size(iatm)
          i=iatm(z)
          if ( i > 0 ) zn_atom(i)=z
       end do

       deallocate( iatm )
       deallocate( katm )
       deallocate( atm  )

    end if

#ifndef _NO_MPI_
    call MPI_BCAST( n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i )
#endif
    if ( .not.allocated(aa_atom) ) then
       allocate( aa_atom(3,n) ) ; aa_atom=0.0d0
       allocate( ki_atom(n)   ) ; ki_atom=0
       allocate( md_atom(n)   ) ; md_atom=0
       allocate( zn_atom(n)   ) ; zn_atom=0
    end if
#ifndef _NO_MPI_
    call MPI_BCAST(aa_atom,3*n,MPI_REAL8,0,MPI_COMM_WORLD,i)
    call MPI_BCAST(ki_atom,n,MPI_INTEGER,0,MPI_COMM_WORLD,i)
    call MPI_BCAST(md_atom,n,MPI_INTEGER,0,MPI_COMM_WORLD,i)
    call MPI_BCAST(zn_atom,n,MPI_INTEGER,0,MPI_COMM_WORLD,i)
#endif

    call write_border( 1, " read_atom_cif(end)" )

  END SUBROUTINE read_atom_cif

  SUBROUTINE find_key( key, unit, flag )
    implicit none
    character(*),intent(IN) :: key
    integer,intent(IN) :: unit
    logical,intent(OUT) :: flag
    character(30) :: cbuf
    flag=.false.
    rewind unit
10  read(unit,*,END=90) cbuf
    if ( cbuf == key ) then
       backspace(unit)
       flag=.true.
       return
    else
       goto 10
    end if
90  continue
    call stop_program_f( "stop@find_key(cif_format_module)" )
  END SUBROUTINE find_key

  SUBROUTINE chr_to_matrix( cbuf, R )
    implicit none
    character(*),intent(IN) :: cbuf
    real(8),intent(OUT) :: R(3,4)
    character(10) :: str(3)
    integer :: i,j
    real(8) :: f
    R=0.0d0
    write(*,*) cbuf
    j=1
    f=1.0d0
    do i=1,len_trim(cbuf)
       if ( cbuf(i:i) == "-" ) then
          f=-1.0d0
       end if
       if ( cbuf(i:i) == "x" ) then
          R(j,1) = f
          f=1.0d0
       end if
       if ( cbuf(i:i) == "y" ) then
          R(j,2) = f
          f=1.0d0
       end if
       if ( cbuf(i:i) == "z" ) then
          R(j,3) = f
          f=1.0d0
       end if
       if ( 49 <= iachar(cbuf(i:i)) .and. iachar(cbuf(i:i)) <= 57 ) then
          if ( R(j,4) == 0.0d0 ) then
             R(j,4) = f*(iachar(cbuf(i:i))-48)
          else
             R(j,4) = R(j,4)/dble( iachar(cbuf(i:i))-48 )
          end if
          f=1.0d0
       end if
       if ( cbuf(i:i) == "," ) then
          j=j+1
          f=1.0d0
       end if
    end do
    do i=1,3
       write(*,'(1x,3f15.10,3x,f15.10)') (R(i,j),j=1,4)
    end do
  END SUBROUTINE chr_to_matrix


  SUBROUTINE shift_aa_coordinates_atom( aa_atom )
    implicit none
    real(8),intent(INOUT) :: aa_atom(:)
    integer :: j
    do j=1,3
10     if ( aa_atom(j) < 0.0d0 ) then
          aa_atom(j)=aa_atom(j)+1.0d0
          goto 10
       else if ( aa_atom(j) >= 1.0d0 ) then
          aa_atom(j)=aa_atom(j)-1.0d0
          goto 10
       else
          cycle
       end if
    end do
  END SUBROUTINE shift_aa_coordinates_atom


END MODULE cif_format_module
