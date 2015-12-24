! assuming the unit of POSCAR is in angstrome,
! and atomic coordinates are given as Cartesian (xyz) coordinates

PROGRAM poscar2aa

  implicit none

  integer,parameter :: u1 = 5,  u2 = 970, u3 = 6
  real(8),parameter :: bohr = 0.529177d0
  real(8),allocatable :: asi(:,:),rsi(:,:)
  real(8) :: ax,aa(3,3),bb(3,3)
  character(30) :: cbuf, cbuf2
  integer :: i,j,k,l2,ne,ichk,natm_tot,natm(10)

  read(u1,*)
  read(u1,*) ax
  read(u1,*) aa(1:3,1)
  read(u1,*) aa(1:3,2)
  read(u1,*) aa(1:3,3)

  read(u1,*)

  read(u1,'(a)') cbuf

  cbuf2 = adjustl(cbuf)
  l2    = len_trim(cbuf2)

  ne=0
  ichk=0
  do i=1,l2
     write(*,*) i,cbuf2(i:i)
     if ( cbuf2(i:i) == " " ) then
        ichk=0
     else
        if ( ichk == 1 ) cycle
        ichk=1
        ne=ne+1
     end if
  end do

  write(*,*) "ne =",ne

  backspace(u1)

  read(u1,*) natm(1:ne)
  read(u1,*)

  natm_tot = sum( natm(1:ne) )

  allocate( rsi(3,natm_tot) ) ; rsi=0.0d0
  allocate( asi(3,natm_tot) ) ; asi=0.0d0

  i=0
  do k=1,ne
     do j=1,natm(k)
        i=i+1
        read(u1,*) rsi(1:3,i)
     end do
  end do

! convert

  ax = ax/bohr

  aa(:,:) = ax*aa(:,:)

  rsi(:,:) = ax*rsi(:,:)

  call calc_aainv( aa, bb )

  asi(:,:) = matmul( bb(:,:), rsi(:,:) )

! output

  write(u3,'("AX",f20.15)') ax
  write(u3,'("A1",3f20.15)') aa(1:3,1)/ax
  write(u3,'("A2",3f20.15)') aa(1:3,2)/ax
  write(u3,'("A3",3f20.15)') aa(1:3,3)/ax
  write(u3,'("AA")')

  write(u2,*) ne, natm_tot, "   /"

  i=0
  do k=1,ne
     do j=1,natm(k)
        i=i+1
        write(u2,'(1x,i4,3f20.15,i4," /")') k, asi(1:3,i), 1
     end do
  end do

CONTAINS

  SUBROUTINE calc_aainv( aa, bb )
    implicit none
    real(8),intent(IN)  :: aa(3,3)
    real(8),intent(OUT) :: bb(3,3)
    real(8) :: va
    va = aa(1,1)*aa(2,2)*aa(3,3)+aa(1,2)*aa(2,3)*aa(3,1) &
        +aa(1,3)*aa(2,1)*aa(3,2)-aa(1,3)*aa(2,2)*aa(3,1) &
        -aa(1,2)*aa(2,1)*aa(3,3)-aa(1,1)*aa(2,3)*aa(3,2)
    bb(1,1) = aa(2,2)*aa(3,3) - aa(3,2)*aa(2,3)
    bb(2,1) = aa(3,2)*aa(1,3) - aa(1,2)*aa(3,3)
    bb(3,1) = aa(1,2)*aa(2,3) - aa(2,2)*aa(1,3)
    bb(1,2) = aa(2,3)*aa(3,1) - aa(3,3)*aa(2,1)
    bb(2,2) = aa(3,3)*aa(1,1) - aa(1,3)*aa(3,1)
    bb(3,2) = aa(1,3)*aa(2,1) - aa(2,3)*aa(1,1)
    bb(1,3) = aa(2,1)*aa(3,2) - aa(3,1)*aa(2,2)
    bb(2,3) = aa(3,1)*aa(1,2) - aa(1,1)*aa(3,2)
    bb(3,3) = aa(1,1)*aa(2,2) - aa(2,1)*aa(1,2)
    bb(:,:) = bb(:,:)/va
    bb(:,:) = transpose( bb )
  END SUBROUTINE calc_aainv

END PROGRAM poscar2aa
