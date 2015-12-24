SUBROUTINE write_border( n, indx )

  implicit none
  integer,intent(IN) :: n
  character(*),intent(IN) :: indx
  character(80) :: axx
  logical :: disp
  integer :: m=80
  integer,save :: u0=6, u1=60

  write(axx,'(i2)') m-len(indx)
  axx=adjustl(axx)
  axx="(a"//axx(1:len_trim(axx))//",a)"

  call check_disp_switch( disp, 0 )
  if ( disp ) then
     if ( n == 0 ) then
        write(u0,axx) repeat("-",m),indx
     else if ( n == 1 ) then
!        open(u1,file="RSDFT_LOG",position="append")
!        write(u1,axx) repeat("-",m),indx
!        close(u1)
     else
        write(u0,axx) repeat("-",m),indx
     end if
     if ( u0 /= u1 ) then
        open(u1,file="RSDFT_LOG",position="append")
        write(u1,axx) repeat("-",m),indx
        close(u1)
     end if
  end if

END SUBROUTINE write_border


SUBROUTINE check_disp_switch( disp_switch, i )
  implicit none
  logical,intent(INOUT) :: disp_switch
  integer,intent(IN)    :: i
  logical,save :: disp_switch_save=.false.
  if ( i == 0 ) then
     disp_switch = disp_switch_save
  else if ( i == 1 ) then
     disp_switch_save = disp_switch
  end if
END SUBROUTINE check_disp_switch


SUBROUTINE check_disp_length( level, i )
  implicit none
  integer,intent(INOUT) :: level
  integer,intent(IN)    :: i
  logical,save :: info_level=0
  if ( i == 0 ) then
     level = info_level
  else if ( i == 1 ) then
     info_level = level
  end if
END SUBROUTINE check_disp_length


SUBROUTINE write_string( string )
  implicit none
  character(*),intent(IN) :: string
  integer,save :: u0=6, u1=60
  logical :: disp
  call check_disp_switch( disp, 0 )
  if ( disp ) write(u0,'(a)') string
END SUBROUTINE write_string


SUBROUTINE write_int_and_real( format_string, m, fi, n, fr )
  implicit none
  character(*),intent(IN) :: format_string
  integer,intent(IN) :: m,n,fi(m)
  real(8),intent(IN) :: fr(n)
  integer,save :: u0=6, u1=60
  integer :: i,j
  logical :: disp
  call check_disp_switch( disp, 0 )
  if ( disp ) write(u0,format_string) ( fi(i), i=1,m ), ( fr(j), j=1,n )
END SUBROUTINE write_int_and_real
