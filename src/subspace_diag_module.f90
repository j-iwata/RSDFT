MODULE subspace_diag_module

  use parallel_module, only: np_band, ir_band, id_band, myrank &
                            ,disp_switch_parallel
  use subspace_diag_variables
  use subspace_diag_la_module
  use subspace_diag_sl_module

  implicit none

  PRIVATE
  PUBLIC :: subspace_diag, init_subspace_diag

CONTAINS


  SUBROUTINE subspace_diag(k,s)
    implicit none
    integer,intent(IN) :: k,s
#ifdef _LAPACK_
    call subspace_diag_la(k,s)
#else
    call subspace_diag_sl(k,s,disp_switch_parallel)
#endif
  END SUBROUTINE subspace_diag


  SUBROUTINE init_subspace_diag( MB_in )
    implicit none
    integer,intent(IN) :: MB_in
    integer :: i,j,mm,ms,me,nme,ne,nn,je,MB

    call write_border( 80, " init_subspace_diag(start)" )

    MB_diag = MB_in

    MB  = MB_diag
    nme = (MB*MB+MB)/2

    call parameter_check(nme,MB)

    if ( .not.allocated(mat_block) ) then
       allocate( mat_block(0:np_band-1,0:4) )
    end if
    mat_block(:,:) = 0

    do i=0,np_band-1
       me=id_band(i)+ir_band(i)
       mm=ir_band(i)
       mat_block(i,0)=(mm*(mm+1))/2
       mat_block(i,1)=MB-me
       mat_block(i,2)=mm
       mat_block(i,3)=mat_block(i,0)+mat_block(i,1)*mat_block(i,2)
    end do

    if ( sum(mat_block(:,3)) /= nme ) then
       write(*,*) sum(mat_block(:,3)),myrank,nme
       stop "stop@init_subspace_diag(1)"
    end if

    if ( np_band>1 ) then

       je = int( (np_band+1)*0.5 )-1
       do j=0,je
          do i=np_band-1,j+1,-1
             mm=ir_band(i)
             if( ((np_band-1)*np_band*0.5/np_band+1)*mm < mat_block(j,1) )then
                mat_block(j,1)=mat_block(j,1)-mm
                mat_block(i,1)=mat_block(i,1)+mm
             end if
          end do
       end do
       mat_block(:,3)=mat_block(:,0)+mat_block(:,1)*mat_block(:,2)

       if ( sum(mat_block(:,3))/=nme ) then
          write(*,*) sum(mat_block(:,3)),myrank,nme
          stop "stop@init_subspace_diag(2)"
       end if

    end if

    do i=0,np_band-1
       mat_block(i,4)=sum( mat_block(0:i,3) )-mat_block(i,3)
    end do

    if ( disp_switch_parallel ) then
       write(*,'(1x,6a10)') "rank_b","tri","m","n","nme","idis"
       do i=0,np_band-1
          write(*,'(1x,6i10)') i,mat_block(i,0:4)
       end do
    end if

    call write_border( 80, " init_subspace_diag(end)" )

  END SUBROUTINE init_subspace_diag

  SUBROUTINE parameter_check(nme,MB)
    implicit none
    integer,intent(IN) :: nme,MB
    real(8) :: d_nme
    d_nme = ( dble(MB)*dble(MB)+dble(MB) )/2.0d0
    if ( abs(d_nme-nme) > 1.d-10 ) then
       write(*,*) "MB,nme,d_nme=",MB,nme,d_nme
       write(*,*) "MB may be too large"
       stop "stop@init_subspace_diag"
    end if
  END SUBROUTINE parameter_check


END MODULE subspace_diag_module
