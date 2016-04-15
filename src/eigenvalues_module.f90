MODULE eigenvalues_module

  use bz_module, only: bzinfo, construct_bzinfo_bz

  implicit none

  PRIVATE
  PUBLIC :: eigv
  PUBLIC :: construct_eigenvalues
  PUBLIC :: write_eigenvalues

  TYPE eigv
     integer :: nband, nkpnt, nspin
     real(8),allocatable :: eigv(:,:,:)
     type(bzinfo) :: bz
  END type eigv

  integer :: unit=20
  character(16) :: file_name="eigenvalues"

CONTAINS


  SUBROUTINE construct_eigenvalues( nb, nk, ns, esp_in, esp )
    implicit none
    integer,intent(IN) :: nb,nk,ns
    real(8),intent(IN) :: esp_in(:,:,:)
    type(eigv),intent(INOUT) :: esp

    esp%nband = nb
    esp%nkpnt = nk
    esp%nspin = ns

    if ( .not.allocated(esp%eigv) ) then
       allocate( esp%eigv(nb,nk,ns) ) ; esp%eigv=0.0d0
       call construct_bzinfo_bz( esp%bz )
    end if

    esp%eigv = 0.0d0

    esp%eigv(1:nb,1:nk,1:ns) = esp_in(1:nb,1:nk,1:ns)


  END SUBROUTINE construct_eigenvalues


  SUBROUTINE write_eigenvalues( esp )
    implicit none
    type(eigv),intent(IN) :: esp
    integer :: n,k,s

    open( unit, file=file_name )

    do k=1,esp%nkpnt

       write(unit,'(1x,"Nband,Nspin,k=",3i8)') esp%nband,esp%nspin,k
       write(unit,'(1x,4f20.15)') dble(esp%bz%kpt(1:3,k))/dble(esp%bz%nk) &
                                 ,esp%bz%weight(k)

       do n=1,esp%nband
          write(unit,'(1x,i6,2(f20.15,1x,g15.7,1x))') &
               n, ( esp%eigv(n,k,s),s=1,esp%nspin )
       end do ! n

    end do ! k

    close( unit )

  END SUBROUTINE write_eigenvalues

END MODULE eigenvalues_module
