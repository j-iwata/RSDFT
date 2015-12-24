MODULE bcset_3_module

  use parallel_module
  use bc_variables
  use watch_module, only: watchb_omp, time_bcfd

  implicit none

  PRIVATE
  PUBLIC :: bcset_3

CONTAINS


  SUBROUTINE bcset_3( ib1, ib2, ndepth, idir_in )

    implicit none
    integer,intent(IN) :: ib1,ib2,ndepth,idir_in
    integer :: a1,a2,a3,b1,b2,b3,nb,ns,ms,mt
    integer :: m,n,ndata,i1,i2,i3,ib,ierr
    integer :: c1,c2,c3,d1,d2,d3,irank,nreq,itags,itagr,ireq(36)
    integer :: i,l,idir,jdir
    integer :: istatus(mpi_status_size,123)
    real(8) :: ttmp(2),ttmp_bc(2)

    !call watchb_omp( ttmp_bc )

    nb = ib2 - ib1 + 1

!$OMP master
    nreq=0
!$OMP end master

    do idir=1,6

       if ( mod(idir,2) == 0 ) then
          jdir = idir - 1
       else
          jdir = idir + 1
       end if

       !call watchb_omp( ttmp )

!$OMP master
       do n=1,n_neighbor(idir)

          ndata = fdinfo_recv(9,n,idir)*nb

          if ( ndata < 1 ) cycle

          irank = fdinfo_recv(8,n,idir)
          itagr = 10*idir
          nreq  = nreq + 1

          call mpi_irecv(rbuf(1,n,idir),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)

       end do ! n
!$OMP end master

       !call watchb_omp( ttmp, time_bcfd(1,1) )

       do n=1,n_neighbor(jdir)
          c1 = fdinfo_send(1,n,jdir)
          d1 = fdinfo_send(2,n,jdir)
          c2 = fdinfo_send(3,n,jdir)
          d2 = fdinfo_send(4,n,jdir)
          c3 = fdinfo_send(5,n,jdir)
          d3 = fdinfo_send(6,n,jdir)
          ndata = fdinfo_send(9,n,jdir)*nb

          if ( ndata < 1 ) cycle

          !call watchb_omp( ttmp )

          do ib=ib1,ib2
!$OMP do
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                i = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
                      + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
                      + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
                sbuf(i,n,jdir)=www(i1,i2,i3,ib)
             end do
             end do
             end do
!$OMP end do
          end do ! ib

          !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
          irank = fdinfo_send(7,n,jdir)
          itags = 10*idir
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,jdir),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master

          !call watchb_omp( ttmp, time_bcfd(1,2) )

       end do ! n

    end do ! idir

    !call watchb_omp( ttmp )

!$OMP master
    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0
!$OMP end master
!$OMP barrier

    !call watchb_omp( ttmp, time_bcfd(1,4) )

    do idir=1,6
       do n=1,n_neighbor(idir)

          c1=fdinfo_recv(1,n,idir)
          d1=fdinfo_recv(2,n,idir)
          c2=fdinfo_recv(3,n,idir)
          d2=fdinfo_recv(4,n,idir)
          c3=fdinfo_recv(5,n,idir)
          d3=fdinfo_recv(6,n,idir)

          if ( fdinfo_recv(9,n,idir)<1 ) cycle

          do ib=ib1,ib2
!$OMP do
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                i = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
                      + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
                      + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
                www(i1,i2,i3,ib)=rbuf(i,n,idir)
             end do
             end do
             end do
!$OMP end do
          end do

       end do ! n
    end do ! idir

    !call watchb_omp( ttmp, time_bcfd(1,5) )
    !call watchb_omp( ttmp_bc, time_bcfd(1,6) )

    return

  END SUBROUTINE bcset_3


END MODULE bcset_3_module
