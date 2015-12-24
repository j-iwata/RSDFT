MODULE bcset_1_module

  use rgrid_module, only: Igrid
  use parallel_module
  use bc_variables
  use watch_module, only: watchb_omp, time_bcfd

  implicit none

  PRIVATE
  PUBLIC :: bcset_1

CONTAINS


  SUBROUTINE bcset_1( ib1, ib2, ndepth, idir )

    implicit none
    integer,intent(IN) :: ib1,ib2,ndepth,idir
    integer :: a1,a2,a3,b1,b2,b3,nb,ns,ms,mt
    integer :: m,n,ndata,i1,i2,i3,ib,ierr
    integer :: c1,c2,c3,d1,d2,d3,irank,nreq,itags,itagr,ireq(36)
    integer :: i,l
    integer :: istatus(mpi_status_size,123)
    real(8) :: ttmp(2),ttmp_bc(2)

    !call watchb_omp( ttmp_bc )

    a1 = Igrid(1,1)
    b1 = Igrid(2,1)
    a2 = Igrid(1,2)
    b2 = Igrid(2,2)
    a3 = Igrid(1,3)
    b3 = Igrid(2,3)

    nb = ib2 - ib1 + 1

!$OMP master
    nreq=0
!$OMP end master

!(1) [sender][2]-->[1][receiver]

    if ( idir == 0 .or. idir == 1 .or. idir == 2 ) then

       !call watchb_omp( ttmp )

!$OMP master
       do n=1,n_neighbor(1)

          if ( ndepth == 1 ) then
             ndata = fdinfo_recv(10,n,1)*nb
          else
             ndata = fdinfo_recv(9,n,1)*nb
          end if

          if ( ndata < 1 ) cycle

          irank = fdinfo_recv(8,n,1)
          itagr = 10
          nreq  = nreq + 1

          call mpi_irecv(rbuf(1,n,1),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)

       end do ! n
!$OMP end master

       !call watchb_omp( ttmp, time_bcfd(1,1) )

       do n=1,n_neighbor(2)

          if ( ndepth == 1 ) then
             c1    = b1
             d1    = b1
             ndata = fdinfo_send(10,n,2)*nb
          else
             c1    = fdinfo_send(1,n,2)
             d1    = fdinfo_send(2,n,2)
             ndata = fdinfo_send(9,n,2)*nb
          end if

          c2 = fdinfo_send(3,n,2)
          d2 = fdinfo_send(4,n,2)
          c3 = fdinfo_send(5,n,2)
          d3 = fdinfo_send(6,n,2)

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
                sbuf(i,n,2)=www(i1,i2,i3,ib)
             end do
             end do
             end do
!$OMP end do
          end do ! ib

          !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
          irank = fdinfo_send(7,n,2)
          itags = 10
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,2),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master

          !call watchb_omp( ttmp, time_bcfd(1,2) )

       end do ! n

    end if

!(3) [sender][4]-->[3][receiver]

    if ( idir == 0 .or. idir == 3 .or. idir == 4 ) then

       !call watchb_omp( ttmp )

!$OMP master
       do n=1,n_neighbor(3)

          if ( ndepth == 1 ) then
             ndata = fdinfo_recv(10,n,3)*nb
          else
             ndata = fdinfo_recv(9,n,3)*nb
          end if

          if ( ndata < 1 ) cycle

          irank = fdinfo_recv(8,n,3)
          itagr = 30
          nreq  = nreq + 1

          call mpi_irecv(rbuf(1,n,3),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)

       end do ! n
!$OMP end master

       !call watchb_omp( ttmp, time_bcfd(1,1) )

       do n=1,n_neighbor(4)

          if ( ndepth == 1 ) then
             c2    = b2
             d2    = b2
             ndata = fdinfo_send(10,n,4)*nb
          else
             c2    = fdinfo_send(3,n,4)
             d2    = fdinfo_send(4,n,4)
             ndata = fdinfo_send(9,n,4)*nb
          end if

          c1 = fdinfo_send(1,n,4)
          d1 = fdinfo_send(2,n,4)
          c3 = fdinfo_send(5,n,4)
          d3 = fdinfo_send(6,n,4)

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
                sbuf(i,n,4)=www(i1,i2,i3,ib)
             end do
             end do
             end do
!$OMP end do

          end do ! ib

          !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
          irank = fdinfo_send(7,n,4)
          itags = 30
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,4),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master

          !call watchb_omp( ttmp, time_bcfd(1,2) )

       end do ! n

    end if

!(5) [sender][6]-->[5][receiver]

    if ( idir == 0 .or. idir == 5 .or. idir == 6 ) then

       !call watchb_omp( ttmp )

!$OMP master
       do n=1,n_neighbor(5)

          if ( ndepth == 1 ) then
             ndata = fdinfo_recv(10,n,5)*nb
          else
             ndata = fdinfo_recv(9,n,5)*nb
          end if

          if ( ndata < 1 ) cycle

          irank = fdinfo_recv(8,n,5)
          itagr = 50
          nreq  = nreq + 1

          call mpi_irecv(rbuf(1,n,5),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)

       end do ! n
!$OMP end master

       !call watchb_omp( ttmp, time_bcfd(1,1) )

       do n=1,n_neighbor(6)

          if ( ndepth == 1 ) then
             c3    = b3
             d3    = b3
             ndata = fdinfo_send(10,n,6)*nb
          else
             c3    = fdinfo_send(5,n,6)
             d3    = fdinfo_send(6,n,6)
             ndata = fdinfo_send(9,n,6)*nb
          end if

          c1 = fdinfo_send(1,n,6)
          d1 = fdinfo_send(2,n,6)
          c2 = fdinfo_send(3,n,6)
          d2 = fdinfo_send(4,n,6)

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
                sbuf(i,n,6)=www(i1,i2,i3,ib)
             end do
             end do
             end do
!$OMP end do
          end do

          !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
          irank = fdinfo_send(7,n,6)
          itags = 50
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,6),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master

          !call watchb_omp( ttmp, time_bcfd(1,2) )

       end do ! n

    end if

!--------------------------------------------

!(2) [receiver][2]<--[1][sender]

    if ( idir == 0 .or. idir == 1 .or. idir == 2 ) then

       !call watchb_omp( ttmp )

!$OMP master
       do n=1,n_neighbor(2)

          if ( ndepth == 1 ) then
             ndata = fdinfo_recv(10,n,2)*nb
          else
             ndata = fdinfo_recv(9,n,2)*nb
          end if

          if ( ndata < 1 ) cycle

          irank = fdinfo_recv(8,n,2)
          itagr = 20
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,2),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)

       end do ! n
!$OMP end master

       !call watchb_omp( ttmp, time_bcfd(1,1) )

       do n=1,n_neighbor(1)

          if ( ndepth == 1 ) then
             c1    = a1
             d1    = a1
             ndata = fdinfo_send(10,n,1)*nb
          else
             c1    = fdinfo_send(1,n,1)
             d1    = fdinfo_send(2,n,1)
             ndata = fdinfo_send(9,n,1)*nb
          end if

          c2 = fdinfo_send(3,n,1)
          d2 = fdinfo_send(4,n,1)
          c3 = fdinfo_send(5,n,1)
          d3 = fdinfo_send(6,n,1)

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
                sbuf(i,n,1)=www(i1,i2,i3,ib)
             end do
             end do
             end do
!$OMP end do
          end do ! ib

          !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
          irank = fdinfo_send(7,n,1)
          itags = 20
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,1),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master

          !call watchb_omp( ttmp, time_bcfd(1,2) )

       end do ! n

    end if

!(4) [receiver][4]<--[3][sender]

    if ( idir == 0 .or. idir == 3 .or. idir == 4 ) then

       !call watchb_omp( ttmp )

!$OMP master
       do n=1,n_neighbor(4)

          if ( ndepth == 1 ) then
             ndata = fdinfo_recv(10,n,4)*nb
          else
             ndata = fdinfo_recv(9,n,4)*nb
          end if

          if ( ndata < 1 ) cycle

          irank = fdinfo_recv(8,n,4)
          itagr = 40
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,4),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)

       end do ! n
!$OMP end master

       !call watchb_omp( ttmp, time_bcfd(1,1) )

       do n=1,n_neighbor(3)

          if ( ndepth == 1 ) then
             c2    = a2
             d2    = a2
             ndata = fdinfo_send(10,n,3)*nb
          else
             c2    = fdinfo_send(3,n,3)
             d2    = fdinfo_send(4,n,3)
             ndata = fdinfo_send(9,n,3)*nb
          end if

          c1 = fdinfo_send(1,n,3)
          d1 = fdinfo_send(2,n,3)
          c3 = fdinfo_send(5,n,3)
          d3 = fdinfo_send(6,n,3)

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
                sbuf(i,n,3)=www(i1,i2,i3,ib)
             end do
             end do
             end do
!$OMP end do
          end do ! ib

          !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
          irank = fdinfo_send(7,n,3)
          itags = 40
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,3),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master

          !call watchb_omp( ttmp, time_bcfd(1,2) )

       end do ! n

    end if

!(6) [receiver][6]<--[5][sender]

    if ( idir == 0 .or. idir == 5 .or. idir == 6 ) then

       !call watchb_omp( ttmp )

!$OMP master
       do n=1,n_neighbor(6)

          if ( ndepth == 1 ) then
             ndata = fdinfo_recv(10,n,6)*nb
          else
             ndata = fdinfo_recv(9,n,6)*nb
          end if

          if ( ndata < 1 ) cycle
          irank = fdinfo_recv(8,n,6)
          itagr = 60
          nreq  = nreq + 1
          call mpi_irecv(rbuf(1,n,6),ndata,TYPE_MAIN,irank,itagr &
                        ,comm_grid,ireq(nreq),ierr)

       end do ! n
!$OMP end master

       !call watchb_omp( ttmp, time_bcfd(1,1) )

       do n=1,n_neighbor(5)

          if ( ndepth == 1 ) then
             c3    = a3
             d3    = a3
             ndata = fdinfo_send(10,n,5)*nb
          else
             c3    = fdinfo_send(5,n,5)
             d3    = fdinfo_send(6,n,5)
             ndata = fdinfo_send(9,n,5)*nb
          end if

          c1 = fdinfo_send(1,n,5)
          d1 = fdinfo_send(2,n,5)
          c2 = fdinfo_send(3,n,5)
          d2 = fdinfo_send(4,n,5)

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
                sbuf(i,n,5)=www(i1,i2,i3,ib)
             end do
             end do
             end do
!$OMP end do
          end do ! ib

          !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
          irank = fdinfo_send(7,n,5)
          itags = 60
          nreq  = nreq + 1
          call mpi_isend(sbuf(1,n,5),ndata,TYPE_MAIN,irank,itags &
                        ,comm_grid,ireq(nreq),ierr)
!$OMP end master

          !call watchb_omp( ttmp, time_bcfd(1,2) )

       end do ! b

    end if

    !call watchb_omp( ttmp )

!$OMP master
    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0
!$OMP end master
!$OMP barrier

    !call watchb_omp( ttmp, time_bcfd(1,4) )

    do m=1,6
       do n=1,n_neighbor(m)

          c1=fdinfo_recv(1,n,m)
          d1=fdinfo_recv(2,n,m)
          c2=fdinfo_recv(3,n,m)
          d2=fdinfo_recv(4,n,m)
          c3=fdinfo_recv(5,n,m)
          d3=fdinfo_recv(6,n,m)

          if ( Md>ndepth .and. ndepth==1 ) then
             if ( fdinfo_recv(10,n,m)<1 ) cycle
             select case(m)
             case(1)
                c1=a1-1
                d1=c1
             case(2)
                c1=b1+1
                d1=c1
             case(3)
                c2=a2-1
                d2=c2
             case(4)
                c2=b2+1
                d2=c2
             case(5)
                c3=a3-1
                d3=c3
             case(6)
                c3=b3+1
                d3=c3
             end select
          else
             if ( fdinfo_recv(9,n,m)<1 ) cycle
          end if

          do ib=ib1,ib2
!$OMP do
             do i3=c3,d3
             do i2=c2,d2
             do i1=c1,d1
                i = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
                      + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
                      + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
                www(i1,i2,i3,ib)=rbuf(i,n,m)
             end do
             end do
             end do
!$OMP end do
          end do

       end do ! n
    end do ! m

    !call watchb_omp( ttmp, time_bcfd(1,5) )
    !call watchb_omp( ttmp_bc, time_bcfd(1,6) )

    return

  END SUBROUTINE bcset_1


END MODULE bcset_1_module
