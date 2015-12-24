MODULE bcset_0_module

  use rgrid_module, only: Igrid
  use parallel_module
  use bc_variables
  use watch_module, only: watchb_omp, time_bcfd
  use omp_variables, only: partition_job_omp

  implicit none

  PRIVATE
  PUBLIC :: bcset_0

CONTAINS


  SUBROUTINE bcset_0( ib1, ib2, ndepth, idir )

    implicit none
    integer,intent(IN) :: ib1,ib2,ndepth,idir
    integer :: a1,a2,a3,b1,b2,b3,nb,ns,ms,mt
    integer :: m,n,ndata,i1,i2,i3,ib,ierr
    integer :: c1,c2,c3,c4,d1,d2,d3,d4,irank,nreq,itags,itagr,ireq(36)
    integer :: i,l,e1,e2,e3,f1,f2,f3,fe1,fe12,fe123,dc1,dc12,dc123
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

    n = 1

!(1) [sender][2]-->[1][receiver]

    !call watchb_omp( ttmp )

!$OMP master
    ndata = fdinfo_recv(10,n,1)*nb
    !if ( ndata < 1 ) write(*,*) "bbb"
    irank = fdinfo_recv(8,n,1)
    itagr = 10
    nreq  = nreq + 1
    call mpi_irecv(rbuf(1,n,1),ndata,TYPE_MAIN,irank,itagr &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,1) )

    c1    = b1
    d1    = b1
    ndata = fdinfo_send(10,n,2)*nb
    !if ( ndata < 1 ) write(*,*) "ccc"

    c2 = fdinfo_send(3,n,2)
    d2 = fdinfo_send(4,n,2)
    c3 = fdinfo_send(5,n,2)
    d3 = fdinfo_send(6,n,2)
    c4 = ib1
    d4 = ib2

    e1 = c1
    f1 = d1
    e2 = c2
    f2 = d2
    e3 = c3
    f3 = d3

    fe1 = f1-e1+1
    fe12 = fe1*(f2-e2+1)
    fe123 = fe12*(f3-e3+1)

    call partition_job_omp( c1,d1,c2,d2,c3,d3,c4,d4 )

    !call watchb_omp( ttmp, time_bcfd(1,7) )

    do ib=c4,d4
       do i3=c3,d3
       do i2=c2,d2
       do i1=c1,d1
          i = 1+(i1-e1)+(i2-e2)*fe1+(i3-e3)*fe12+(ib-ib1)*fe123
          sbuf(i,n,2)=www(i1,i2,i3,ib)
       end do
       end do
       end do
    end do ! ib
!$OMP barrier

    !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
    irank = fdinfo_send(7,n,2)
    itags = 10
    nreq  = nreq + 1
    call mpi_isend(sbuf(1,n,2),ndata,TYPE_MAIN,irank,itags &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,2) )

!(3) [sender][4]-->[3][receiver]

!$OMP master
    ndata = fdinfo_recv(10,n,3)*nb
    !if ( ndata < 1 ) write(*,*) "bbb"
    irank = fdinfo_recv(8,n,3)
    itagr = 30
    nreq  = nreq + 1
    call mpi_irecv(rbuf(1,n,3),ndata,TYPE_MAIN,irank,itagr &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,1) )

    c2    = b2
    d2    = b2
    ndata = fdinfo_send(10,n,4)*nb
    !if ( ndata < 1 ) write(*,*) "ccc"

    c1 = fdinfo_send(1,n,4)
    d1 = fdinfo_send(2,n,4)
    c3 = fdinfo_send(5,n,4)
    d3 = fdinfo_send(6,n,4)
    c4 = ib1
    d4 = ib2

    e1 = c1
    f1 = d1
    e2 = c2
    f2 = d2
    e3 = c3
    f3 = d3

    fe1 = f1-e1+1
    fe12 = fe1*(f2-e2+1)
    fe123 = fe12*(f3-e3+1)

    call partition_job_omp( c1,d1,c2,d2,c3,d3,c4,d4 )

    !call watchb_omp( ttmp, time_bcfd(1,7) )

    do ib=c4,d4
       do i3=c3,d3
       do i2=c2,d2
       do i1=c1,d1
          i = 1+(i1-e1)+(i2-e2)*fe1+(i3-e3)*fe12+(ib-ib1)*fe123
          sbuf(i,n,4)=www(i1,i2,i3,ib)
       end do
       end do
       end do
    end do ! ib
!$OMP barrier

    !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
    irank = fdinfo_send(7,n,4)
    itags = 30
    nreq  = nreq + 1
    call mpi_isend(sbuf(1,n,4),ndata,TYPE_MAIN,irank,itags &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,2) )

!(5) [sender][6]-->[5][receiver]

!$OMP master
    ndata = fdinfo_recv(10,n,5)*nb
    !if ( ndata < 1 ) write(*,*) "bbb"
    irank = fdinfo_recv(8,n,5)
    itagr = 50
    nreq  = nreq + 1
    call mpi_irecv(rbuf(1,n,5),ndata,TYPE_MAIN,irank,itagr &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,1) )

    c3    = b3
    d3    = b3
    ndata = fdinfo_send(10,n,6)*nb
    !if ( ndata < 1 ) write(*,*) "ccc"

    c1 = fdinfo_send(1,n,6)
    d1 = fdinfo_send(2,n,6)
    c2 = fdinfo_send(3,n,6)
    d2 = fdinfo_send(4,n,6)
    c4 = ib1
    d4 = ib2

    e1 = c1
    f1 = d1
    e2 = c2
    f2 = d2
    e3 = c3
    f3 = d3

    fe1 = f1-e1+1
    fe12 = fe1*(f2-e2+1)
    fe123 = fe12*(f3-e3+1)

    call partition_job_omp( c1,d1,c2,d2,c3,d3,c4,d4 )

    !call watchb_omp( ttmp, time_bcfd(1,7) )

    do ib=c4,d4
       do i3=c3,d3
       do i2=c2,d2
       do i1=c1,d1
          i = 1+(i1-e1)+(i2-e2)*fe1+(i3-e3)*fe12+(ib-ib1)*fe123
          sbuf(i,n,6)=www(i1,i2,i3,ib)
       end do
       end do
       end do
    end do ! ib
!$OMP barrier

    !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
    irank = fdinfo_send(7,n,6)
    itags = 50
    nreq  = nreq + 1
    call mpi_isend(sbuf(1,n,6),ndata,TYPE_MAIN,irank,itags &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,2) )

!--------------------------------------------

!(2) [receiver][2]<--[1][sender]

!$OMP master
    ndata = fdinfo_recv(10,n,2)*nb
    !if ( ndata < 1 ) write(*,*) "bbb"
    irank = fdinfo_recv(8,n,2)
    itagr = 20
    nreq  = nreq + 1
    call mpi_irecv(rbuf(1,n,2),ndata,TYPE_MAIN,irank,itagr &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,1) )

    c1    = a1
    d1    = a1
    ndata = fdinfo_send(10,n,1)*nb
    !if ( ndata < 1 ) write(*,*) "ccc"

    c2 = fdinfo_send(3,n,1)
    d2 = fdinfo_send(4,n,1)
    c3 = fdinfo_send(5,n,1)
    d3 = fdinfo_send(6,n,1)
    c4 = ib1
    d4 = ib2

    e1 = c1
    f1 = d1
    e2 = c2
    f2 = d2
    e3 = c3
    f3 = d3

    fe1 = f1-e1+1
    fe12 = fe1*(f2-e2+1)
    fe123 = fe12*(f3-e3+1)

    call partition_job_omp( c1,d1,c2,d2,c3,d3,c4,d4 )

    !call watchb_omp( ttmp, time_bcfd(1,7) )

    do ib=c4,d4
       do i3=c3,d3
       do i2=c2,d2
       do i1=c1,d1
          i = 1+(i1-e1)+(i2-e2)*fe1+(i3-e3)*fe12+(ib-ib1)*fe123
          sbuf(i,n,1)=www(i1,i2,i3,ib)
       end do
       end do
       end do
    end do ! ib
!$OMP barrier

    !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
    irank = fdinfo_send(7,n,1)
    itags = 20
    nreq  = nreq + 1
    call mpi_isend(sbuf(1,n,1),ndata,TYPE_MAIN,irank,itags &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,2) )

!(4) [receiver][4]<--[3][sender]

!$OMP master
    ndata = fdinfo_recv(10,n,4)*nb
    !if ( ndata < 1 ) write(*,*) "bbb"
    irank = fdinfo_recv(8,n,4)
    itagr = 40
    nreq  = nreq + 1
    call mpi_irecv(rbuf(1,n,4),ndata,TYPE_MAIN,irank,itagr &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,1) )

    c2    = a2
    d2    = a2
    ndata = fdinfo_send(10,n,3)*nb
    !if ( ndata < 1 ) write(*,*) "ccc"

    c1 = fdinfo_send(1,n,3)
    d1 = fdinfo_send(2,n,3)
    c3 = fdinfo_send(5,n,3)
    d3 = fdinfo_send(6,n,3)
    c4 = ib1
    d4 = ib2

    e1 = c1
    f1 = d1
    e2 = c2
    f2 = d2
    e3 = c3
    f3 = d3

    fe1 = f1-e1+1
    fe12 = fe1*(f2-e2+1)
    fe123 = fe12*(f3-e3+1)

    call partition_job_omp( c1,d1,c2,d2,c3,d3,c4,d4 )

    !call watchb_omp( ttmp, time_bcfd(1,7) )

    do ib=c4,d4
       do i3=c3,d3
       do i2=c2,d2
       do i1=c1,d1
          i = 1+(i1-e1)+(i2-e2)*fe1+(i3-e3)*fe12+(ib-ib1)*fe123
          sbuf(i,n,3)=www(i1,i2,i3,ib)
       end do
       end do
       end do
    end do ! ib
!$OMP barrier

    !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
    irank = fdinfo_send(7,n,3)
    itags = 40
    nreq  = nreq + 1
    call mpi_isend(sbuf(1,n,3),ndata,TYPE_MAIN,irank,itags &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,2) )

!(6) [receiver][6]<--[5][sender]

!$OMP master
    ndata = fdinfo_recv(10,n,6)*nb
    !if ( ndata < 1 ) write(*,*) "bbb"
    irank = fdinfo_recv(8,n,6)
    itagr = 60
    nreq  = nreq + 1
    call mpi_irecv(rbuf(1,n,6),ndata,TYPE_MAIN,irank,itagr &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,1) )

    c3    = a3
    d3    = a3
    ndata = fdinfo_send(10,n,5)*nb
    !if ( ndata < 1 ) write(*,*) "ccc"

    c1 = fdinfo_send(1,n,5)
    d1 = fdinfo_send(2,n,5)
    c2 = fdinfo_send(3,n,5)
    d2 = fdinfo_send(4,n,5)
    c4 = ib1
    d4 = ib2

    e1 = c1
    f1 = d1
    e2 = c2
    f2 = d2
    e3 = c3
    f3 = d3

    fe1 = f1-e1+1
    fe12 = fe1*(f2-e2+1)
    fe123 = fe12*(f3-e3+1)

    call partition_job_omp( c1,d1,c2,d2,c3,d3,c4,d4 )

    !call watchb_omp( ttmp, time_bcfd(1,7) )

    do ib=c4,d4
       do i3=c3,d3
       do i2=c2,d2
       do i1=c1,d1
          i = 1+(i1-e1)+(i2-e2)*fe1+(i3-e3)*fe12+(ib-ib1)*fe123
          sbuf(i,n,5)=www(i1,i2,i3,ib)
       end do
       end do
       end do
    end do ! ib
!$OMP barrier

    !call watchb_omp( ttmp, time_bcfd(1,3) )

!$OMP master
    irank = fdinfo_send(7,n,5)
    itags = 60
    nreq  = nreq + 1
    call mpi_isend(sbuf(1,n,5),ndata,TYPE_MAIN,irank,itags &
                  ,comm_grid,ireq(nreq),ierr)
!$OMP end master

    !call watchb_omp( ttmp, time_bcfd(1,2) )

!$OMP master
    call mpi_waitall(nreq,ireq,istatus,ierr) ; nreq=0
!$OMP end master
!$OMP barrier

    !call watchb_omp( ttmp, time_bcfd(1,4) )

    do m=1,6

       c1=fdinfo_recv(1,n,m)
       d1=fdinfo_recv(2,n,m)
       c2=fdinfo_recv(3,n,m)
       d2=fdinfo_recv(4,n,m)
       c3=fdinfo_recv(5,n,m)
       d3=fdinfo_recv(6,n,m)

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

       dc1 = d1-c1+1
       dc12 = dc1*(d2-c2+1)
       dc123 = dc12*(d3-c3+1)

       do ib=ib1,ib2
!$OMP do
          do i3=c3,d3
          do i2=c2,d2
          do i1=c1,d1
!             i = 1 + i1-c1 + (i2-c2)*(d1-c1+1) &
!                   + (i3-c3)*(d1-c1+1)*(d2-c2+1) &
!                   + (ib-ib1)*(d1-c1+1)*(d2-c2+1)*(d3-c3+1) 
             i=1+(i1-c1)+(i2-c2)*dc1+(i3-c3)*dc12+(ib-ib1)*dc123
             www(i1,i2,i3,ib)=rbuf(i,n,m)
          end do
          end do
          end do
!$OMP end do
       end do ! ib

    end do ! m

    !call watchb_omp( ttmp, time_bcfd(1,5) )
    !call watchb_omp( ttmp_bc, time_bcfd(1,6) )

    return

  END SUBROUTINE bcset_0


END MODULE bcset_0_module
