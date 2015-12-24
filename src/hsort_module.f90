MODULE hsort_module

  implicit none

  PRIVATE
  PUBLIC :: indexx

CONTAINS


  SUBROUTINE indexx( n, arr, indx )
    implicit none
    integer,intent(IN)  :: n
    real(8),intent(IN)  :: arr(n)
    integer,intent(OUT) :: indx(n)
    call indexx_0( n, arr, indx )
    !call indexx_1( n, arr, indx )
    !call ascending_sort_indexx_2( n, arr, indx )
  END SUBROUTINE indexx


  SUBROUTINE indexx_0(n,arr,indx)
    implicit none
    integer,intent(IN)  :: n
    real(8),intent(IN)  :: arr(n)
    integer,intent(OUT) :: indx(n)
    integer,parameter :: m=7,nstck=50
    integer :: i,j,k,ll,jstck,jndx,lndx,ir,istck(nstck),loop0,loop1
    real(8) :: arrj,arrl

    do i=1,n
       indx(i)=i
    end do
    jstck=0
    ll=1
    ir=n
    do loop0=1,n
    if ( ir-ll < m ) then
       do j=ll+1,ir
          jndx=indx(j)
          arrj=arr(jndx)
          do i=j-1,1,-1
             if ( arr(indx(i)) <= arrj ) exit
             indx(i+1)=indx(i)
          end do
          if ( i < 1 ) i=0
          indx(i+1)=jndx
       end do
       if ( jstck == 0 ) return
       ir=istck(jstck)
       ll=istck(jstck-1)
       jstck=jstck-2
    else
       k=(ll+ir)/2
       call iswap( indx(k), indx(ll+1) )
       if ( arr(indx(ll+1)) > arr(indx(ir)) ) then
          call iswap( indx(ll+1), indx(ir) )
       end if
       if ( arr(indx(ll)) > arr(indx(ir)) ) then
          call iswap( indx(ll), indx(ir) )
       end if
       if ( arr(indx(ll+1)) > arr(indx(ll)) ) then
          call iswap( indx(ll+1), indx(ll) )
       end if
       i=ll+1
       j=ir
       lndx=indx(ll)
       arrl=arr(lndx)
       do loop1=1,n
          i=i+1
          if ( arr(indx(i)) < arrl ) cycle
          j=j-1
          do while ( arr(indx(j)) > arrl )
             j=j-1
          end do
          if ( j < i ) exit
          call iswap( indx(i), indx(j) )
       end do ! loop1
       indx(ll)=indx(j)
       indx(j )=lndx
       jstck=jstck+2
       if ( jstck > nstck ) call stop_program( "stop@indexx_0" )
       if ( ir-i+1 >= j-ll ) then
          istck(jstck)  =ir
          istck(jstck-1)=i
          ir=j-1
       else
          istck(jstck)  =j-1
          istck(jstck-1)=ll
          ll=i
       end if
    end if
    end do ! loop0
  END SUBROUTINE indexx_0


  SUBROUTINE indexx_1( ndat, arry, indx )
    implicit none
    integer,intent(IN)  :: ndat
    real(8),intent(IN)  :: arry(ndat)
    integer,intent(OUT) :: indx(ndat)
    integer :: ilayer,iboss,jlayer,jboss,kboss
    integer :: i,i1,i2,max_layer,loop,ndat_max
    integer,allocatable :: irslt(:)
    logical :: flag_boss_change
    real(8) :: arry_max

    arry_max = maxval( arry )

    max_layer = log( dble(ndat+1) )/log( 2.0d0 )

    ndat_max = 2**max_layer - 1
    if ( ndat_max < ndat ) max_layer = max_layer + 1
    ndat_max = 2**max_layer - 1
    if ( ndat_max < ndat ) stop "stop(1)"

    do i=1,ndat
       indx(i)=i
    end do

! ---

    do ilayer=2,max_layer

       loop_boss : do iboss=2**(ilayer-2),2**(ilayer-1)-1

          if ( iboss > ndat ) exit loop_boss

          flag_boss_change = .false.

          do i=2*iboss,2*iboss+1

             if ( i > ndat ) exit loop_boss

             if ( arry(indx(i)) < arry(indx(iboss)) ) then
                call iswap( indx(iboss), indx(i) )
                flag_boss_change = .true.
             end if

             if ( flag_boss_change ) then
                kboss = iboss
                do jlayer=ilayer-1,2,-1
                   do jboss=2**(jlayer-2),2**(jlayer-1)-1
                      if ( arry(indx(kboss)) < arry(indx(jboss)) ) then
                         call iswap( indx(jboss), indx(kboss) )
                         kboss = jboss
                         exit
                      end if
                   end do ! jboss
                end do ! jlayer
                if ( arry(indx(kboss)) < arry(indx(1)) ) then
                   call iswap( indx(1), indx(kboss) )
                end if
             end if

          end do ! i

       end do loop_boss

    end do ! ilayer

! ---

    allocate( irslt(ndat) ) ; irslt=0

    do loop=1,ndat

       irslt(loop)=indx(1)
       indx(1)=-1

       jboss=1
       do ilayer=2,max_layer

          i1=2*jboss
          i2=2*jboss+1

          if ( i1 > ndat ) then
             exit
          else if ( i1 <= ndat .and. i2 > ndat ) then
             call iswap( indx(jboss), indx(i1) )
             exit
          end if

          if ( indx(i1) < 0 .and. indx(i2) < 0 ) then
             exit
          else if ( indx(i1) > 0 .and. indx(i2) < 0 ) then
             call iswap( indx(jboss), indx(i1) )
             jboss=i1
          else if ( indx(i1) < 0 .and. indx(i2) > 0 ) then
             call iswap( indx(jboss), indx(i2) )
             jboss=i2
          else if ( arry(indx(i1)) < arry(indx(i2)) ) then
             call iswap( indx(jboss), indx(i1) )
             jboss=i1
          else
             call iswap( indx(jboss), indx(i2) )
             jboss=i2
          end if

       end do ! ilayer

    end do ! loop

    indx(:)=irslt(:)

    deallocate( irslt )

  END SUBROUTINE indexx_1

  SUBROUTINE iswap( i, j )
    implicit none
    integer,intent(INOUT) :: i,j
    integer :: k
    k=i
    i=j
    j=k
  END SUBROUTINE iswap


  SUBROUTINE ascending_sort_indexx_2( ndat, arry, indx )
    implicit none
    integer,intent(IN)  :: ndat
    real(8),intent(IN)  :: arry(ndat)
    integer,intent(OUT) :: indx(ndat)
    integer :: ipos, jpos,i,j,ndati,ndatj,loop,ipivo
    real(8) :: apivo

    do i=1,ndat
       indx(i) = i
    end do

    if ( ndat <= 1 ) then
       return
    else if ( ndat == 2 ) then
       if ( arry(2) < arry(1) ) then
          indx(1) = 2
          indx(2) = 1
       end if
       return
    end if

    apivo = arry(1)
    if ( apivo < arry(2) ) apivo = arry(2)

    ipos = 1
    jpos = ndat

! ---

    do loop=1,ndat

    if ( ipos == jpos ) then
       if ( arry(indx(ipos)) >= apivo ) then
          ipos = ipos - 1
       else if ( arry(indx(jpos)) < apivo ) then
          jpos = jpos + 1
       end if
       exit
    end if

    do i=ipos,jpos-1
       if ( arry(indx(i)) >= apivo ) then
          ipos = i
          exit
       end if
    end do

    do j=jpos,ipos+1,-1
       if ( arry(indx(j)) < apivo ) then
          jpos = j
          exit
       end if
    end do

    if ( i /= ipos .and. j == jpos ) then
       ipos = jpos
       jpos = jpos + 1
       exit
    else if ( i == ipos .and. j /= jpos ) then
       ipos = ipos - 1
       jpos = ipos + 1
       exit
    else if ( i /= ipos .and. j /= jpos ) then
       jpos = ipos + 1
       exit
    end if

    if ( arry(indx(ipos)) > arry(indx(jpos)) ) then
       i=indx(ipos)
       indx(ipos)=indx(jpos)
       indx(jpos)=i
    end if

    if ( ipos+1 == jpos ) exit

    ipos = ipos + 1
    jpos = jpos - 1

    end do

    if ( ipos+1 == jpos ) then
       ndati = ipos
       ndatj = ndat-jpos+1
       ipivo = 0
       call indexx_2_sub( ipivo, ndat, ndati, arry, indx )
       call indexx_2_sub( ipivo, ndat, ndatj, arry, indx(jpos) )
    else
       stop "error!:stop@ascending_sort_indexx_2"
    end if

  END SUBROUTINE ascending_sort_indexx_2

  RECURSIVE SUBROUTINE indexx_2_sub( ipivo_in, ndat, ndat1, arry, indx )
    implicit none
    integer,intent(IN) :: ipivo_in, ndat, ndat1
    real(8),intent(IN) :: arry(ndat)
    integer,intent(INOUT) :: indx(ndat1)
    integer :: i,j,ipos,jpos,ndati,ndatj,loop,ipivo
    real(8) :: apivo,apivo0,apivo1,err

    if ( ndat1 <= 1 ) then
       return
    else if ( ndat1 == 2 ) then
       if ( arry(indx(2)) < arry(indx(1)) ) then
          i=indx(1)
          indx(1) = indx(2)
          indx(2) = i
       end if
       return
    end if

    if ( ipivo_in == 0 ) then
       apivo = arry(indx(1))
       if ( apivo < arry(indx(2)) ) apivo = arry(indx(2))
    else if ( ipivo_in == 1 ) then
       apivo0 = 1.d100
       apivo1 =-1.d100
       do i=1,ndat1
          apivo = arry(indx(i))
          apivo0 = min( apivo, apivo0 )
          apivo1 = max( apivo, apivo1 )
       end do
       apivo = 0.5d0*( apivo0 + apivo1 )
    else
       stop "error!: stop@indexx_sub_2(1)"
    end if
       
    ipos = 1
    jpos = ndat1

! ---

    do loop=1,ndat1

    if ( ipos == jpos ) then
       if ( arry(indx(ipos)) >= apivo ) then
          ipos = ipos - 1
       else if ( arry(indx(jpos)) < apivo ) then
          jpos = jpos + 1
       end if
       exit
    end if

    do i=ipos,jpos-1
       if ( arry(indx(i)) >= apivo ) then
          ipos = i
          exit
       end if
    end do

    do j=jpos,ipos+1,-1
       if ( arry(indx(j)) < apivo ) then
          jpos = j
          exit
       end if
    end do

    if ( i /= ipos .and. j == jpos ) then
       ipos = jpos
       jpos = jpos + 1
       exit
    else if ( i == ipos .and. j /= jpos ) then
       ipos = ipos - 1
       jpos = ipos + 1
       exit
    else if ( i /= ipos .and. j /= jpos ) then
       jpos = ipos + 1
       exit
    end if

    if ( arry(indx(ipos)) > arry(indx(jpos)) ) then
       i=indx(ipos)
       indx(ipos)=indx(jpos)
       indx(jpos)=i
    end if

    if ( ipos+1 == jpos ) exit

    ipos = ipos + 1
    jpos = jpos - 1

    end do

    if ( ipos+1 == jpos ) then

       ndati = ipos
       ndatj = ndat1-jpos+1
       ipivo = 0

       if ( ipos == 0 .or. jpos == 0 ) then
          do i=1,ndat1
             err=abs(arry(indx(i))-apivo)
             if ( err > 1.d-12 ) exit
          end do
          if ( i > ndat1 ) return
          ipivo=1
       end if

       if ( ipos > 1 ) then
          call indexx_2_sub( ipivo,ndat,ndati,arry,indx )
       end if

       if ( jpos >= 1 .and. jpos < ndat1 ) then
          call indexx_2_sub( ipivo,ndat,ndatj,arry,indx(jpos) )
       end if

    else

       stop "error!: stop@indexx_2_sub"

    end if

  END SUBROUTINE indexx_2_sub


END MODULE hsort_module
