MODULE force_sub_module

  use parallel_module, only: myrank
  use polint_module
  use spline_module

  implicit none

  real(8),allocatable :: SH_Y1(:,:,:,:),SH_Y2(:,:,:,:),SH_Y3(:,:,:,:)
  real(8),allocatable :: C_ijLM(:,:,:,:,:,:)
  logical :: isFirstSHY=.true.,isFirstC=.true.
  logical :: isFirst_irad=.true.

  integer :: a1b,a2b,a3b,ab1,ab2,ab3
  integer :: ML1,ML2,ML3
  real(8) :: c1,c2,c3
  real(8),allocatable :: y2a(:,:,:),y2b(:,:,:)
#ifndef _SPLINE_
  integer,allocatable :: irad(:,:)
  integer :: M_irad
#endif
  logical,allocatable :: isAtomInThisNode(:)
  logical :: noAtomHere
  integer :: iatom,iL,im,iorb,ikind
  integer :: ik1,ik2,ikk1,lma1,lma2,l_1,l_2,m_1,m_2,iorb1,iorb2
  real(8) :: Rx,Ry,Rz
  real(8) :: x,y,z,r
  real(8) :: tmp0
  real(8) :: maxerr

CONTAINS
!-------------------------------------------------------------------------
  SUBROUTINE setLocalIndexForBoundary(Igrid)
    implicit none
    integer,intent(IN) :: Igrid(1:2,0:3)
    a1b=Igrid(1,1)
    a2b=Igrid(1,2)
    a3b=Igrid(1,3)
    ab1=Igrid(2,1)-Igrid(1,1)+1
    ab2=Igrid(2,2)-Igrid(1,2)+1
    ab3=Igrid(2,3)-Igrid(1,3)+1
  END SUBROUTINE setLocalIndexForBoundary
!-------------------------------------------------------------------------
  SUBROUTINE setConstGridWidth(Ngrid)
    implicit none
    integer,intent(IN) :: Ngrid(0:3)
    ML1 = Ngrid(1)
    ML2 = Ngrid(2)
    ML3 = Ngrid(3)
    c1=1.d0/ML1
    c2=1.d0/ML2
    c3=1.d0/ML3
  END SUBROUTINE setConstGridWidth
!-------------------------------------------------------------------------
#ifndef _SPLINE_
  SUBROUTINE setIndexForAtomCenteredGrid(Nelement_,NRps,rad1)
    implicit none
    integer,intent(IN) :: Nelement_
    integer,allocatable,intent(IN) :: NRps(:,:)
    real(8),allocatable,intent(IN) :: rad1(:,:)
    integer :: ik,ir,i,m
    integer :: NRc
    if (.not.isFirst_irad) then
      return
    endif
    isFirst_irad=.false.
!$OMP single
    allocate( irad(0:3000,Nelement_) )
    irad=0
    M_irad=0
    do ik=1,Nelement_
      NRc=maxval( NRps(:,ik) )
      NRc=min( 3000, NRc )
      m=0
      irad(0,ik)=1
      do ir=1,NRc
        m=int(100.d0*rad1(ir,ik))+1
        irad( m,ik )=ir
      end do
      ir=irad(0,ik)
      do i=1,m
        if ( irad(i,ik)==0 ) then
          irad(i,ik)=ir
          cycle
        end if
        ir=irad(i,ik)
      end do
      irad(m+1:,ik)=ir
      M_irad=max(M_irad,m)
    end do
!$OMP end single
  END SUBROUTINE setIndexForAtomCenteredGrid
#endif
!-------------------------------------------------------------------------
  SUBROUTINE setLocalAtoms(maxAtom,aa_atom,Igrid)
    implicit none
    integer,intent(IN) :: maxAtom
    real(8),intent(IN) :: aa_atom(3,maxAtom)
    integer,intent(IN) :: Igrid(2,0:3)
    integer :: ia
    integer :: i1,i2,i3
    real(8) :: k1,k2,k3
    if (.not. allocated(isAtomInThisNode)) then
      allocate( isAtomInThisNode(maxAtom) ) ; isAtomInThisNode=.false.
    endif
!$OMP parallel
!$OMP do private( i1,i2,i3,k1,k2,k3 )
    do ia=1,maxAtom
      i1 = nint( aa_atom(1,ia)*ML1 )
      i2 = nint( aa_atom(2,ia)*ML2 )
      i3 = nint( aa_atom(3,ia)*ML3 )
      k1 = i1/ML1 ; if ( i1<0 ) k1=(i1+1)/ML1-1
      k2 = i2/ML2 ; if ( i2<0 ) k2=(i2+1)/ML2-1
      k3 = i3/ML3 ; if ( i3<0 ) k3=(i3+1)/ML3-1
      i1 = i1 - k1*ML1
      i2 = i2 - k2*ML2
      i3 = i3 - k3*ML3
      if ( Igrid(1,1) <= i1 .and. i1 <= Igrid(2,1) .and. &
          Igrid(1,2) <= i2 .and. i2 <= Igrid(2,2) .and. &
          Igrid(1,3) <= i3 .and. i3 <= Igrid(2,3) ) then
        isAtomInThisNode(ia)=.true.
      end if
    end do
!$OMP end do
!$OMP end parallel
  END SUBROUTINE setLocalAtoms
!-------------------------------------------------------------------------
  SUBROUTINE getAtomInfoFrom_lma(lma)
    use ps_nloc2_variables, only: amap,lmap,mmap,iorbmap
    use atom_module, only: ki_atom
    implicit none
    integer,intent(IN) :: lma
    noAtomHere=.false.
    iatom = amap(lma)
    if (iatom<=0) then
      noAtomHere=.true.
      return
    endif
    iL    = lmap(lma)
    im    = mmap(lma)
    iorb  = iorbmap(lma)
    ikind = ki_atom(iatom)
    return
  END SUBROUTINE getAtomInfoFrom_lma
!-------------------------------------------------------------------------
  SUBROUTINE getAtomInfoFrom_iqr(iqr)
    use ps_nloc2_variables, only: amap,lmap,iorbmap,mmap
    use var_ps_member_g, only: k1_to_k2, nzqr_pair, k1map
    use atom_module, only: ki_atom
    implicit none
    integer,intent(IN) :: iqr
    lma1  = nzqr_pair(iqr,1)
    lma2  = nzqr_pair(iqr,2)
    iatom = amap(lma1)
    if ( iatom <= 0 ) then
      noAtomHere = .true.
      return
    else
       noAtomHere = .false.
    end if
    ikind = ki_atom(iatom)
    ik1   = k1map(iqr)
    ik2   = k1_to_k2(ik1,ikind)
    l_1   = lmap(lma1)
    m_1   = mmap(lma1)
    l_2   = lmap(lma2)
    m_2   = mmap(lma2)
    iorb1 = iorbmap(lma1)
    iorb2 = iorbmap(lma2)
    return
  END SUBROUTINE getAtomInfoFrom_iqr
!-------------------------------------------------------------------------
  SUBROUTINE getAtomPosition_real
    use atom_module, only: aa_atom
    use aa_module, only: aa
    implicit none
    Rx=aa(1,1)*aa_atom(1,iatom)+aa(1,2)*aa_atom(2,iatom)+aa(1,3)*aa_atom(3,iatom)
    Ry=aa(2,1)*aa_atom(1,iatom)+aa(2,2)*aa_atom(2,iatom)+aa(2,3)*aa_atom(3,iatom)
    Rz=aa(3,1)*aa_atom(1,iatom)+aa(3,2)*aa_atom(2,iatom)+aa(3,3)*aa_atom(3,iatom)
  END SUBROUTINE getAtomPosition_real
!-------------------------------------------------------------------------
  SUBROUTINE getAtomCenteredPositionFrom_lma(lma,j)
    use ps_nloc2_variables, only: JJ_MAP
    use aa_module, only: aa
    implicit none
    integer,intent(IN) :: lma,j
    real(8) :: d1,d2,d3
    d1=c1*JJ_MAP(1,j,lma)+JJ_MAP(4,j,lma)
    d2=c2*JJ_MAP(2,j,lma)+JJ_MAP(5,j,lma)
    d3=c3*JJ_MAP(3,j,lma)+JJ_MAP(6,j,lma)
    x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
    y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
    z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
    r = sqrt(x*x+y*y+z*z)
    return
  END SUBROUTINE getAtomCenteredPositionFrom_lma
!-------------------------------------------------------------------------
  SUBROUTINE getAtomCenteredPositionFrom_iqr(iqr,j)
    use var_para_ps_nloc_g, only: JJ_MAP_Q
    use aa_module, only: aa
    implicit none
    integer,intent(IN) :: iqr,j
    real(8) :: d1,d2,d3
    d1=c1*JJ_MAP_Q(1,j,iqr)+JJ_MAP_Q(4,j,iqr)
    d2=c2*JJ_MAP_Q(2,j,iqr)+JJ_MAP_Q(5,j,iqr)
    d3=c3*JJ_MAP_Q(3,j,iqr)+JJ_MAP_Q(6,j,iqr)
    x = aa(1,1)*d1+aa(1,2)*d2+aa(1,3)*d3-Rx
    y = aa(2,1)*d1+aa(2,2)*d2+aa(2,3)*d3-Ry
    z = aa(3,1)*d1+aa(3,2)*d2+aa(3,3)*d3-Rz
    r = sqrt(x*x+y*y+z*z)
    return
  END SUBROUTINE getAtomCenteredPositionFrom_iqr
!-------------------------------------------------------------------------
  SUBROUTINE interpolate_dviod(lm1,ir,NRc)
    use var_ps_member, only: rad1,dviod
    implicit none
    integer,intent(IN) :: lm1
    integer,intent(IN) :: ir,NRc
    real(8) :: tmp1,err0,err
    integer :: im,m1,m2
    real(8),parameter :: ep=1.d-8
#ifdef _SPLINE_
    if ( r < rad1(2,ikind) ) then
      tmp0=dviod(2,lm1,ikind)/(rad1(2,ikind)**2)
    else
      call splint(rad1(1,ikind),dviod(1,lm1,ikind),y2b(1,lm1,ikind),NRc,r,tmp0)
      tmp0=tmp0/(r*r)
    end if
#else
    if ( ir <= 2 ) then
      err0=0.d0
      tmp0=dviod(2,lm1,ikind)/(rad1(2,ikind)**2)
      if ( ir < 1 ) stop "calc_force_ps_nloc_uspp"
    else if ( ir <= NRc ) then
      err0=1.d10
      do im=1,20
        m1=max(1,ir-im)
        m2=min(ir+im,NRc)
        call polint(rad1(m1,ikind),dviod(m1,lm1,ikind),m2-m1+1,r,tmp1,err)
        if ( abs(err)<err0 ) then
          tmp0=tmp1
          err0=abs(err)
          if ( err0<ep ) exit
        end if
      end do
      tmp0=tmp0/(r*r)
    else
      write(*,*) "force_ps_nloc_uspp",ir,NRc
      stop
    end if
    maxerr=max(maxerr,err0)
#endif
  END SUBROUTINE interpolate_dviod
!-------------------------------------------------------------------------
  SUBROUTINE interpolate_dqrL(ll3,cJ,ir,NRc)
    use var_ps_member, only: rad1
    use var_ps_member_g, only: dqrL
    implicit none
    integer,intent(IN) :: ll3,cJ
    integer,intent(IN) :: ir,NRc
    real(8) :: tmp1,err0,err
    integer :: im,m1,m2
    real(8),parameter :: ep=1.d-8
#ifdef _SPLINE_
    if ( r < rad1(2,ikind) ) then
      tmp0=dqrL(2,ll3,ik2,ikind,cJ)/(rad1(2,ikind)**3)
    else
      call splint(rad1(1,ikind),dqrL(1,ll3,ik2,ikind,cJ),y2b,NRc,r,tmp0)
      tmp0=tmp0/(r*r)
    end if
#else
    if ( ir <= 2 ) then
      err0=0.d0
      tmp0=dqrL(2,ll3,ik2,ikind,cJ)/(rad1(2,ikind)**3)
      if ( ir < 1 ) stop "calc_force_ps_Q"
    else if ( ir <= NRc ) then
      err0=1.d10
      do im=1,20
        m1=max(1,ir-im)
        m2=min(ir+im,NRc)
        call polint(rad1(m1,ikind),dqrL(m1,ll3,ik2,ikind,cJ),m2-m1+1,r,tmp1,err)
        if ( abs(err)<err0 ) then
          tmp0=tmp1
          err0=abs(err)
          if ( err0<ep ) exit
        end if
      end do
      tmp0=tmp0/(r*r*r)
    else
      write(*,*) "force_ps_Q",ir,NRc
      stop
    end if
    maxerr=max(maxerr,err0)
#endif
  END SUBROUTINE interpolate_dqrL

!-------------------------------------------------------------------------
  SUBROUTINE getSHY
    use var_ps_member_g, only: max_Lref
    implicit none
    integer :: L,Lp1
    if (.not. isFirstSHY) then
      return
    endif

    isFirstSHY=.false.

    L=2*max_Lref
    Lp1=L+1
    allocate(SH_Y1(0:L,-L:L,-Lp1:Lp1,-Lp1:Lp1)) ; SH_Y1=0.d0
    allocate(SH_Y2(0:L,-L:L,-Lp1:Lp1,-Lp1:Lp1)) ; SH_Y2=0.d0
    allocate(SH_Y3(0:L,-L:L,-Lp1:Lp1,-Lp1:Lp1)) ; SH_Y3=0.d0

    SH_Y1(  0,  0,  1,  1)=   0.282094791773878d0
    SH_Y2(  0,  0,  1, -1)=   0.282094791773878d0
    SH_Y3(  0,  0,  1,  0)=   0.282094791773878d0

    if (L>=2) then
      SH_Y1(  1,  1,  0,  0)=   0.282094791773878d0
      SH_Y1(  1, -1,  2, -2)=  -0.218509686118416d0
      SH_Y1(  1,  1,  2,  0)=  -0.126156626101008d0
      SH_Y1(  1,  0,  2,  1)=   0.218509686118416d0
      SH_Y1(  1,  1,  2,  2)=   0.218509686118416d0
      SH_Y1(  2, -2,  1, -1)=  -0.218509686118416d0
      SH_Y1(  2,  1,  1,  0)=   0.218509686118416d0
      SH_Y1(  2,  0,  1,  1)=  -0.126156626101008d0
      SH_Y1(  2,  2,  1,  1)=   0.218509686118416d0
      SH_Y1(  2, -2,  3, -3)=  -0.226179013159540d0
      SH_Y1(  2, -1,  3, -2)=  -0.184674390922372d0
      SH_Y1(  2, -2,  3, -1)=   0.058399170081902d0
      SH_Y1(  2,  1,  3,  0)=  -0.143048168102669d0
      SH_Y1(  2,  0,  3,  1)=   0.202300659403421d0
      SH_Y1(  2,  2,  3,  1)=  -0.058399170081902d0
      SH_Y1(  2,  1,  3,  2)=   0.184674390922372d0
      SH_Y1(  2,  2,  3,  3)=   0.226179013159540d0

      SH_Y2(  1, -1,  0,  0)=   0.282094791773878d0
      SH_Y2(  1,  1,  2, -2)=  -0.218509686118416d0
      SH_Y2(  1,  0,  2, -1)=   0.218509686118416d0
      SH_Y2(  1, -1,  2,  0)=  -0.126156626101008d0
      SH_Y2(  1, -1,  2,  2)=  -0.218509686118416d0
      SH_Y2(  2,  0,  1, -1)=  -0.126156626101008d0
      SH_Y2(  2,  2,  1, -1)=  -0.218509686118416d0
      SH_Y2(  2, -1,  1,  0)=   0.218509686118416d0
      SH_Y2(  2, -2,  1,  1)=  -0.218509686118416d0
      SH_Y2(  2,  2,  3, -3)=   0.226179013159540d0
      SH_Y2(  2,  1,  3, -2)=  -0.184674390922372d0
      SH_Y2(  2,  0,  3, -1)=   0.202300659403421d0
      SH_Y2(  2,  2,  3, -1)=   0.058399170081902d0
      SH_Y2(  2, -1,  3,  0)=  -0.143048168102669d0
      SH_Y2(  2, -2,  3,  1)=   0.058399170081902d0
      SH_Y2(  2, -1,  3,  2)=  -0.184674390922372d0
      SH_Y2(  2, -2,  3,  3)=   0.226179013159540d0

      SH_Y3(  1,  0,  0,  0)=   0.282094791773878d0
      SH_Y3(  1, -1,  2, -1)=   0.218509686118416d0
      SH_Y3(  1,  0,  2,  0)=   0.252313252202016d0
      SH_Y3(  1,  1,  2,  1)=   0.218509686118416d0
      SH_Y3(  2, -1,  1, -1)=   0.218509686118416d0
      SH_Y3(  2,  0,  1,  0)=   0.252313252202016d0
      SH_Y3(  2,  1,  1,  1)=   0.218509686118416d0
      SH_Y3(  2, -2,  3, -2)=   0.184674390922372d0
      SH_Y3(  2, -1,  3, -1)=   0.233596680327607d0
      SH_Y3(  2,  0,  3,  0)=   0.247766695083476d0
      SH_Y3(  2,  1,  3,  1)=   0.233596680327607d0
      SH_Y3(  2,  2,  3,  2)=   0.184674390922372d0
    endif

    if (L>=4) then
      SH_Y1(  3, -3,  2, -2)=  -0.226179013159540d0
      SH_Y1(  3, -1,  2, -2)=   0.058399170081902d0
      SH_Y1(  3, -2,  2, -1)=  -0.184674390922372d0
      SH_Y1(  3,  1,  2,  0)=   0.202300659403421d0
      SH_Y1(  3,  0,  2,  1)=  -0.143048168102669d0
      SH_Y1(  3,  2,  2,  1)=   0.184674390922372d0
      SH_Y1(  3,  1,  2,  2)=  -0.058399170081902d0
      SH_Y1(  3,  3,  2,  2)=   0.226179013159540d0
      SH_Y1(  3, -3,  4, -4)=  -0.230329432980890d0
      SH_Y1(  3, -2,  4, -3)=  -0.199471140200716d0
      SH_Y1(  3, -3,  4, -2)=   0.043528171377568d0
      SH_Y1(  3, -1,  4, -2)=  -0.168583882836184d0
      SH_Y1(  3, -2,  4, -1)=   0.075393004386513d0
      SH_Y1(  3,  1,  4,  0)=  -0.150786008773027d0
      SH_Y1(  3,  0,  4,  1)=   0.194663900273006d0
      SH_Y1(  3,  2,  4,  1)=  -0.075393004386513d0
      SH_Y1(  3,  1,  4,  2)=   0.168583882836184d0
      SH_Y1(  3,  3,  4,  2)=  -0.043528171377568d0
      SH_Y1(  3,  2,  4,  3)=   0.199471140200716d0
      SH_Y1(  3,  3,  4,  4)=   0.230329432980890d0
      SH_Y1(  4, -4,  3, -3)=  -0.230329432980890d0
      SH_Y1(  4, -2,  3, -3)=   0.043528171377568d0
      SH_Y1(  4, -3,  3, -2)=  -0.199471140200716d0
      SH_Y1(  4, -1,  3, -2)=   0.075393004386513d0
      SH_Y1(  4, -2,  3, -1)=  -0.168583882836184d0
      SH_Y1(  4,  1,  3,  0)=   0.194663900273006d0
      SH_Y1(  4,  0,  3,  1)=  -0.150786008773027d0
      SH_Y1(  4,  2,  3,  1)=   0.168583882836184d0
      SH_Y1(  4,  1,  3,  2)=  -0.075393004386513d0
      SH_Y1(  4,  3,  3,  2)=   0.199471140200716d0
      SH_Y1(  4,  2,  3,  3)=  -0.043528171377568d0
      SH_Y1(  4,  4,  3,  3)=   0.230329432980890d0
      SH_Y1(  4, -4,  5, -5)=  -0.232932108055429d0
      SH_Y1(  4, -3,  5, -4)=  -0.208340811101706d0
      SH_Y1(  4, -4,  5, -3)=   0.034723468516951d0
      SH_Y1(  4, -2,  5, -3)=  -0.183739324706867d0
      SH_Y1(  4, -3,  5, -2)=   0.060142811686378d0
      SH_Y1(  4, -1,  5, -2)=  -0.159122922870344d0
      SH_Y1(  4, -2,  5, -1)=   0.085054779966126d0
      SH_Y1(  4,  1,  5,  0)=  -0.155288072036953d0
      SH_Y1(  4,  0,  5,  1)=   0.190188269815546d0
      SH_Y1(  4,  2,  5,  1)=  -0.085054779966126d0
      SH_Y1(  4,  1,  5,  2)=   0.159122922870344d0
      SH_Y1(  4,  3,  5,  2)=  -0.060142811686378d0
      SH_Y1(  4,  2,  5,  3)=   0.183739324706867d0
      SH_Y1(  4,  4,  5,  3)=  -0.034723468516951d0
      SH_Y1(  4,  3,  5,  4)=   0.208340811101706d0
      SH_Y1(  4,  4,  5,  5)=   0.232932108055429d0

      SH_Y2(  3,  1,  2, -2)=   0.058399170081902d0
      SH_Y2(  3,  3,  2, -2)=   0.226179013159540d0
      SH_Y2(  3,  0,  2, -1)=  -0.143048168102669d0
      SH_Y2(  3,  2,  2, -1)=  -0.184674390922372d0
      SH_Y2(  3, -1,  2,  0)=   0.202300659403421d0
      SH_Y2(  3, -2,  2,  1)=  -0.184674390922372d0
      SH_Y2(  3, -3,  2,  2)=   0.226179013159540d0
      SH_Y2(  3, -1,  2,  2)=   0.058399170081902d0
      SH_Y2(  3,  3,  4, -4)=  -0.230329432980890d0
      SH_Y2(  3,  2,  4, -3)=   0.199471140200716d0
      SH_Y2(  3,  1,  4, -2)=  -0.168583882836184d0
      SH_Y2(  3,  3,  4, -2)=  -0.043528171377568d0
      SH_Y2(  3,  0,  4, -1)=   0.194663900273006d0
      SH_Y2(  3,  2,  4, -1)=   0.075393004386513d0
      SH_Y2(  3, -1,  4,  0)=  -0.150786008773027d0
      SH_Y2(  3, -2,  4,  1)=   0.075393004386513d0
      SH_Y2(  3, -3,  4,  2)=  -0.043528171377568d0
      SH_Y2(  3, -1,  4,  2)=  -0.168583882836184d0
      SH_Y2(  3, -2,  4,  3)=   0.199471140200716d0
      SH_Y2(  3, -3,  4,  4)=  -0.230329432980890d0
      SH_Y2(  4,  2,  3, -3)=  -0.043528171377568d0
      SH_Y2(  4,  4,  3, -3)=  -0.230329432980890d0
      SH_Y2(  4,  1,  3, -2)=   0.075393004386513d0
      SH_Y2(  4,  3,  3, -2)=   0.199471140200716d0
      SH_Y2(  4,  0,  3, -1)=  -0.150786008773027d0
      SH_Y2(  4,  2,  3, -1)=  -0.168583882836184d0
      SH_Y2(  4, -1,  3,  0)=   0.194663900273006d0
      SH_Y2(  4, -2,  3,  1)=  -0.168583882836184d0
      SH_Y2(  4, -3,  3,  2)=   0.199471140200716d0
      SH_Y2(  4, -1,  3,  2)=   0.075393004386513d0
      SH_Y2(  4, -4,  3,  3)=  -0.230329432980890d0
      SH_Y2(  4, -2,  3,  3)=  -0.043528171377568d0
      SH_Y2(  4,  4,  5, -5)=   0.232932108055429d0
      SH_Y2(  4,  3,  5, -4)=  -0.208340811101706d0
      SH_Y2(  4,  2,  5, -3)=   0.183739324706867d0
      SH_Y2(  4,  4,  5, -3)=   0.034723468516951d0
      SH_Y2(  4,  1,  5, -2)=  -0.159122922870344d0
      SH_Y2(  4,  3,  5, -2)=  -0.060142811686378d0
      SH_Y2(  4,  0,  5, -1)=   0.190188269815546d0
      SH_Y2(  4,  2,  5, -1)=   0.085054779966126d0
      SH_Y2(  4, -1,  5,  0)=  -0.155288072036953d0
      SH_Y2(  4, -2,  5,  1)=   0.085054779966126d0
      SH_Y2(  4, -3,  5,  2)=  -0.060142811686378d0
      SH_Y2(  4, -1,  5,  2)=  -0.159122922870344d0
      SH_Y2(  4, -4,  5,  3)=   0.034723468516951d0
      SH_Y2(  4, -2,  5,  3)=   0.183739324706867d0
      SH_Y2(  4, -3,  5,  4)=  -0.208340811101706d0
      SH_Y2(  4, -4,  5,  5)=   0.232932108055429d0

      SH_Y3(  3, -2,  2, -2)=   0.184674390922372d0
      SH_Y3(  3, -1,  2, -1)=   0.233596680327607d0
      SH_Y3(  3,  0,  2,  0)=   0.247766695083476d0
      SH_Y3(  3,  1,  2,  1)=   0.233596680327607d0
      SH_Y3(  3,  2,  2,  2)=   0.184674390922372d0
      SH_Y3(  3, -3,  4, -3)=   0.162867503967640d0
      SH_Y3(  3, -2,  4, -2)=   0.213243618622923d0
      SH_Y3(  3, -1,  4, -1)=   0.238413613504448d0
      SH_Y3(  3,  0,  4,  0)=   0.246232521229829d0
      SH_Y3(  3,  1,  4,  1)=   0.238413613504448d0
      SH_Y3(  3,  2,  4,  2)=   0.213243618622923d0
      SH_Y3(  3,  3,  4,  3)=   0.162867503967640d0
      SH_Y3(  4, -3,  3, -3)=   0.162867503967640d0
      SH_Y3(  4, -2,  3, -2)=   0.213243618622923d0
      SH_Y3(  4, -1,  3, -1)=   0.238413613504448d0
      SH_Y3(  4,  0,  3,  0)=   0.246232521229829d0
      SH_Y3(  4,  1,  3,  1)=   0.238413613504448d0
      SH_Y3(  4,  2,  3,  2)=   0.213243618622923d0
      SH_Y3(  4,  3,  3,  3)=   0.162867503967640d0
      SH_Y3(  4, -4,  5, -4)=   0.147319200327922d0
      SH_Y3(  4, -3,  5, -3)=   0.196425600437230d0
      SH_Y3(  4, -2,  5, -2)=   0.225033795607689d0
      SH_Y3(  4, -1,  5, -1)=   0.240571246745510d0
      SH_Y3(  4,  0,  5,  0)=   0.245532000546537d0
      SH_Y3(  4,  1,  5,  1)=   0.240571246745510d0
      SH_Y3(  4,  2,  5,  2)=   0.225033795607689d0
      SH_Y3(  4,  3,  5,  3)=   0.196425600437230d0
      SH_Y3(  4,  4,  5,  4)=   0.147319200327922d0
    endif

    return
  END SUBROUTINE getSHY


  SUBROUTINE getCijLM
    implicit none
    if (.not. isFirstC) then
      return
    endif

    isFirstC=.false.

    allocate(C_ijLM(0:4,-4:4,0:2,-2:2,0:2,-2:2))
    C_ijLM=0.d0

    C_ijLM(  0,  0,  0,  0,  0,  0  )=    0.282094791773878d0
    C_ijLM(  1, -1,  1, -1,  0,  0  )=    0.282094791773878d0
    C_ijLM(  1,  0,  1,  0,  0,  0  )=    0.282094791773878d0
    C_ijLM(  1,  1,  1,  1,  0,  0  )=    0.282094791773878d0
    C_ijLM(  2, -2,  2, -2,  0,  0  )=    0.282094791773878d0
    C_ijLM(  2, -1,  2, -1,  0,  0  )=    0.282094791773878d0
    C_ijLM(  2,  0,  2,  0,  0,  0  )=    0.282094791773878d0
    C_ijLM(  2,  1,  2,  1,  0,  0  )=    0.282094791773878d0
    C_ijLM(  2,  2,  2,  2,  0,  0  )=    0.282094791773878d0
    C_ijLM(  1, -1,  0,  0,  1, -1  )=    0.282094791773878d0
    C_ijLM(  0,  0,  1, -1,  1, -1  )=    0.282094791773878d0
    C_ijLM(  2, -2,  1,  1,  1, -1  )=   -0.218509686118416d0
    C_ijLM(  2, -1,  1,  0,  1, -1  )=    0.218509686118416d0
    C_ijLM(  2,  0,  1, -1,  1, -1  )=   -0.126156626101008d0
    C_ijLM(  2,  2,  1, -1,  1, -1  )=   -0.218509686118416d0
    C_ijLM(  1, -1,  2,  0,  1, -1  )=   -0.126156626101008d0
    C_ijLM(  1, -1,  2,  2,  1, -1  )=   -0.218509686118416d0
    C_ijLM(  1,  0,  2, -1,  1, -1  )=    0.218509686118416d0
    C_ijLM(  1,  1,  2, -2,  1, -1  )=   -0.218509686118416d0
    C_ijLM(  3, -3,  2,  2,  1, -1  )=    0.226179013159540d0
    C_ijLM(  3, -2,  2,  1,  1, -1  )=   -0.184674390922372d0
    C_ijLM(  3, -1,  2,  0,  1, -1  )=    0.202300659403421d0
    C_ijLM(  3, -1,  2,  2,  1, -1  )=    0.058399170081902d0
    C_ijLM(  3,  0,  2, -1,  1, -1  )=   -0.143048168102669d0
    C_ijLM(  3,  1,  2, -2,  1, -1  )=    0.058399170081902d0
    C_ijLM(  3,  2,  2, -1,  1, -1  )=   -0.184674390922372d0
    C_ijLM(  3,  3,  2, -2,  1, -1  )=    0.226179013159540d0
    C_ijLM(  1,  0,  0,  0,  1,  0  )=    0.282094791773878d0
    C_ijLM(  0,  0,  1,  0,  1,  0  )=    0.282094791773878d0
    C_ijLM(  2, -1,  1, -1,  1,  0  )=    0.218509686118416d0
    C_ijLM(  2,  0,  1,  0,  1,  0  )=    0.252313252202016d0
    C_ijLM(  2,  1,  1,  1,  1,  0  )=    0.218509686118416d0
    C_ijLM(  1, -1,  2, -1,  1,  0  )=    0.218509686118416d0
    C_ijLM(  1,  0,  2,  0,  1,  0  )=    0.252313252202016d0
    C_ijLM(  1,  1,  2,  1,  1,  0  )=    0.218509686118416d0
    C_ijLM(  3, -2,  2, -2,  1,  0  )=    0.184674390922372d0
    C_ijLM(  3, -1,  2, -1,  1,  0  )=    0.233596680327607d0
    C_ijLM(  3,  0,  2,  0,  1,  0  )=    0.247766695083476d0
    C_ijLM(  3,  1,  2,  1,  1,  0  )=    0.233596680327607d0
    C_ijLM(  3,  2,  2,  2,  1,  0  )=    0.184674390922372d0
    C_ijLM(  1,  1,  0,  0,  1,  1  )=    0.282094791773878d0
    C_ijLM(  0,  0,  1,  1,  1,  1  )=    0.282094791773878d0
    C_ijLM(  2, -2,  1, -1,  1,  1  )=   -0.218509686118416d0
    C_ijLM(  2,  0,  1,  1,  1,  1  )=   -0.126156626101008d0
    C_ijLM(  2,  1,  1,  0,  1,  1  )=    0.218509686118416d0
    C_ijLM(  2,  2,  1,  1,  1,  1  )=    0.218509686118416d0
    C_ijLM(  1, -1,  2, -2,  1,  1  )=   -0.218509686118416d0
    C_ijLM(  1,  0,  2,  1,  1,  1  )=    0.218509686118416d0
    C_ijLM(  1,  1,  2,  0,  1,  1  )=   -0.126156626101008d0
    C_ijLM(  1,  1,  2,  2,  1,  1  )=    0.218509686118416d0
    C_ijLM(  3, -3,  2, -2,  1,  1  )=   -0.226179013159540d0
    C_ijLM(  3, -2,  2, -1,  1,  1  )=   -0.184674390922372d0
    C_ijLM(  3, -1,  2, -2,  1,  1  )=    0.058399170081902d0
    C_ijLM(  3,  0,  2,  1,  1,  1  )=   -0.143048168102669d0
    C_ijLM(  3,  1,  2,  0,  1,  1  )=    0.202300659403421d0
    C_ijLM(  3,  1,  2,  2,  1,  1  )=   -0.058399170081902d0
    C_ijLM(  3,  2,  2,  1,  1,  1  )=    0.184674390922372d0
    C_ijLM(  3,  3,  2,  2,  1,  1  )=    0.226179013159540d0
    C_ijLM(  2, -2,  0,  0,  2, -2  )=    0.282094791773878d0
    C_ijLM(  1, -1,  1,  1,  2, -2  )=   -0.218509686118416d0
    C_ijLM(  1,  1,  1, -1,  2, -2  )=   -0.218509686118416d0
    C_ijLM(  3, -3,  1,  1,  2, -2  )=   -0.226179013159540d0
    C_ijLM(  3, -2,  1,  0,  2, -2  )=    0.184674390922372d0
    C_ijLM(  3, -1,  1,  1,  2, -2  )=    0.058399170081902d0
    C_ijLM(  3,  1,  1, -1,  2, -2  )=    0.058399170081902d0
    C_ijLM(  3,  3,  1, -1,  2, -2  )=    0.226179013159540d0
    C_ijLM(  0,  0,  2, -2,  2, -2  )=    0.282094791773878d0
    C_ijLM(  2, -2,  2,  0,  2, -2  )=   -0.180223751572869d0
    C_ijLM(  2, -1,  2,  1,  2, -2  )=   -0.156078347227440d0
    C_ijLM(  2,  0,  2, -2,  2, -2  )=   -0.180223751572869d0
    C_ijLM(  2,  1,  2, -1,  2, -2  )=   -0.156078347227440d0
    C_ijLM(  4, -4,  2,  2,  2, -2  )=    0.238413613504448d0
    C_ijLM(  4, -3,  2,  1,  2, -2  )=   -0.168583882836184d0
    C_ijLM(  4, -2,  2,  0,  2, -2  )=    0.156078347227440d0
    C_ijLM(  4, -1,  2,  1,  2, -2  )=    0.063718718434028d0
    C_ijLM(  4,  0,  2, -2,  2, -2  )=    0.040299255967697d0
    C_ijLM(  4,  1,  2, -1,  2, -2  )=    0.063718718434028d0
    C_ijLM(  4,  3,  2, -1,  2, -2  )=    0.168583882836184d0
    C_ijLM(  4,  4,  2, -2,  2, -2  )=   -0.238413613504448d0
    C_ijLM(  2, -1,  0,  0,  2, -1  )=    0.282094791773878d0
    C_ijLM(  1, -1,  1,  0,  2, -1  )=    0.218509686118416d0
    C_ijLM(  1,  0,  1, -1,  2, -1  )=    0.218509686118416d0
    C_ijLM(  3, -2,  1,  1,  2, -1  )=   -0.184674390922372d0
    C_ijLM(  3, -1,  1,  0,  2, -1  )=    0.233596680327607d0
    C_ijLM(  3,  0,  1, -1,  2, -1  )=   -0.143048168102669d0
    C_ijLM(  3,  2,  1, -1,  2, -1  )=   -0.184674390922372d0
    C_ijLM(  0,  0,  2, -1,  2, -1  )=    0.282094791773878d0
    C_ijLM(  2, -2,  2,  1,  2, -1  )=   -0.156078347227440d0
    C_ijLM(  2, -1,  2,  0,  2, -1  )=    0.090111875786434d0
    C_ijLM(  2, -1,  2,  2,  2, -1  )=   -0.156078347227440d0
    C_ijLM(  2,  0,  2, -1,  2, -1  )=    0.090111875786434d0
    C_ijLM(  2,  1,  2, -2,  2, -1  )=   -0.156078347227440d0
    C_ijLM(  2,  2,  2, -1,  2, -1  )=   -0.156078347227440d0
    C_ijLM(  4, -3,  2,  2,  2, -1  )=    0.168583882836184d0
    C_ijLM(  4, -2,  2,  1,  2, -1  )=   -0.180223751572869d0
    C_ijLM(  4, -1,  2,  0,  2, -1  )=    0.220728115441823d0
    C_ijLM(  4, -1,  2,  2,  2, -1  )=    0.063718718434028d0
    C_ijLM(  4,  0,  2, -1,  2, -1  )=   -0.161197023870787d0
    C_ijLM(  4,  1,  2, -2,  2, -1  )=    0.063718718434028d0
    C_ijLM(  4,  2,  2, -1,  2, -1  )=   -0.180223751572869d0
    C_ijLM(  4,  3,  2, -2,  2, -1  )=    0.168583882836184d0
    C_ijLM(  2,  0,  0,  0,  2,  0  )=    0.282094791773878d0
    C_ijLM(  1, -1,  1, -1,  2,  0  )=   -0.126156626101008d0
    C_ijLM(  1,  0,  1,  0,  2,  0  )=    0.252313252202016d0
    C_ijLM(  1,  1,  1,  1,  2,  0  )=   -0.126156626101008d0
    C_ijLM(  3, -1,  1, -1,  2,  0  )=    0.202300659403421d0
    C_ijLM(  3,  0,  1,  0,  2,  0  )=    0.247766695083476d0
    C_ijLM(  3,  1,  1,  1,  2,  0  )=    0.202300659403421d0
    C_ijLM(  0,  0,  2,  0,  2,  0  )=    0.282094791773878d0
    C_ijLM(  2, -2,  2, -2,  2,  0  )=   -0.180223751572869d0
    C_ijLM(  2, -1,  2, -1,  2,  0  )=    0.090111875786434d0
    C_ijLM(  2,  0,  2,  0,  2,  0  )=    0.180223751572869d0
    C_ijLM(  2,  1,  2,  1,  2,  0  )=    0.090111875786434d0
    C_ijLM(  2,  2,  2,  2,  2,  0  )=   -0.180223751572869d0
    C_ijLM(  4, -2,  2, -2,  2,  0  )=    0.156078347227440d0
    C_ijLM(  4, -1,  2, -1,  2,  0  )=    0.220728115441823d0
    C_ijLM(  4,  0,  2,  0,  2,  0  )=    0.241795535806181d0
    C_ijLM(  4,  1,  2,  1,  2,  0  )=    0.220728115441823d0
    C_ijLM(  4,  2,  2,  2,  2,  0  )=    0.156078347227440d0
    C_ijLM(  2,  1,  0,  0,  2,  1  )=    0.282094791773878d0
    C_ijLM(  1,  0,  1,  1,  2,  1  )=    0.218509686118416d0
    C_ijLM(  1,  1,  1,  0,  2,  1  )=    0.218509686118416d0
    C_ijLM(  3, -2,  1, -1,  2,  1  )=   -0.184674390922372d0
    C_ijLM(  3,  0,  1,  1,  2,  1  )=   -0.143048168102669d0
    C_ijLM(  3,  1,  1,  0,  2,  1  )=    0.233596680327607d0
    C_ijLM(  3,  2,  1,  1,  2,  1  )=    0.184674390922372d0
    C_ijLM(  0,  0,  2,  1,  2,  1  )=    0.282094791773878d0
    C_ijLM(  2, -2,  2, -1,  2,  1  )=   -0.156078347227440d0
    C_ijLM(  2, -1,  2, -2,  2,  1  )=   -0.156078347227440d0
    C_ijLM(  2,  0,  2,  1,  2,  1  )=    0.090111875786434d0
    C_ijLM(  2,  1,  2,  0,  2,  1  )=    0.090111875786434d0
    C_ijLM(  2,  1,  2,  2,  2,  1  )=    0.156078347227440d0
    C_ijLM(  2,  2,  2,  1,  2,  1  )=    0.156078347227440d0
    C_ijLM(  4, -3,  2, -2,  2,  1  )=   -0.168583882836184d0
    C_ijLM(  4, -2,  2, -1,  2,  1  )=   -0.180223751572869d0
    C_ijLM(  4, -1,  2, -2,  2,  1  )=    0.063718718434028d0
    C_ijLM(  4,  0,  2,  1,  2,  1  )=   -0.161197023870787d0
    C_ijLM(  4,  1,  2,  0,  2,  1  )=    0.220728115441823d0
    C_ijLM(  4,  1,  2,  2,  2,  1  )=   -0.063718718434028d0
    C_ijLM(  4,  2,  2,  1,  2,  1  )=    0.180223751572869d0
    C_ijLM(  4,  3,  2,  2,  2,  1  )=    0.168583882836184d0
    C_ijLM(  2,  2,  0,  0,  2,  2  )=    0.282094791773878d0
    C_ijLM(  1, -1,  1, -1,  2,  2  )=   -0.218509686118416d0
    C_ijLM(  1,  1,  1,  1,  2,  2  )=    0.218509686118416d0
    C_ijLM(  3, -3,  1, -1,  2,  2  )=    0.226179013159540d0
    C_ijLM(  3, -1,  1, -1,  2,  2  )=    0.058399170081902d0
    C_ijLM(  3,  1,  1,  1,  2,  2  )=   -0.058399170081902d0
    C_ijLM(  3,  2,  1,  0,  2,  2  )=    0.184674390922372d0
    C_ijLM(  3,  3,  1,  1,  2,  2  )=    0.226179013159540d0
    C_ijLM(  0,  0,  2,  2,  2,  2  )=    0.282094791773878d0
    C_ijLM(  2, -1,  2, -1,  2,  2  )=   -0.156078347227440d0
    C_ijLM(  2,  0,  2,  2,  2,  2  )=   -0.180223751572869d0
    C_ijLM(  2,  1,  2,  1,  2,  2  )=    0.156078347227440d0
    C_ijLM(  2,  2,  2,  0,  2,  2  )=   -0.180223751572869d0
    C_ijLM(  4, -4,  2, -2,  2,  2  )=    0.238413613504448d0
    C_ijLM(  4, -3,  2, -1,  2,  2  )=    0.168583882836184d0
    C_ijLM(  4, -1,  2, -1,  2,  2  )=    0.063718718434028d0
    C_ijLM(  4,  0,  2,  2,  2,  2  )=    0.040299255967697d0
    C_ijLM(  4,  1,  2,  1,  2,  2  )=   -0.063718718434028d0
    C_ijLM(  4,  2,  2,  0,  2,  2  )=    0.156078347227440d0
    C_ijLM(  4,  3,  2,  1,  2,  2  )=    0.168583882836184d0
    C_ijLM(  4,  4,  2,  2,  2,  2  )=    0.238413613504448d0

    return
  End SUBROUTINE getCijLM

END MODULE force_sub_module
