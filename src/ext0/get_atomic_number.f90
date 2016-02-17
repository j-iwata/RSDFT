SUBROUTINE get_atomic_number( element_name, z )

  implicit none
  character(*),intent(IN) :: element_name
  integer,intent(OUT) :: z
  character(2),save :: name(112)
  integer :: i

  data name/ "H" ,"He", &
             "Li","Be","B" ,"C" ,"N" ,"O" ,"F" ,"Ne", &
             "Na","Mg","Al","Si","P" ,"S" ,"Cl","Ar", &
             "K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co", &
             "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", &
             "Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh", &
             "Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe", &
             "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu", &
             "Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", &
             "Hf","Ta","W" ,"Re","Os","Ir", &
             "Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", &
             "Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am", &
             "Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf", &
             "Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn" /

  z=0

  do i=1,size(name)
     if ( element_name == name(i) ) then
        z=i
        exit
     end if
  end do

END SUBROUTINE get_atomic_number

