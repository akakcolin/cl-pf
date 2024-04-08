module permu 
  use accuracy
  use constants

  implicit none
  private
  public :: sym_permu

contains
  
  ! determines the loop-structure of the element IN
  ! considered as the permutation of G elements
  ! numl = number of loops
  ! loopl(K) = length of Kth loop
  ! lpstr(K) = loopstructure, parentheses omitted

  !The loop of each group element, i.e. the successive powers of each group element is
  !calculated. The loop length of each loop is stored in loopl(1:numl), where numl is the
  !number of loops. The group elements are ordered in the vector lpstr(1:G) as successive
  !loops.
  subroutine sym_permu(loopl, lpstr, numl, multab, G, inel, IN, steer)
    integer, intent(out) :: loopl(:)
    integer, intent(out) :: lpstr(:)
    integer, intent(out) :: numl

    integer, intent(in) :: multab(:,:)
    integer, intent(in) :: G 
    integer, intent(in) :: inel(:)
    integer, intent(in) :: IN
    integer, intent(in) :: steer(:)

    integer :: K
    integer :: L1, LT, K1, K2, N
    integer, allocatable :: flip(:)
    
    allocate(flip(G))
   
    flip(:) = 0
    K = inel(IN)
    K2 = 1
    numl = 0
    L1 = 1
    LT = 1
    
    do K1=1, G
       if(flip(K1) .ne. 1) then
          flip(K1) = 1
          K2 = K1
          lpstr(LT) = K2
          N = multab(K, K2)
          do while(N .ne. K1)
             LT = LT + 1
             L1 = L1 + 1
             lpstr(LT) = N
             flip(N) = 1
             K2 = N
             N = multab(K, K2)
             write(*,*) multab(K,K2)
          end do
          numl = numl + 1
          loopl(numl) = L1
          if (LT .eq. G) then
             exit
          end if
          LT = LT + 1
          L1 = 1
       end if
    end do
    if (steer(4) .ne. 0) then
        write(*,*) "Loopstructure of group element ", IN
        write(*,*) "Loopstructure"
        write(*,*) lpstr(1:G)
        write(*,*) "Looplength"
        write(*,*) loopl(1:numl)
        write(*,*)"numl", numl
    end if
    deallocate(flip)
  end subroutine sym_permu
end module permu
