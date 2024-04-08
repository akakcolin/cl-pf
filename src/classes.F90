module classes
  use accuracy
  use constants
  implicit none
  private
  public :: sym_classes

contains
  ! The group element indices of class I are put into classl(nfirst(I)) to
  ! classl(nfirst(I) + h(I) - 1). Within the class, these numbers are in
  ! increasing order. h(I) is the order of class I. ncl is the number of
  ! classes for the group.
  subroutine sym_classes( nfirst, h, classl, ncl, G, mtab2, inel)
    integer, intent(out) :: nfirst(:)
    integer, intent(out) :: h(:)
    integer, intent(out) :: classl(:)
    integer, intent(inout) :: ncl

    integer, intent(in) :: G
    integer, intent(in) :: mtab2(:,:)
    integer, intent(in) :: inel(:)

    integer :: L, M, N1, N2, I, NT
    integer :: I1, J, L1, M1, N, NIP
    integer, allocatable :: KL(:)
    allocate(KL(G))
    h(1) = 1
    classl(1) = 1
    ncl = 1
    nfirst(1) = 1

    if(G > 1) then
       NT = 2
       L  = 2
       ncl = 2
       do while (L <= G)
          h(ncl) = 1
          N2 = L
          classl(NT) = L
          nfirst(ncl) = NT

          do M = 1, G
             L1 = mtab2(M, L)
             M1 = inel(M)
             KL(M) = mtab2(L1, M1)
          end do

          do M = 1, G
             do I = 1, G
                if (KL(I) == N2) then
                   KL(I) = 0
                end if
             end do

             if (KL(M) .ne. 0) then
                NT = NT + 1
                h(ncl) = h(ncl) + 1
                ! a new class element is registered , The other KL(I) equal to this one
                ! are skipped
                N2 = KL(M)
                classl(NT) = N2
             end if
          end do

          ! NT is equal to the total number of elements that have been registered in class so far
          N = nfirst(ncl)
          if( N <= NT) then
             do I = N, NT
                I1 = I + 1
                if( I1 <= NT) then
                   do J = I1, NT
                      if( classl(I) > classl(J)) then
                         NIP = classl(J)
                         classl(J) = classl(I)
                         classl(I) = NIP
                      end if
                   end do
                end if
             end do
          end if

          L = L + 1
          I = 2
          do while( I <= NT)
             if (L == classl(I)) then
                L = L+ 1
                I = 1
             end if
             I = I + 1
          end do
          if (L <= G) then
             ncl = ncl + 1
             NT = NT + 1
          end if
       end do
    end if
    deallocate(KL)
  end subroutine sym_classes

end module classes

