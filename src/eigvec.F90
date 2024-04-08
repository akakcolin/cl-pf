module eigvec 
  use accuracy
  use constants
  use permu
  implicit none
  private
  public :: sym_eigvec 
contains

  ! calculates the eigenvectors of group element, corresponding th eigenvalue eval and projects these
  ! eigenvectors on the Jth irreducible subspace. The eigenvectors are stored in fi(1:G, 1:nvec)

  !The eigenvectors of group element IN, with eigenvalue lab are calculated. Function
  !permu is called to create eigenvectors, using the loop structure of element IN. The
  !eigenvectors are projected on the Jth irreducible eigenspace, using the projection
  !operator Sj. The resulting eigenvectors are orthonormalized to each other. If no
  !eigenvector corresponding to eigenvalue lab is found, nvec is set equal to 0.
  
  subroutine sym_eigvec(fi, nvec, elem, eval, J, kloop, inel, cind, ch, multab, G, steer)

    complex(dp), intent(out) :: fi(:,:)
    integer, intent(out) :: nvec

    integer, intent(in) :: elem
    complex(dp), intent(in) :: eval !! non-degenerate eigenvalue 
    integer, intent(in) :: J
    integer, intent(in) :: kloop
    integer, intent(in) :: inel(:)
    integer, intent(in) :: cind(:)
    complex(dp), intent(in) :: ch(:,:)
    integer, intent(in) :: multab(:,:)
    integer, intent(in) :: G
    integer, intent(in) :: steer(:)


    integer :: II, I3 
    integer :: K2, K3, K4, K5, K6
    integer :: nvr, IND
    integer :: nml, numl
    integer ::  LPE
    integer ::  N
    real(dp) :: rnorm
    complex(dp) :: P

    integer, allocatable :: loopl(:)
    integer, allocatable :: lpstr(:)
    real(dp), allocatable :: vec(:)
    real(dp), allocatable :: vec2(:)

    allocate(loopl(G))
    allocate(lpstr(G))
    allocate(vec(G))
    allocate(vec2(G))

    if(kloop .ne. 1) then
       call sym_permu(loopl, lpstr, numl, multab, G, inel, elem, steer)
    end if
    vec(:) = 0

    IND = 1
    nvr = 0
    ! numl  is the number of loops
    ! lpstr loopstructure
    ! loopl length of Kth loop
    do nml= 1, numl
       if(nml >1) then
          IND = IND + loopl(nml-1)
       end if

       nvr = nvr + 1
       vec(:) = 0.0
       K2 = lpstr(IND)
       vec(K2) = 1
       P = 1
       P = P*eval
       LPE = IND + loopl(nml)
       do N = IND+1, LPE-1
          K2 = lpstr(N)
          vec(K2) = P
          P = P*eval
       end do
       write(*,*) "degen create " , nvr, " eigenvector ", vec
       ! Projection SJ*vec
       do K4 = 1, G
          fi(K4, nvr) = 0
       end do
       do K4 = 1, G
          do K5 = 1, G
             K6 = inel(K5)
             K6 = multab(K4, K6)
             K6 = cind(K6)
             fi(K4, nvr) = fi(K4, nvr) + conjg(ch(J,K6))*vec(K5)
          end do
       end do
    end do

    ! the first non-zero eigencolum is normalized
    nvec = nvr
    nvr = 1

    do while(nvr <= nvec)
       ! in case thera are no eigenvectors with eigenvalue eval, belonging to the Jth irreducible
       ! subspace, set nvec= 0 and return
       II = 1
       do while(II <= G)
          if(abs(fi(II, nvr)) >= 0.001) then
             exit
          end if
          II = II + 1
       end do
       if( II <= G) then
          !rnorm = norm2(fi(1:G, nvr))
          do I3 = 1, G
             rnorm = rnorm + fi(I3, nvr)*fi(I3, nvr)
          end do
          !rnorm = sum(abs(fi(1:G, nvr))*abs(fi(1:G,nvr)))
          rnorm = 1/ sqrt(rnorm)

          fi(1:G, nvr) = fi(1:G, nvr)*rnorm
          if(nvr .gt. 1) then
             ! orthogonalize the eigencolumns of elements elem with eigenvalue eval to each other
             K3 = nvr - 1
             vec(1:K3) = matmul(transpose(fi(1:G, 1:K3)), fi(1:G, nvr))

             do I3 = 1, G
                vec2(I3) = fi(I3, nvr)
                do K2=1,K3
                   vec2(I3) = vec2(I3) - vec(K2)*fi(I3, K2)
                end do
             end do
             I3 = 1
             do while (I3 <= G)
                if(abs(vec2(I3)) >= 0.001) then
                   exit
                end if
                I3 = I3 + 1
             end do
             if(I3 <= G) then
                rnorm = norm2(vec2(1:G))
                rnorm = 1/rnorm
                fi(1:G, nvr) = vec2(1:G) * rnorm
                nvr = nvr + 1
             else
                nvec = nvec + 1
                fi(1:G, nvr:nvec) = fi(1:G, (nvr+1):(nvec+1))
             end if
          else
             nvr = nvr + 1
          end if
       else
          nvec = nvec - 1
          fi(1:G, nvr:nvec) = fi(1:G,(nvr+1):(nvec+1))
       end if
    end do
    !
    deallocate(loopl)
    deallocate(lpstr)
    deallocate(vec)
    deallocate(vec2)
  end subroutine sym_eigvec

end module eigvec
