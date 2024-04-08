program test_proj
   use iso_c_binding
  !!use symproj
  integer(C_INT) :: num_kpoints=1
  integer(C_INT) :: order=96
  character(len=11), parameter :: filename="s123_3_proj"
  double complex, allocatable:: dmat(:,:)
  integer :: ii,jj
  allocate(dmat(1:order,1:order))
  dmat=(0,0)
  call parse_SYMPROJ("symproj.bin", %VAL(num_kpoints), %VAL(order), dmat)
  !!write(*,*) "print dmat"
  
  !!dmat=transpose(dmat)
  do ii=1,order
     do jj=1, order

  write (*,*) dmat(ii,jj)
  end do
  end do

  deallocate(dmat)
       
end program test_proj
