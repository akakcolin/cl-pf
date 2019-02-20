subroutine print_coords(coords) bind(C, name="print_coords")
  implicit none
  ! indexes are reversed compared to C code
  double precision, dimension(3,2) :: coords
  integer :: i, j

  write (*,*) "Called print_coords()"
  write (*,*) "SIZE:", size(coords, dim=1), size(coords, dim=2)
  do i = 1, 2
     do j = 1, 3
        ! indexes are reversed compared to C code
        write (*,*) "COORDS()", i, j, coords(j, i)
     enddo
  enddo
  
end subroutine print_coords
