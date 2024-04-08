module symproj
  interface
     
     function parse_SYMPROJ_real(fname, num_kpoints, dftb_matrix_order, dmatrixs) BIND(C, name='parse_SYMPROJ_real_')
       use iso_c_binding
       
       type (C_PTR), value :: fname
       integer(C_INT), value :: num_kpoints
       integer(C_INT), value :: dftb_matrix_order
       type(C_PTR), value :: dmatrixs
     end function parse_SYMPROJ_real
  end interface
end module symproj
