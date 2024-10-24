! A fortran95 program for G95
! By WQY
program main

  real r, area, pi
  write (*,*) "Give radius r:"
  read (*,*) r
  pi = atan(1.0e0)*4.0e0
  area = pi*r*r
  write (*,*) "Area = ", area

end
