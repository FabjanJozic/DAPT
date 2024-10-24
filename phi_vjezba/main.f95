! A fortran95 program for G95
! By WQY
program main
  real phi(50), x, pi
  integer i
  parameter( pi=3.1415927 )
  x = -3.0
  i=1
  do while (x .le. 3.2)
    phi(i) = 1.0/SQRT(2.0*pi)*EXP(-x**2.0/2.0)
    write(*,'("For x = ",F4.1,", phi(i) = ",F6.4)') x, phi(i)
    x = x + 0.2
  enddo
end
