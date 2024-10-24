! A fortran95 program for G95
! By WQY
program main
  real r(1001), F(1000)
  integer j
  call force(r, F)
  open (10, file='sila.dat', status='new')
  do j=1, 1000
    write (10,'(F9.5," "F40.5)') r(j), F(j)
  end do
  close (10)
end


subroutine force(r,F)
    implicit none
    real ep, sig, dr
    integer i
    real r(1001), U(1001), F(1000)
    ep = 0.010323
    sig = 3.405
    dr = 0.01
    r(1) = dr
    U(1) = 4*ep*((sig/r(1))**12-(sig/r(1))**6)
    i = 2
    do while (i.le.1001)
        r(i) = dr*i
        U(i) = 4*ep*((sig/r(i))**12-(sig/r(i))**6)
        F(i-1) = -(U(i)-U(i-1))/dr
        i = i+1
    end do
end subroutine
