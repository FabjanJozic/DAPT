! A fortran95 program for G95
! By WQY
program main
  real t(6001), x(6001), v(6001)
  integer j
  call kutija(t, x, v)
  open (10, file='box.dat', status='new')
  write (10, '(A5," "A6," "A5)') "t", "x(t)", "v(t)"
  do j=1, 6001
    write (10,'(F5.2," "F6.2," "F5.2)') t(j), x(j), v(j)
  end do
  close (10)
end

subroutine kutija(t, x, v)
    implicit none
    real t(6001), x(6001), v(6001)
    real x0, v0, L, t0, dt
    integer i
    parameter (x0=0.0, v0=1.0, L=10.0, t0=0.0, dt=0.01)
    t(1) = t0
    x(1) = x0
    v(1) = v0
    i = 2
    do while (i.le.6001)
        t(i) = (i-1)*dt
        if (x(i-1).lt.0.0 .or. x(i-1).gt.L) then
            v(i) = -v(i-1)
        else
            v(i) = v(i-1)
        end if
        x(i) = x(i-1)+v(i)*dt
        i = i+1
    end do
end subroutine
