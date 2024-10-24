! A fortran95 program for G95
! By WQY
program main
  real t(101), x(101), y(101), vx(101), vy(101)
  integer j
  call oscilator(t, x, y, vx, vy)
  open (10, file='osci.dat', status='new')
  write (10, '(A8," "A8," "A8," "A8," "A8)') "t", "x(t)", "y(t)", "vx(t)", "vy(t)"
  do j=1, 101
    write (10,'(F8.5," "F8.5," "F8.5," "F8.5," "F8.5)') t(j), x(j), y(j), vx(j), vy(j)
  end do
  close (10)
end

subroutine oscilator(t, x, y, vx, vy)
    implicit none
    real t(101), x(101), y(101), vx(101), vy(101)
    real dt, w1, w2, A
    integer i
    dt = 0.1
    w1 = 3.0
    w2 = 5.0
    A = 1.0
    i = 1
    do while (i.le.101)
        t(i) = (i-1)*dt
        x(i) = A*cos(w1*t(i))
        y(i) = A*sin(w2*t(i))
        vx(i) = -w1*A*sin(w1*t(i))
        vy(i) = w2*A*cos(w2*t(i))
        i = i+1
    end do
end subroutine
