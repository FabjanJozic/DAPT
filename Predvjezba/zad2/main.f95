! A fortran95 program for G95
! By WQY
program main
  real v0, h0, g, t0, v, h, t
  g = 9.81
  write (*,*) "The object is free-falling, or moving along z-axis. Up is +z and down is -z."
  write (*,*) "Input initial velocity: "
  read (*,*) v0
  write (*,*) "Input initial hight: "
  read (*,*) h0
  write (*,*) "Input time: "
  read (*,*) t0
  t = v0/g+SQRT((v0**2)+(2*h0)/g)
  if (t0.lt.t) then
    v = -g*t0+v0
    h = -0.5*g*(t0**2)+v0*t0+h0
  else
    v = -g*t+v0
    h = 0.0
  end if
  write (*, '("v = ",F7.2," , h = ",F6.2)') v, h
end
