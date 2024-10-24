! A fortran95 program for G95
! By WQY
program main
  real a, b, c, x1, x2
  real realni, kompleksni, korijen
  write (*,*) "Square equation is given as: ax^2 + bx + c = 0, where a is not equal to 0."
  write (*,*) "Insert a:"
  read (*,*) a
  write (*,*) "Insert b:"
  read (*,*) b
  write (*,*) "Insert c:"
  read (*,*) c
  realni = -b/(2*a)
  korijen = (b**2)-4*a*c
  if (korijen.ge.0.0) then
    x1 = realni-SQRT(korijen)/(2*a)
    x2 = realni+SQRT(korijen)/(2*a)
    write (*, '("x1 = ",F6.2," , x2 = ",F6.2)') x1, x2
  else
    kompleksni = SQRT(-korijen)/(2*a)
    write (*, '("x1 = ",F6.2," +i",F6.2," , x2 = ",F6.2," -i",F6.2)') &
    + realni, kompleksni, realni, kompleksni
  end if
end
