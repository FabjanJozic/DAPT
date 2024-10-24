! A fortran95 program for G95
! By WQY
program main
  integer n
  real x, rez
  x = 1.1
  n = 5
  call Horner(n, x, rez)
  write (*,*) rez
end

subroutine Horner(n, x, suma)
    implicit none
    integer n, i
    real x, suma, ai
    i = n
    write (*,*) "Input polynomial coefficient value:"
    read (*,*) ai
    suma = ai
    do while (i.ge.0)
        suma = suma*x
        write (*,*) "Input polynomial coefficient value:"
        read (*,*) ai
        suma = suma + ai
        i = i-1
    enddo
end subroutine
