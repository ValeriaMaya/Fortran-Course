program principal

  use, intrinsic :: iso_fortran_env, only: WP => REAL64
  use roots
  implicit none
  real(WP) ::   a, b, tol, raiz
  integer :: n, clave

  a = 1.0_WP
  b = 4.0_WP
  tol = 0.5E-6_WP
  n = 100

  Call biseccion(f,a,b,n,tol,raiz,clave)
  write(*,*) ' ---Biseccion--- '
  if (clave == 0) then
    write(*,*) ' Raíz =  ', raiz
    write(*,*) ' Numero de iteraciones realizadas = ',n
  else
    write(*,*) ' Error = ',clave
  end if

  a = 1.0_WP
  b = 4.0_WP
  tol = 0.5E-6_WP
  n = 100

  Call newton(f,df,a,n,tol,raiz,clave)
  write(*,*) ' ---Newton-Raphson--- '
  if (clave == 0) then
    write(*,*) ' raiz =  ', raiz
    write(*,*) ' Numero de iteraciones realizadas = ',n
  else
    write(*,*) ' Error = ',clave
  end if

  a = 1.0_WP
  b = 4.0_WP
  tol = 0.5E-6_WP
  n = 100


  Call REGFAL(a, b, f, tol, raiz, n)
  write(*,*) ' ---Regla falsa--- '
  if (clave == 0) then
    write(*,*) ' raiz =  ', raiz
    write(*,*) ' Numero de iteraciones realizadas = ',n
  else
    write(*,*) ' Error = ',clave
  end if

  a = 1.0_WP
  b = 4.0_WP
  tol = 0.5E-6_WP
  n = 100


  Call secante(f,a,b,n,tol,raiz,clave)
  write(*,*) ' ---Metodo de la secante--- '
  if (clave == 0) then
    write(*,*) ' raiz =  ', raiz
    write(*,*) ' Numero de iteraciones realizadas = ',n
  else
    write(*,*) ' Error = ',clave
  end if


  a = 1.0_WP
  b = 4.0_WP
  tol = 0.5E-6_WP
  n = 100


  CALL modifinewton(f,df,d2f,a,n,tol,raiz,clave)
  write(*,*) ' ---Newton-Raphson (Mejorado)--- '
    if (clave == 0) then
      write(*,*) ' raiz =  ', raiz
      write(*,*) ' Numero de iteraciones realizadas = ',n
    else
      write(*,*) ' Error = ',clave
    end if

contains

  function f(x)
    implicit none
    Real(KIND(1D0)), INTENT (in) :: x
    Real(KIND(1D0)) :: f
    f = (x**2)-x-6
  end function

  function df(x)
    implicit none
    Real(KIND(1D0)), INTENT (in) :: x
    Real(KIND(1D0)) :: df
    df = (2*x)-1
  end function

  function d2f(x)
    implicit none
    Real(KIND(1D0)), INTENT (in) :: x
    Real(KIND(1D0)) :: d2f
    d2f = 2
  end function

end program principal
