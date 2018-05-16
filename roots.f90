MODULE roots
USE, intrinsic :: iso_fortran_env, ONLY: WP => REAL64
IMPLICIT NONE
CONTAINS

SUBROUTINE biseccion(f,a,b,n,tol,raiz,clave)
! ---------------------------------------------------
! METODO DE BISECCION para encontrar una soluciÃ³n
! de f(x)=0 dada la funciÃ³n continua f en el intervalo
! [a,b] donde f(a) y f(b) tienen signos opuestos.
! ---------------------------------------------------
! Bloque de declaraciÃ³n de argumentos
! ---------------------------------------------------
INTERFACE
REAL(WP) FUNCTION f(x) ! FunciÃ³n que define la ecuaciÃ³n
IMPORT :: WP
IMPLICIT NONE
REAL(WP), INTENT(IN) :: x
END FUNCTION f
END INTERFACE
REAL(WP), INTENT(IN) :: a ! Extremo izquierdo del intervalo inicial
REAL(WP), INTENT(IN) :: b ! Extremo derecho del intervalo inicial

INTEGER, INTENT(INOUT) :: n ! LÃ­mite de iteraciones/iteraciones realizadas
REAL(WP), INTENT(IN) :: tol ! Tolerancia para el error absoluto
REAL(WP), INTENT(OUT) :: raiz ! AproximaciÃ³n a la raiz
INTEGER, INTENT(OUT) :: clave ! Clave de Ã©xito:
! 0 : Ã©xito
! >0 : iteraciones excedidas
! <0 : no se puede proceder (f de igual signo
! en a y b)
! ---------------------------------------------------
! Bloque de declaraciÃ³n de variables locales
! ---------------------------------------------------
INTEGER :: i
REAL(WP) :: xl, xr, error
! ---------------------------------------------------
! Bloque de procesamiento
! ---------------------------------------------------
clave = 1
xl = a
xr = b
IF (SIGN(1.0_WP,f(xl))*SIGN(1.0_WP,f(xr)) > 0.0_WP) THEN
clave = -1
RETURN
ENDIF
DO i=1,n
error = (xr-xl)*0.5_WP
raiz = xl + error
IF (error < tol) THEN
clave = 0
n = i
EXIT
ENDIF
IF ( SIGN(1.0_WP,f(xl))* SIGN(1.0_WP,f(raiz)) > 0.0_WP) THEN
xl = raiz
ELSE
xr = raiz
ENDIF
ENDDO
! ---------------------------------------------------
END SUBROUTINE biseccion








SUBROUTINE newton(f,df,x0,n,tol,raiz,clave)
! ---------------------------------------------------
! Metodo DE NEWTON-RAPHSON para encontrar una
! soluciÃ³n de f(x)=0 dada la funciÃ³n derivable
! f y una aproximaciÃ³n inicial x0.
! ---------------------------------------------------
! Bloque de declaraciÃ³n de argumentos
! ---------------------------------------------------
INTERFACE
REAL(WP) FUNCTION f(x) ! FunciÃ³n que define la ecuaciÃ³n
IMPORT :: WP
IMPLICIT NONE
REAL(WP), INTENT(IN) :: x
END FUNCTION f
REAL(WP) FUNCTION df(x) ! Derivada de la funciÃ³n
IMPORT :: WP ! que define a la ecuaciÃ³n
IMPLICIT NONE
REAL(WP), INTENT(IN) :: x
END FUNCTION df
END INTERFACE
REAL(WP), INTENT(IN) :: x0 ! AproximaciÃ³n inicial a la raÃ­z
INTEGER, INTENT(INOUT) :: n ! LÃ­mite de iteraciones/iteraciones realizadas
REAL(WP), INTENT(IN) :: tol ! Tolerancia para el error relativo
REAL(WP), INTENT(OUT) :: raiz ! AproximaciÃ³n a la raÃ­z
INTEGER, INTENT(OUT) :: clave ! Clave de Ã©xito:
! 0 : Ã©xito
! >0 : iteraciones excedidas
! ---------------------------------------------------
! DeclaraciÃ³n de variables locales
! ---------------------------------------------------
INTEGER :: i
REAL(WP) :: xx0
! ---------------------------------------------------
! Bloque de procesamiento
! ---------------------------------------------------
clave = 1
xx0 = x0
DO i=1,n
raiz = xx0 - f(xx0)/df(xx0)
IF (ABS(raiz-xx0) < tol*ABS(raiz) ) THEN
clave = 0
n = i
EXIT
ENDIF
xx0 = raiz
END DO
! ---------------------------------------------------
END SUBROUTINE newton










!     Regula falsi (false position) method
SUBROUTINE REGFAL(A, B, FUNC, TOLER, X, ITER)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
  DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
  DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
  INTEGER, INTENT (out)          :: ITER   !     Number of iterations needed
  INTEGER, EXTERNAL :: FUNC


  DOUBLE PRECISION :: X1, X2, F1, F2, FX

  X1 = A
  X2 = B

  F1 = FUNC(X1)
  F2 = FUNC(X2)

  IF (F1 * F2 >= 0.0D0) THEN
     PRINT *, "Limits are wrong"
     RETURN
  ENDIF

  ITER = 0

  DO
     X = X1 - F1 * (X2-X1) / (F2-F1)
     FX = FUNC(X)

     !     Check IF root is found
     IF (ABS(FX) < TOLER) EXIT

     !     NEW limits
     IF (F1 * FX > 0.D0) THEN
        X1 = X
        F1 = FX
     ELSE
        X2 = X
        F2 = FX
     END IF

     ITER = ITER + 1
  ENDDO

END SUBROUTINE REGFAL









SUBROUTINE secante(f,x0,x1,n,tol,raiz,clave)
! ---------------------------------------------------
! ALGORITMO DE LA SECANTE para encontrar una soluciÃ³n
! de f(x)=0, siendo f una funciÃ³n continua, dada las
! aproximaciones iniciales x0 y x1.
! ---------------------------------------------------
! Bloque de declaraciÃ³n de argumentos
! ---------------------------------------------------
INTERFACE
REAL(WP) FUNCTION f(x) ! FunciÃ³n que define la ecuaciÃ³n
IMPORT :: WP
IMPLICIT NONE
REAL(WP), INTENT(IN) :: x
END FUNCTION f
END INTERFACE
REAL(WP), INTENT(IN) :: x0,x1 ! Aproximaciones iniciales a la raÃ­z
INTEGER, INTENT(INOUT):: n ! LÃ­mite de iteraciones/iteraciones realizadas
REAL(WP), INTENT(IN) :: tol ! Tolerancia para el error relativo
REAL(WP), INTENT(OUT) :: raiz ! AproximaciÃ³n a la raiz
INTEGER, INTENT(OUT) :: clave ! Clave de Ã©xito:
! 0 : Ã©xito
! >0 : iteraciones excedidas
! ---------------------------------------------------
! Bloque de declaraciÃ³n de variables locales
! ---------------------------------------------------
INTEGER :: i
REAL(WP):: xx0, xx1, fx0, fx1
! ---------------------------------------------------
! Bloque de procesamiento
! ---------------------------------------------------
clave = 1
xx0 = x0
xx1 = x1
fx0 = f(x0)
fx1 = f(x1)
DO i= 2,n
raiz = xx1 - fx1*((xx1-xx0)/(fx1-fx0))
IF (ABS(raiz-xx1) < tol*ABS(raiz)) THEN
clave = 0
n = i
EXIT
ENDIF
xx0 = xx1
fx0 = fx1
xx1 = raiz
fx1 = f(raiz)
END DO
! ---------------------------------------------------
END SUBROUTINE secante



END MODULE roots
