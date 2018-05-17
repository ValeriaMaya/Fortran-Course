module roots

  USE, intrinsic :: iso_fortran_env, ONLY: WP => REAL64
  IMPLICIT NONE
CONTAINS

  SUBROUTINE biseccion(f,a,b,n,tol,raiz,clave)
    INTERFACE
      REAL(WP) FUNCTION f(x) ! Funci贸n que define la ecuaci贸n
        IMPORT :: WP
        IMPLICIT NONE
        REAL(WP), INTENT(IN) :: x
      END FUNCTION f
    END INTERFACE
    REAL(WP), INTENT(IN) :: a ! Extremo izquierdo del intervalo inicial
    REAL(WP), INTENT(IN) :: b ! Extremo derecho del intervalo inicial

    INTEGER, INTENT(INOUT) :: n ! L铆mite de iteraciones/iteraciones realizadas
    REAL(WP), INTENT(IN) :: tol ! Tolerancia para el error absoluto
    REAL(WP), INTENT(OUT) :: raiz ! Aproximaci贸n a la raiz
    INTEGER, INTENT(OUT) :: clave ! Clave de 茅xito:

    INTEGER :: i
    REAL(WP) :: xl, xr, error
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
  END SUBROUTINE biseccion




SUBROUTINE newton(f,df,a,n,tol,raiz,clave)
        INTERFACE
            REAL(WP) FUNCTION f(x) ! Funci贸n que define la ecuaci贸n
                IMPORT :: WP
                IMPLICIT NONE
                REAL(WP), INTENT(IN) :: x
            END FUNCTION f
            REAL(WP) FUNCTION df(x) ! Derivada de la funci贸n
                IMPORT :: WP ! que define a la ecuaci贸n
                IMPLICIT NONE
                REAL(WP), INTENT(IN) :: x
            END FUNCTION df
        END INTERFACE
        REAL(WP), INTENT(IN) :: a
        INTEGER, INTENT(INOUT) :: n ! L铆mite de iteraciones/iteraciones realizadas
        REAL(WP), INTENT(IN) :: tol ! Tolerancia para el error relativo
        REAL(WP), INTENT(OUT) :: raiz ! Aproximaci贸n a la ra铆z
        INTEGER, INTENT(OUT) :: clave ! Clave de 茅xito:
        INTEGER :: i
        REAL(WP) :: xx0
        clave = 1
        xx0 = a
        DO i=1,n
            raiz = xx0 - f(xx0)/df(xx0)
            IF (ABS(raiz-xx0) < tol*ABS(raiz) ) THEN
                clave = 0
                n = i
                EXIT
            ENDIF
            xx0 = raiz
        END DO
    END SUBROUTINE newton




  
  SUBROUTINE REGFAL(A, B, f, TOLER, RAIZ, ITER)
        IMPLICIT NONE
        Real(KIND(1D0)) :: f
        DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
        DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
        DOUBLE PRECISION, INTENT (out) :: RAIZ      !     Final approximation
        INTEGER, INTENT (out)          :: ITER   !     Number of iterations needed
        DOUBLE PRECISION :: X1, X2, F1, F2, FX
        X1 = A
        X2 = B
        F1 = f(X1)
        F2 = f(X2)
        IF (F1 * F2 >= 0.0D0) THEN
            PRINT *, "Limits are wrong"
            RETURN
        END IF
        ITER = 0
        DO
            RAIZ = X1 - F1 * (X2-X1) / (F2-F1)
            FX = f(RAIZ)
            !     Check IF root is found
            IF (ABS(FX) < TOLER) EXIT
            !     NEW limits
            IF (F1 * FX > 0.D0) THEN
                X1 = RAIZ
                F1 = FX
            ELSE
                X2 = RAIZ
                F2 = FX
            END IF
            ITER = ITER + 1
        END DO
    END SUBROUTINE REGFAL







  SUBROUTINE secante(f,x0,x1,n,tol,raiz,clave)
    INTERFACE
      REAL(WP) FUNCTION f(x) ! Funci贸n que define la ecuaci贸n
        IMPORT :: WP
        IMPLICIT NONE
        REAL(WP), INTENT(IN) :: x
      END FUNCTION f
    END INTERFACE
    REAL(WP), INTENT(IN) :: x0,x1 ! Aproximaciones iniciales a la ra铆z
    INTEGER, INTENT(INOUT):: n ! L铆mite de iteraciones/iteraciones realizadas
    REAL(WP), INTENT(IN) :: tol ! Tolerancia para el error relativo
    REAL(WP), INTENT(OUT) :: raiz ! Aproximaci贸n a la raiz
    INTEGER, INTENT(OUT) :: clave ! Clave de 茅xito:
    INTEGER :: i
    REAL(WP):: xx0, xx1, fx0, fx1
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
  END SUBROUTINE secante

  SUBROUTINE NEWTONMEJ(f,df,d2f,a,n,tol,raiz,clave)
    INTERFACE
      REAL(WP) FUNCTION f(x) ! Funci贸n que define la ecuaci贸n
        IMPORT :: WP
        IMPLICIT NONE
        REAL(WP), INTENT(IN) :: x
      END FUNCTION f
      REAL(WP) FUNCTION df(x) ! Derivada de la funci贸n
        IMPORT :: WP ! que define a la ecuaci贸n
        IMPLICIT NONE
        REAL(WP), INTENT(IN) :: x
      END FUNCTION df
      REAL(WP) FUNCTION d2f(x) ! Derivada de la funci贸n
        IMPORT :: WP ! que define a la ecuaci贸n
        IMPLICIT NONE
        REAL(WP), INTENT(IN) :: x
      END FUNCTION d2f
    END INTERFACE
    REAL(WP), INTENT(IN) :: a
    INTEGER, INTENT(INOUT) :: n ! L铆mite de iteraciones/iteraciones realizadas
    REAL(WP), INTENT(IN) :: tol ! Tolerancia para el error relativo
    REAL(WP), INTENT(OUT) :: raiz ! Aproximaci贸n a la ra铆z
    INTEGER, INTENT(OUT) :: clave ! Clave de 茅xito:
    INTEGER :: i
    REAL(WP) :: xx0
    clave = 1
    xx0 = a
    DO i=1,n
      raiz = xx0 - (( f(xx0) * df(xx0) )/( (df(xx0))**2 - f(xx0) * d2f(xx0) ))
      IF (ABS(raiz-xx0) < tol*ABS(raiz) ) THEN
        clave = 0
        n = i
        EXIT
      ENDIF
      xx0 = raiz
    END DO
  END SUBROUTINE NEWTONMEJ

SUBROUTINE modifinewton(f,df,ddf,x0,n,tol,raiz,clave)
! ---------------------------------------------------
! Metodo DE NEWTON-RAPHSON para encontrar una
! soluci髇 de f(x)=0 dada la funci髇 derivable
! f y una aproximaci髇 inicial x0.
! ---------------------------------------------------
! Bloque de declaraci髇 de argumentos
! ---------------------------------------------------
INTERFACE
REAL(WP) FUNCTION f(x) ! Funci髇 que define la ecuaci髇
IMPORT :: WP
IMPLICIT NONE
REAL(WP), INTENT(IN) :: x
END FUNCTION f
REAL(WP) FUNCTION df(x) ! Derivada de la funci髇
IMPORT :: WP ! que define a la ecuaci髇
IMPLICIT NONE
REAL(WP), INTENT(IN) :: x
END FUNCTION df
REAL(WP) FUNCTION ddf(x)
IMPORT:: WP
IMPLICIT NONE
REAL(WP), INTENT(IN)::x
END FUNCTION ddf
END INTERFACE
REAL(WP), INTENT(IN) :: x0 ! Aproximaci髇 inicial a la ra韟
INTEGER, INTENT(INOUT) :: n ! L韒ite de iteraciones/iteraciones realizadas
REAL(WP), INTENT(IN) :: tol ! Tolerancia para el error relativo
REAL(WP), INTENT(OUT) :: raiz ! Aproximaci髇 a la ra韟
INTEGER, INTENT(OUT) :: clave ! Clave de 閤ito:
! 0 : 閤ito
! >0 : iteraciones excedidas
! ---------------------------------------------------
! Declaraci髇 de variables locales
! ---------------------------------------------------
INTEGER :: i
REAL(WP) :: xx0

! Bloque de procesamiento
! ---------------------------------------------------

xx0 = x0
DO i=1,n
raiz = xx0 - (f(xx0)*df(xx0))/((df(xx0))**2-(f(xx0)*ddf(xx0)))
IF (ABS(raiz-xx0) < tol*ABS(raiz) ) THEN
clave = 0
n = i
Exit
ENDIF
xx0 = raiz
if (n==i)clave=1

END DO
END SUBROUTINE modifinewton

end module roots
