MODULE pwl_integ
   use misc
   IMPLICIT NONE
   
CONTAINS   
   
   FUNCTION Pinteg(a, s, eps) RESULT(res)
      !             a
      !            /     -s
      !  ibpfunc = | dx x
      !            /
      !             1
      IMPLICIT NONE
      
      DOUBLEPRECISION :: res
      DOUBLEPRECISION, INTENT(IN) :: a, s, eps
      
      IF ((1d0 / 6d0) * LN3(a, eps) * (s - 1d0)**2 .gt. eps) THEN
         res = (1d0 - a**(1d0 - s)) / (s - 1d0)
      ELSE
         res = LN1(a, eps) - 5d-1 * LN2(a, eps) * (s - 1d0)
      ENDIF
      
      RETURN
   END FUNCTION Pinteg
   
   
   FUNCTION Qinteg(a, s, eps) RESULT(res)
      !             a
      !            /     -s
      !  ibqfunc = | dx x   log(x)
      !            /
      !            1
      IMPLICIT NONE
      
      DOUBLEPRECISION :: res
      DOUBLEPRECISION, INTENT(IN) :: a, s, eps
      
      IF (0.125d0 * LN4(a, eps) * (s - 1d0)**2.gt.eps) THEN
         res = (1d0 - a**(1d0 - s) * (1d0 + (s - 1d0) * LN1(a, eps))) / (s - 1d0)**2
      ELSE
         res = 5d-1 * LN2(a, eps) - (1d0 / 3d0) * LN3(a, eps) * (s - 1d0)
      ENDIF
      
      RETURN
   END FUNCTION Qinteg
   
   function Q2integ(a, s, eps) result(res)
      !              a
      !             /     -s    2
      !  ibq2func = | dx x   log (x)
      !             /
      !             1
      implicit none
      
      double precision, intent(in) :: a,s,eps
      double precision :: res
      
      if (0.1d0 * LN5(a, eps) * (s - 1d0)**2 > eps) then
         res = (2d0 * a**s + a * ( (s - 1d0) * LN1(a, eps) * ( LN1(a, eps) - s * LN1(a, eps) - 2d0 ) - 2d0 )) / (a**s * (s - 1d0)**3)
      else
         res = LN3(a, eps) / 3d0 - 0.25d0 * (s - 1d0) * LN4(a, eps)
      end if
      
      return
   end function Q2integ
   
   
   FUNCTION LN1(x, eps) RESULT(res)
      IMPLICIT NONE
      
      DOUBLEPRECISION :: res
      DOUBLEPRECISION, INTENT(IN) :: x, eps
      
      IF (5d-1 * (x - 1d0)**2.gt.eps) THEN
         res = dlog(x)
      ELSE
         res = x - 1d0
      ENDIF
      
      RETURN
   END FUNCTION LN1
   
   
   FUNCTION LN2(x, eps) RESULT(res)
      IMPLICIT NONE
      
      DOUBLEPRECISION :: res
      DOUBLEPRECISION, INTENT(IN) :: x, eps
      
      IF ((x - 1d0)**3.gt.eps) THEN
         res = dlog(x)**2
      ELSE
         res = (x - 1d0)**2
      ENDIF
      
      RETURN
   END FUNCTION LN2
   
   FUNCTION LN3(x, eps) RESULT(res)
      IMPLICIT NONE
      
      DOUBLEPRECISION :: res
      DOUBLEPRECISION, INTENT(IN) :: x, eps
      
      IF (1.5d0 * (x - 1d0)**4.gt.eps) THEN
         res = dlog(x)**3
      ELSE
         res = (x - 1d0)**3
      ENDIF
      
      RETURN
   END FUNCTION LN3
   
   FUNCTION LN4(x, eps) RESULT(res)
      IMPLICIT NONE
      
      DOUBLEPRECISION :: res
      DOUBLEPRECISION, INTENT(IN) :: x, eps
      
      IF (2d0 * (x - 1d0)**5.gt.eps) THEN
         res = dlog(x)**4
      ELSE
         res = (x - 1d0)**4
      ENDIF
      
      RETURN
   END FUNCTION LN4
   
   function LN5(x, eps) result(res)
      implicit none
      
      double precision :: res
      double precision, intent(in) :: x, eps
      
      if (2.5d0 * (x - 1d0)**6 > eps) then
         res = dlog(x)**5
      else
         res = (x - 1d0)**5
      end if
      
   end function LN5
   
END MODULE pwl_integ
