module pwl_integ
   use data_types
   ! use misc
   implicit none
   
contains   
   
   function Pinteg(a, s, eps) result(res)
      !             a
      !            /     -s
      !  ibpfunc = | dx x
      !            /
      !             1
      implicit none
      real(dp) :: res
      real(dp), intent(in) :: a, s, eps
      if ((1d0 / 6d0) * LN3(a, eps) * (s - 1d0)**2 > eps) then
         res = (1d0 - a**(1d0 - s)) / (s - 1d0)
      else
         res = LN1(a, eps) - 5d-1 * LN2(a, eps) * (s - 1d0)
      endif
   end function Pinteg


   function Qinteg(a, s, eps) result(res)
      !             a
      !            /     -s
      !  ibqfunc = | dx x   log(x)
      !            /
      !            1
      implicit none
      real(dp) :: res
      real(dp), intent(in) :: a, s, eps
      if (0.125d0 * LN4(a, eps) * (s - 1d0)**2 > eps) then
         res = (1d0 - a**(1d0 - s) * (1d0 + (s - 1d0) * LN1(a, eps))) / (s - 1d0)**2
      else
         res = 5d-1 * LN2(a, eps) - (1d0 / 3d0) * LN3(a, eps) * (s - 1d0)
      endif
   end function Qinteg


   function Q2integ(a, s, eps) result(res)
      !              a
      !             /     -s    2
      !  ibq2func = | dx x   log (x)
      !             /
      !             1
      implicit none
      real(dp), intent(in) :: a,s,eps
      real(dp) :: res
      if (0.1d0 * LN5(a, eps) * (s - 1d0)**2 > eps) then
         res = (2d0 * a**s + a * ( (s - 1d0) * LN1(a, eps) * ( LN1(a, eps) - s * LN1(a, eps) - 2d0 ) - 2d0 )) / (a**s * (s - 1d0)**3)
      else
         res = LN3(a, eps) / 3d0 - 0.25d0 * (s - 1d0) * LN4(a, eps)
      end if
   end function Q2integ


   function LN1(x, eps) result(res)
      implicit none
      real(dp) :: res
      real(dp), intent(in) :: x, eps
      if (5d-1 * (x - 1d0)**2 > eps) then
         res = dlog(x)
      else
         res = x - 1d0
      endif
   end function LN1


   function LN2(x, eps) result(res)
      implicit none
      real(dp) :: res
      real(dp), intent(in) :: x, eps
      if ((x - 1d0)**3 > eps) then
         res = dlog(x)**2
      else
         res = (x - 1d0)**2
      endif
   end function LN2
   
   function LN3(x, eps) result(res)
      implicit none
      real(dp) :: res
      real(dp), intent(in) :: x, eps
      if (1.5d0 * (x - 1d0)**4 > eps) then
         res = dlog(x)**3
      else
         res = (x - 1d0)**3
      endif
   end function LN3


   function LN4(x, eps) result(res)
      implicit none
      real(dp) :: res
      real(dp), intent(in) :: x, eps
      if (2d0 * (x - 1d0)**5 > eps) then
         res = dlog(x)**4
      else
         res = (x - 1d0)**4
      endif
   end function LN4


   function LN5(x, eps) result(res)
      implicit none
      real(dp) :: res
      real(dp), intent(in) :: x, eps
      if (2.5d0 * (x - 1d0)**6 > eps) then
         res = dlog(x)**5
      else
         res = (x - 1d0)**5
      end if
   end function LN5

end module pwl_integ
