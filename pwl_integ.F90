module pwl_integ
   use data_types
   implicit none

contains

   function Pinteg(a, s, eps) result(res)
      !             a
      !            /     -s
      !  P(a, s) = | dx x
      !            /
      !             1
      implicit none
      real(dp) :: res
      real(dp), intent(in) :: a, s, eps
      if ( LN3(a, eps) * (s - 1d0)**2 > 6d0 * eps) then
         res = (1d0 - a**(1d0 - s)) / (s - 1d0)
      else
         res = LN1(a, eps) - 0.5d0 * LN2(a, eps) * (s - 1d0)
      endif
   end function Pinteg


   function Qinteg(a, s, eps) result(res)
      !             a
      !            /     -s
      !  Q(a, s) = | dx x   log(x)
      !            /
      !            1
      implicit none
      real(dp) :: res
      real(dp), intent(in) :: a, s, eps
      if ( 0.125d0 * LN4(a, eps) * (s - 1d0)**2 > eps ) then
         res = ( 1d0 - a**(1d0 - s) * (1d0 + (s - 1d0) * LN1(a, eps)) ) / (s - 1d0)**2
      else
         res = 0.5d0 * LN2(a, eps) - LN3(a, eps) * (s - 1d0) / 3d0
      endif
   end function Qinteg


   function Q2integ(a, s, eps) result(res)
      !              a
      !             /     -s    2
      !  Q2(a, s) = | dx x   log (x)
      !             /
      !             1
      implicit none
      real(dp), intent(in) :: a, s, eps
      real(dp) :: res
      if ( 0.1d0 * LN5(a, eps) * (s - 1d0)**2 > eps ) then
         res = (2d0 - a**(1d0 - s)) * ( 2d0 + (s - 1d0) * (2d0 + (s - 1d0) * LN1(a, eps)) * LN1(a, s) ) / (s - 1d0)**3
      else
         res = LN3(a, eps) / 3d0 - 0.25d0 * (s - 1d0) * LN4(a, eps)
      end if
   end function Q2integ

   !     Eq. (3.19) of Rueda-Becerril (2017)
   function sscR(a, b, c, d, alpha, beta) result(res)
      implicit none
      real(dp), intent(in) :: a, b, c, d, alpha, beta
      real(dp) :: res
      real(dp) :: eps = 1d-9
      res = c**(beta + alpha + 2d0) * a**(alpha + 1d0) * Pinteg(d / c, -(beta + alpha + 1d0), eps) * Qinteg(b / a, -alpha, eps)
   end function sscR
   
   !     Eq. (3.21) of Rueda-Becerril (2017)
   function sscS(a, b, c, d, alpha, beta) result(res)
      implicit none
      real(dp), intent(in) :: a, b, c, d, alpha, beta
      real(dp) :: res, b_ac, d_c
      real(dp) :: eps = 1d-9
      b_ac = b / (a * c)
      d_c = d / c
      if ( LN3(b / (a * d), eps) * (alpha + 1d0)**2 > eps * 6d0 ) then
         res = c**(beta + 1d0) * ( b**(alpha + 1d0) * Pinteg(d_c, -beta, eps) - &
            (a * c)**(alpha + 1d0) * Pinteg(d_c, -(alpha + beta + 1d0), eps) ) / (alpha + 1d0)
      else
         res = c**(beta + 1d0) * b**(alpha + 1d0) * ( LN1(b_ac, eps) * &
            Pinteg(d_c, -beta, eps) - Qinteg(d_c, -beta, eps) + 0.5d0 * &
            (alpha + 1d0) * ( LN2(b_ac, eps) * Pinteg(d_c, -beta, eps) - &
            2d0 * LN1(b_ac, eps) * Qinteg(d_c, -beta, eps) + &
            Q2integ(d_c, -beta, eps) ) )
      endif
   end function sscS
   
   !     Eq. 2.106 of Mimica (2004)
   function sscG1ISO(a, b, c, d, alpha, beta) result(res)
      implicit none
      real(dp), intent(in) :: a, b, c, d, alpha, beta
      real(dp) :: res
      res = sscR(a, b, c, d, alpha, beta) - sscR(a, b, c, d, alpha + 1d0, beta)
   end function sscG1ISO
   
   !     Eq. 2.107 of Mimica (2004)
   function sscG2ISO(a, b, c, d, alpha, beta) result(res)
      implicit none
      real(dp), intent(in) :: a, b, c, d, alpha, beta
      real(dp) :: res
      res = sscS(a, b, c, d, alpha, beta) - sscS(a, b, c, d, alpha + 1d0, beta)
   end function sscG2ISO


   !
   !   -----{   Power-law distribution Normalization   }-----
   !
   function pwl_norm(norm, index, e1, e2) result(k)
      implicit none
      real(dp), intent(in) :: norm, index, e1, e2
      real(dp), parameter :: eps = 1e-9
      real(dp) :: k, integ
      integ = e1**(1.0 - index) * Pinteg(e2 / e1, index, eps)
      k = 1d0 / (norm * integ)
   end function pwl_norm



   ! ===========================================================================
   !  #       ####   ####     #####   ####  #    # ###### #####   ####
   !  #      #    # #    #    #    # #    # #    # #      #    # #
   !  #      #    # #         #    # #    # #    # #####  #    #  ####
   !  #      #    # #  ###    #####  #    # # ## # #      #####       #
   !  #      #    # #    #    #      #    # ##  ## #      #   #  #    #
   !  ######  ####   ####     #       ####  #    # ###### #    #  ####
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
   ! ===========================================================================

end module pwl_integ
