! *************************************************************************
!  Miscellaneous functions and routines
!  ====================================
!
!  - polint      : subroutine
!  - an_error    : subroutine
!  - arth        : function
!  - locate      : function
!  - chebev      : function
!  - char2int    : function
!  - char2double : function
!  - dbesseljn   : function
!!!TODO: Update list
!
! *************************************************************************
module misc
   use data_types
   implicit none

   !> De-allocate and reallocate array
   interface realloc
      module procedure realloc_1d
      module procedure realloc_2d
   end interface

contains

   ! #    # ##### # #      # ##### # ######  ####
   ! #    #   #   # #      #   #   # #      #
   ! #    #   #   # #      #   #   # #####   ####
   ! #    #   #   # #      #   #   # #           #
   ! #    #   #   # #      #   #   # #      #    #
   !  ####    #   # ###### #   #   # ######  ####

   !> De-allocate and reallocate 1D array
   subroutine realloc_1d(array, size)
      implicit none
      integer, intent(in) :: size
      real(dp), intent(inout), allocatable, dimension(:) :: array
      if ( allocated(array) ) deallocate(array)
      allocate(array(size))
   end subroutine realloc_1d

   !> De-allocate and reallocate 2D array
   subroutine realloc_2d(array, size1, size2)
      implicit none
      integer, intent(in) :: size1, size2
      real(dp), intent(inout), allocatable, dimension(:, :) :: array
      if ( allocated(array) ) deallocate(array)
      allocate(array(size1, size2))
   end subroutine realloc_2d

   !> This produces an error message
   subroutine an_error(string)
      implicit none
      character(len=*), intent(IN) :: string
      write (*, "('ERROR ', A)") string
      stop 'program terminated by an_error'
   end subroutine an_error

   !> This produces a warning message
   subroutine a_warning(string)
      implicit none
      character(len=*), intent(in) :: string
      write(*,*) 'warging message: '//trim(string)
   end subroutine a_warning

   !> This creates an arithmetic progression array
   function arth(first, increment, n)
      implicit none
      integer, parameter :: NPAR_ARTH=16, NPAR2_ARTH=8
      real(dp), intent(in) :: first, increment
      integer, intent(in) :: n
      real(dp), dimension(n) :: arth
      integer :: k, k2
      real(dp) :: temp
      if (n > 0) arth(1) = first
      if (n <= NPAR_ARTH) then
         do k=2,n
            arth(k) = arth(k - 1) + increment
         end do
      else
         do k=2,NPAR2_ARTH
            arth(k) = arth(k - 1) + increment
         end do
         temp = increment * NPAR2_ARTH
         k = NPAR2_ARTH
         do
            if (k >= n) exit
            k2 = k + k
            arth(k + 1:min(k2, n)) = temp + arth(1:min(k, n - k))
            temp = temp + temp
            k = k2
         end do
      end if
   end function arth

   !> Find closest element in an array to a value
   function locate(xx, x, in_bounds)
      implicit none
      real(dp), intent(in) :: x
      real(dp), intent(in), dimension(:) :: xx
      logical, intent(in), optional :: in_bounds
      integer :: locate
      integer :: n, jl, jm, ju
      logical :: ascnd, bounds
      !!!NOTE: The flag 'bounds' tells the function to return the index of the
      !!!      lower value of the interval of 'xx' in which 'x' is
      if ( present(in_bounds) ) then
         bounds = in_bounds
      else
         bounds = .false.
      end if
      n = size(xx)
      ascnd = ( xx(n) >= xx(1) )
      jl = 0
      ju = n + 1
      do
         if ( ju - jl <= 1 ) exit
         jm = (ju + jl) / 2
         if ( ascnd .eqv. (x >= xx(jm)) ) then
            jl = jm
         else
            ju = jm
         end if
      end do
      if ( dabs(x - xx(1)) < 1d-9 ) then
         locate = 1
         return
      else if ( dabs(x - xx(n)) < 1d-9 ) then
         locate = n - 1
         return
      else
         if ( jl < 1 .and. bounds ) then
            locate = 1
            return
         else if ( jl >= n .and. bounds ) then
            locate = n - 1
            return
         else
            !call a_warning('value out of array bounds in locate function.')
            locate = jl
            return
         end if
      end if
   end function locate

   !> Checking equality of integers
   function assert_eq(nn, string)
      implicit none
      character(len=*), intent(in) :: string
      integer, intent(in), dimension(:) :: nn
      integer :: assert_eq
      if (all(nn(2:) == nn(1))) then
         assert_eq = nn(1)
      else
         write(*,*) 'ERROR: an assert_eq failed with this tag:', string
         stop 'program terminated by assert_eq'
      end if
   end function assert_eq

   !> Convert a integer to a string value using an internal read.
   function int2char(i) result(c)
      integer, intent(in) :: i
      character(len=60) :: c
      write (c,*) i
      c = adjustl(c)
   end function int2char

   !> Convert a string to a integer value using an internal read.
   function char2int(c) result(i)
      character(len=*), intent(in) :: c
      integer :: i
      read (c,'(I5)') i
   end function char2int

   !> Convert a string to a double value using an internal read.
   function char2double(c) result(d)
      character(len=*), intent(in) :: c
      real(dp) :: d
      read(c, *) d
   end function char2double

   !> A zeros (or almost zero) valued array
   !! @parameter n integer size of the array
   !! @parameter small logical 1d-200 (true) or 0d0 (false)
   function zeros1D(n, small) result(a)
      implicit none
      integer, intent(in) :: n
      logical, intent(in) :: small
      real(dp), dimension(n) :: a
      if (small) then
         a = 1d-200
      else
         a = 0d0
      end if
   end function zeros1D

   !> Count lines of a text file
   function count_lines(filename) result(nlines)
      implicit none
      character(len=*), intent(in) :: filename
      integer :: nlines, io
      open(10, file=filename, iostat=io, status="old")
      if ( io /= 0 ) stop "Cannot open file!"
      nlines = 0
      do
         read(10, *, iostat=io)
         if (io/=0) exit
         nlines = nlines + 1
      end do
      close(10)
   end function count_lines


   ! #    # #    # #    # ###### #####  #  ####   ####
   ! ##   # #    # ##  ## #      #    # # #    # #
   ! # #  # #    # # ## # #####  #    # # #       ####
   ! #  # # #    # #    # #      #####  # #           #
   ! #   ## #    # #    # #      #   #  # #    # #    #
   ! #    #  ####  #    # ###### #    # #  ####   ####
   !
   !> Polinomial interpolation
   !! @parameter xa input array
   !! @parameter ya input array
   !! @parameter x input value where y is interpolated
   !! @parameter y output interpolated value at y
   !! @parameter dy output error
   subroutine polint(xa, ya, x, y, dy)
      implicit none
      real(dp), intent(in) :: x
      real(dp), dimension(:), intent(in) :: xa, ya
      real(dp), intent(out) :: y, dy
      integer :: m, n, ns
      integer, dimension(1) :: iminloc
      real(dp), dimension(size(xa)) :: c, d, den, ho
      if ( size(xa) == size(ya) ) then
         n = size(xa)
      else
         stop 'polint: xa and ya have different size'
      end if
      c = ya
      d = ya
      ho = xa - x
      iminloc = minloc(abs(x - xa))
      ns = iminloc(1)
      y = ya(ns)
      ns = ns - 1
      do m = 1, n - 1
         den(1:n - m) = ho(1:n - m) - ho(1 + m:n)
         if ( any( dabs(den(1:n - m)) == 0d0)) then
            write(*, *) x, xa
            call an_error('polint: calculation failure')
         end if
         den(1:n - m) = (c(2:n - m + 1) - d(1:n - m)) / den(1:n - m)
         d(1:n - m) = ho(1 + m:n) * den(1:n - m)
         c(1:n - m) = ho(1:n - m) * den(1:n - m)
         if (2 * ns .lt. n - m) then
            dy = c(ns + 1)
         else
            dy = d(ns)
            ns = ns - 1
         end if
         y = y + dy
      end do
   end subroutine polint


   !> Linear interpolation
   subroutine linint(x1, x2, x, y1, y2, y)
   implicit none
   real(dp), intent(out) :: y
   real(dp), intent(in) :: y1, y2, x1, x2, x
   real(dp) :: dx1, dx2, dx, dy1, dy2, dy3
   dx = abs(x2 - x1)
   dx1 = abs(x1 - x)
   dx2 = abs(x2 - x)
   if (dx1 + dx2 > dx) then
      write(*, *) x1, x2, x
      call an_error('linint: x not between x1 and x2')
   end if
   y = y1 + ( y2 - y1 ) * ( x - x1 ) / dx
   dy1 = abs(y2 - y1)
   dy2 = abs(y1 - y)
   dy3 = abs(y2 - y)
   if(dy3 + dy2 > dy1) then
      write(*, *) y1, y2, y
      call an_error('linint: Interpolation not between y1 and y2')
   end if
   end subroutine linint

   !> Chebychev evaluation
   function chebev(x, coef, num_coefs, xmin, xmax) result(res)
      implicit none
      integer, intent(in) :: num_coefs
      real(dp), intent(in) :: x,xmin,xmax
      real(dp), intent(in), dimension(:) :: coef
      integer :: j
      real(dp) :: d,dd,y,y2,sv,res
      if ( (x - xmin) * (x - xmax) > 0d0 ) then
         print*, x, xmin, xmax
         write(*, *) 'x is not in rage in chebev'
         stop
      end if
      d = 0d0
      dd = 0d0
      y = (2d0 * x - xmin - xmax) / (xmax - xmin)
      y2 = 2d0 * y
      do j=num_coefs, 2, -1
         sv = dd
         dd = d
         d = y2 * dd - sv + coef(j)
      end do
      res = y * d - dd + 0.5d0 * coef(1)
   end function chebev

   !> Solver of a tridiagonal matrix. Based on the code in "Numberical Recipes".
   subroutine tridag_ser(a, b, c, r, u)
      implicit none
      real(dp), dimension(:), intent(in) :: a, b, c, r
      real(dp), dimension(:), intent(out) :: u
      real(dp), dimension(size(b)) :: gam
      integer :: n, j
      real(dp) :: bet
      n = assert_eq((/ size(a) + 1, size(b), size(c) + 1, size(r), size(u) /), 'tridag_ser')
      bet = b(1)
      if ( bet == 0.0d0 ) call an_error('tridag_ser: error at code stage 1')
      u(1) = r(1) / bet
      do j=2, n
         gam(j) = c(j - 1) / bet
         bet = b(j) - a(j - 1) * gam(j)
         if ( bet == 0.0d0 ) call an_error('tridag_ser: error at code stage 2')
         u(j) = (r(j) - a(j - 1) * u(j - 1)) /bet
      end do
      do j = n-1, 1, -1
         u(j) = u(j) - gam(j + 1) * u(j + 1)
      end do
   end subroutine tridag_ser

   !> Trapezoidal rule
   subroutine trapzd_w2arg(func, a, b, s, n, p)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: a, b, p
      real(dp), intent(inout) :: s
      interface
         function func(x, pp) result(res)
            use data_types
            real(dp), intent(in) :: pp
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(size(x, dim=1)) :: res
         end function func
      end interface
      integer :: it
      real(dp) :: del, fsum
      if ( n == 1 ) then
         s = 0.5d0 * (b - a)*  sum(func( (/ a, b /), p) )
      else
         it = 2**(n - 2)
         del =(b - a) / it
         fsum = sum( func(arth(a + 0.5d0 * del, del, it), p) )
         s = 0.5d0 * (s + del * fsum)
      end if
   end subroutine trapzd_w2arg

   !> Romberg integrator
   function qromb_w2arg(func, a, b, p) result(qromb)
      implicit none
      real(dp), intent(in) :: a, b, p
      real(dp) :: qromb
      interface
         function func(x, pp) result(res)
            use data_types
            real(dp), intent(in) :: pp
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(size(x, dim=1)) :: res
         end function func
      end interface
      integer, parameter :: jmax=20, jmaxp=jmax+1, k=5, km=k-1
      real(dp), parameter :: eps=1d-5
      integer :: j
      real(dp) :: dqromb
      real(dp), dimension(jmaxp) :: h, s
      h(1) = 1d0
      do j=1, jmax
         call trapzd_w2arg(func, a, b, s(j), j, p)
         if ( j >= k ) then
            call polint(h(j - km:j), s(j - km:j), 0.0d0, qromb, dqromb)
            if ( dabs(dqromb) <= eps * dabs(qromb) ) return
         end if
         s(j + 1)=s(j)
         h(j + 1) = 0.25d0 * h(j)
      end do
      call an_error('qromb: too many steps')
   end function qromb_w2arg

   !> Trapezoidal rule with arrays
   subroutine trapzd_arr(x, a, b, f, s, n)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: a, b
      real(dp), intent(in), dimension(:) :: f, x
      real(dp), intent(inout) :: s
      integer :: it, i
      real(dp) :: del, fsum, fa, fb, df, xx, fx
      if ( n == 1 ) then
         call polint(x, f, a, fa, df)
         call polint(x, f, b, fb, df)
         s = 0.5d0 * (b - a) * (fa + fb)
      else
         it = 2**(n - 2)
         del = (b - a) / dble(it)
         xx = a + 0.5d0 * del
         fsum = 0d0
         do i=1, it
            call polint(x, f, xx, fx, df)
            fsum = fsum + fa
            xx = xx + del
         end do
         s = 0.5d0 * (s + del * fsum)
      end if
   end subroutine trapzd_arr

   !> Romberg integrator with arrays
   function qromb_arr(x, a, b, f) result(res)
      implicit none
      real(dp), intent(in) :: a, b
      real(dp), intent(in), dimension(:) :: x, f
      integer, parameter :: jmax=30, jmaxp=jmax+1, k=10, km=k-1
      real(dp), parameter :: eps=1d-5
      integer :: j
      real(dp) :: dqromb, res
      real(dp), dimension(jmaxp) :: h, s
      h(1)=1d0
      do j=1, jmax
         call trapzd_arr(x, a, b, f, s(j), j)
         if ( j >= k ) then
            call polint(h(j-km:j), s(j-km:j), 0.0d0, res, dqromb)
            if ( dabs(dqromb) <= eps * dabs(res) ) return
         end if
         s(j + 1) = s(j)
         h(j + 1) = 0.25d0 * h(j)
      end do
      call an_error('qromb_flux: too many steps')
   end function qromb_arr

   !> 2nd order Runge-Kutta for arrays
   subroutine rk2_arr(y, dydx, h, yout)
      implicit none
      real(dp), intent(in) :: h, y
      real(dp), dimension(:), intent(in) :: dydx
      real(dp), intent(out) :: yout
      integer :: ndum
      ndum = size(dydx)
      yout = y + 0.5d0 * h * (dydx(ndum-1) + dydx(ndum))
   end subroutine rk2_arr

   !> 4th order Runge-Kutta
   subroutine rk4(y, dydx, x, h, yout, derivs)
      implicit none
      real(dp), dimension(:), intent(in) :: y, dydx
      real(dp), intent(in) :: x, h
      real(dp), dimension(:), intent(out) :: yout
      interface
         subroutine derivs(x, y, dydx)
            use data_types
            implicit none
            real(dp), intent(in) :: x
            real(dp), dimension(:), intent(in) :: y
            real(dp), dimension(:), intent(out) :: dydx
         end subroutine derivs
      end interface
      integer :: ndum
      real(dp) :: h6, hh, xh
      real(dp), dimension(size(y)) :: dym, dyt, yt
      ndum = assert_eq((/ size(y), size(dydx), size(yout) /), 'rk4')
      hh = h * 0.5d0
      h6 = h / 6.0d0
      xh = x + hh
      yt = y + hh * dydx
      call derivs(xh, yt, dyt)
      yt = y + hh * dyt
      call derivs(xh, yt, dym)
      yt = y + h * dym
      dym = dyt + dym
      call derivs(x + h, yt, dyt)
      yout = y + h6 * (dydx + dyt + 2.0d0 * dym)
   end subroutine rk4

   !> 5th order Cash-Karp Runge-Kutta
   subroutine rkck(y, dydx, x, h, yout, yerr, derivs)
      implicit none
      real(dp), dimension(:), intent(in) :: y, dydx
      real(dp), intent(in) :: x, h
      real(dp), dimension(:), intent(out) :: yout, yerr
      interface
         subroutine derivs(x, y, dydx)
            use data_types
            implicit none
            real(dp), intent(in) :: x
            real(dp), dimension(:), intent(in) :: y
            real(dp), dimension(:), intent(out) :: dydx
         end subroutine derivs
      end interface
      integer :: ndum
      real(dp), dimension(size(y)) :: ak2, ak3, ak4, ak5, ak6, ytemp
      real(dp), parameter :: A2=0.2d0, A3=0.3d0, A4=0.6d0, A5=1.0d0, &
            A6=0.875d0, B21=0.2d0, B31=3.0d0/40d0, B32=9.0d0/40d0, B41=0.3d0, &
            B42=-0.9d0, B43=1.2d0, B51=-11.0d0/54d0, B52=2.5d0, B53=-70d0/27d0,&
            B54=35.0d0/27.0d0, B61=1631d0/55296d0, B62=175d0/512d0, &
            B63=575d0/13824d0, B64=44275d0/110592d0, B65=253d0/4096d0, &
            C1=37d0/378d0, C3=250d0/621d0, C4=125d0/594d0, C6=512d0/1771d0, &
            DC1=C1-2825d0/27648d0, DC3=C3-18575d0/48384d0, &
            DC4=C4-13525d0/55296d0, DC5=-277d0/14336d0, DC6=C6-0.25d0
      ndum = assert_eq((/ size(y), size(dydx), size(yout), size(yerr) /), 'rkck')
      ytemp = y + B21 * h * dydx
      call derivs(x + A2 * h, ytemp, ak2)
      ytemp = y + h * (B31 * dydx + B32 * ak2)
      call derivs(x + A3 * h, ytemp, ak3)
      ytemp = y + h * (B41 * dydx + B42 * ak2 + B43 * ak3)
      call derivs(x + A4 * h, ytemp, ak4)
      ytemp = y + h * (B51 * dydx + B52 * ak2 + B53 * ak3 + B54 * ak4)
      call derivs(x + A5 * h, ytemp, ak5)
      ytemp = y + h * (B61 * dydx + B62 * ak2 + B63 * ak3 + B64 * ak4 + B65 * ak5)
      call derivs(x + A6 * h, ytemp, ak6)
      yout = y + h * (C1 * dydx + C3 * ak3 + C4 * ak4 + C6 * ak6)
      yerr = h * (DC1 * dydx + DC3 * ak3 + DC4 * ak4 + DC5 * ak5 + DC6 * ak6)
   end subroutine rkck

   !returns 1 if value is nan
   function check_isnan_s(x) result(res)
     implicit none
     real(dp), intent(in) :: x
     integer :: res
     if ( isnan(x) ) then
       res = 1
     else
       res = 0
     end if
   end function check_isnan_s

   !returns 1 if any value is nan
   function check_isnan_v(x) result(res)
     implicit none
     real(dp),dimension(:), intent(in) :: x
     integer :: res,i
     res = 0
     do i = 1, size(x)
       if ( isnan(x(i)) ) then
         res = 1
         exit
       end if
     end do
   end function check_isnan_v

#if 0
   ! ====================================================================
   !  First derivative of the Bessel function of the first kind of order
   !  nu of z
   ! ====================================================================
   function dbesseljn(nu, z) result(dbj)
      implicit none
      integer, intent(in) :: nu
      real(dp), intent(in) :: z
      real(dp) :: dbj
      !dbj = fgsl_sf_bessel_jcn(nu - 1, z) - nu * fgsl_sf_bessel_jcn(nu, z) / z
      !dbj = dbesjn(nu - 1, z) - dble(nu) * dbesjn(nu, z) / z
      if ( dabs(z) < 1d-9 ) then
         dbj = 0.5d0 * (dbesjn(nu-1, z) - dbesjn(nu+1, z))
      else
         dbj = - dbesjn(nu + 1, z) + dble(nu) * dbesjn(nu, z) / z
      end if
   end function dbesseljn


#endif
end module misc
