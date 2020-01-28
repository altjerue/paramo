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
!
! *************************************************************************
module misc
   use data_types
   implicit none

   interface realloc
      module procedure realloc_1d
      module procedure realloc_2d
   end interface

contains

   subroutine realloc_1d(array, size)
      implicit none
      integer, intent(in) :: size
      real(dp), intent(inout), allocatable, dimension(:) :: array
      if ( allocated(array) ) deallocate(array)
      allocate(array(size))
   end subroutine realloc_1d

   subroutine realloc_2d(array, size1, size2)
      implicit none
      integer, intent(in) :: size1, size2
      real(dp), intent(inout), allocatable, dimension(:, :) :: array
      if ( allocated(array) ) deallocate(array)
      allocate(array(size1, size2))
   end subroutine realloc_2d


   !
   !  ::::  Polinomial interpolation  ::::
   !
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
         if ( any( dabs(den(1:n - m)) == 0d0 ) ) then !&
            print*,x,xa
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


   !
   !     :::::   Linear interpolation   :::::
   !
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
   if(dy1 + dy2 > dy3) then
      write(*, *) y1, y2, y
      call an_error('linint: Interpolation not between y1 and y2')
   end if
   end subroutine linint


   ! :::: This produces an error message ::::
   subroutine an_error(string)
      implicit none
      character(len=*), intent(IN) :: string
      write (*, "('ERROR ', A)") string
      stop 'program terminated by an_error'
   end subroutine an_error


   !       This produces a warning message
   subroutine a_warning(string)
      implicit none
      character(len=*), intent(in) :: string
      write(*,*) 'warging message: '//trim(string)
   end subroutine a_warning


   !      This creates an arithmetic progression array
   function arth(first, increment, n)
      implicit none
      integer, parameter :: NPAR_ARTH=16,NPAR2_ARTH=8
      real(dp), intent(in) :: first,increment
      integer, intent(in) :: n
      real(dp), dimension(n) :: arth
      integer :: k,k2
      real(dp) :: temp
      if (n > 0) arth(1)=first
      if (n <= NPAR_ARTH) then
         do k=2,n
            arth(k)=arth(k-1)+increment
         end do
      else
         do k=2,NPAR2_ARTH
            arth(k)=arth(k-1)+increment
         end do
         temp=increment*NPAR2_ARTH
         k=NPAR2_ARTH
         do
            if (k >= n) exit
            k2=k+k
            arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
            temp=temp+temp
            k=k2
         end do
      end if
   end function arth


   !
   !     Find closest element in an array to a value
   !
   function locate(xx, x, in_bounds)
      implicit none
      real(dp), intent(in) :: x
      real(dp), intent(in), dimension(:) :: xx
      logical, intent(in), optional :: in_bounds
      integer :: locate
      integer :: n, jl, jm, ju
      logical :: ascnd, bounds
      !! NOTE: The flag 'bounds' tells the function to return the index of the
      !!       lower value of the interval of 'xx' in which 'x' is
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


   !
   !   -----{  Chebychev evaluation  }-----
   !
   function chebev(x,coef,num_coefs,xmin,xmax) result(res)
      implicit none
      integer, intent(in) :: num_coefs
      real(dp), intent(in) :: x,xmin,xmax
      real(dp), intent(in), dimension(:) :: coef
      integer :: j
      real(dp) :: d,dd,y,y2,sv,res
      if ((x-xmin)*(x-xmax) > 0d0) then
         print*,x,xmin,xmax
         write(*,*)'x is not in rage in chebev'
         stop
      end if
      d = 0d0
      dd = 0d0
      y = (2d0 * x - xmin - xmax) / (xmax - xmin)
      y2 = 2d0 * y
      do j=num_coefs,2,-1
         sv = dd
         dd = d
         d = y2 * dd - sv + coef(j)
      end do
      res = y * d - dd + 0.5d0 * coef(1)
   end function chebev


   !
   !   -----{  Tridiagonal matrix solver  }-----
   !
   subroutine tridag_ser(a, b, c, r, u)
      !  Description:
      !     Solver of a tridiagonal matrix. Based on the code in
      !     "Numberical Recipes".
      !
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

      do j = 2, n
         gam(j) = c(j - 1) / bet
         bet = b(j) - a(j - 1) * gam(j)
         if ( bet == 0.0d0 ) call an_error('tridag_ser: error at code stage 2')
         u(j) = (r(j) - a(j - 1) * u(j - 1)) / bet
      end do

      do j = n - 1, 1, -1
         u(j) = u(j) - gam(j + 1) * u(j + 1)
      end do

   end subroutine tridag_ser


   !
   !   -----{  Trapezoidal rule  }-----
   !
   subroutine trapzd(func, a, b, s, n)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: a, b
      real(dp), intent(inout) :: s
      interface
         function func(x) result(res)
            use data_types
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(size(x, dim=1)) :: res
         end function func
      end interface
      real(dp) :: del, fsum
      integer :: it
      if (n == 1) then
         s = 0.5d0 * (b - a) * sum( func((/ a, b /)) )
      else
         it = 2**(n - 2)
         del = (b - a) / it
         fsum = sum( func(arth(a + 0.5d0 * del, del, it)) )
         s = 0.5d0 * (s + del * fsum)
      end if
   end subroutine trapzd


   !
   !   -----{  checking equality of integers  }-----
   !
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


   !                              #####
   !  ####  #    #   ##   #####  #     # # #    # #####
   ! #    # #    #  #  #  #    #       # # ##   #   #
   ! #      ###### #    # #    #  #####  # # #  #   #
   ! #      #    # ###### #####  #       # #  # #   #
   ! #    # #    # #    # #   #  #       # #   ##   #
   !  ####  #    # #    # #    # ####### # #    #   #
   !
   function char2int(c) result(i)
      !  Description:
      !     Convert a string to a numeric value using an internal read.
      !
      character(len=*), intent(in) :: c
      integer :: i
      read (c,'(I5)') i
   end function char2int

   function char2double(c) result(d)
      character(len=*), intent(in) :: c
      real(dp) :: d
      read(c, *) d
   end function char2double

   function zeros1D(n) result(a)
      implicit none
      integer, intent(in) :: n
      real(dp), dimension(n) :: a
      a = 0d0
   end function zeros1D

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
