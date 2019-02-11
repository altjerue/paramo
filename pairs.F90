module pairs
  use data_types
  use constants
  use misc
  implicit none

  interface Rgg_MK95
     module procedure Rgg_MK95_s
     module procedure Rgg_MK95_v
  end interface Rgg_MK95

contains

   !
   ! Reaction rate fit as in Mastichiadis & Kirk (1995)
   !
   function Rgg_MK95_s(w) result(Rgg)
      implicit none
      real(dp), intent(in) :: w
      real(dp) :: Rgg
      if ( w - 1d0 > 0d0 ) then
         Rgg = 0.652d0 * (w**2 - 1d0) * dlog(w) / w**3
      else
         Rgg = 0d0
      end if
   end function Rgg_MK95_s

   function Rgg_MK95_v(w) result(Rgg)
      implicit none
      real(dp), intent(in), dimension(:) :: w
      integer :: i
      real(dp), dimension(size(w)) :: Rgg
      do i = 1, size(w)
         if ( w(i) - 1d0 > 0d0 ) then
            Rgg(i) = 0.652d0 * (w(i)**2 - 1d0) * dlog(w(i)) / w(i)**3
         else
            Rgg(i) = 0d0
         end if
      end do
   end function Rgg_MK95_v


   subroutine pairs_loss_rate(nu, n, Lgg)
      implicit none
      real(dp), intent(in), dimension(:) :: nu, n
      integer :: i, Nf
      real(dp), dimension(size(nu)) :: x, Lgg
      Nf = size(nu)
      x = hPlanck * nu / (mass_e * cLight**2)
      do i = 1, Nf
         Lgg(i) = x(1) * n(1) * Rgg_MK95(x(i) * x(1)) + x(Nf) * n(Nf) * Rgg_MK95(x(i) * x(Nf))
         Lgg(i) = 0.5d0 * dlog(x(Nf) / x(1)) * (Lgg(i) + 2d0 * &
            sum(x(2:Nf - 1) * n(2:Nf - 1) * Rgg_MK95(x(i) * x(2:Nf)))) / &
            dble(Nf - 1)
         Lgg(i) = n(i) * Lgg(i)
      end do
   end subroutine pairs_loss_rate


   subroutine electron_pairs_inj(g, nu, n, Qgg)
      implicit none
      real(dp), intent(in), dimension(:) :: nu, n, g
      integer :: k, Nf, Ng
      real(dp) :: nn, dn
      real(dp), dimension(size(g)) :: x, Qgg
      real(dp), dimension(size(nu)) :: xp
      Nf = size(nu)
      Ng = size(g)
      x = 2d0 * g
      xp = hPlanck * nu / (mass_e * cLight**2)
      gloop: do k = 1, Ng
         if ( x(k) > nu(Nf) .or. x(k) < nu(1) ) then
            Qgg(k) = 0d0
            cycle gloop
         end if
         call polint(nu, n, x(k), nn, dn)
         Qgg(k) = xp(1) * n(1) * Rgg_MK95(x(k) * xp(1)) + xp(Nf) * n(Nf) * Rgg_MK95(x(k) * xp(Nf))
         Qgg(k) = 0.5d0 * dlog(xp(Nf) / xp(1)) * ( Qgg(k) + 2d0 * &
            sum(xp(2:Nf - 1) * n(2:Nf - 1) * Rgg_MK95(x(k) * xp(2:Nf))) ) / &
            dble(Nf - 1)
         Qgg(k) = 4d0 * nn * Qgg(k)
      end do gloop
   end subroutine electron_pairs_inj


end module pairs
