module blastwave
   use data_types
   use constants
   use transformers
   use misc
   use SRtoolkit
   implicit none

contains

   !   #####  ######  #     #  #####   #####
   !  #     # #     # ##    # #     # #     #
   !  #       #     # # #   # #     # #     #
   !   #####  ######  #  #  #  ######  #####
   !        # #       #   # #       # #     #
   !  #     # #       #    ## #     # #     #
   !   #####  #       #     #  #####   #####
   !
   !> Evolution model of a blast-wave as in eqs. (9) and (10) of Sari, Piran &
   !! Narayan (1998)
   subroutine blastwave_approx_SPN98(G0, E0, n, tobs, Gshk, Rshk, adiabatic)
      implicit none
      real(dp), intent(in) :: G0, E0, n, tobs
      logical, intent(in) :: adiabatic
      real(dp), intent(out) :: Rshk, Gshk
      real(dp) :: M, L
      M = E0 / (G0 * cLight**2)
      L = (17d0 * M / (16d0 * pi * mass_p * n))**(1d0 / 3d0)
      adiab: if ( adiabatic ) then
         Rshk = (17d0 * E0 * tobs / (fourpi * mass_p * n * cLight))**(0.25d0)
         Gshk = (17d0 * E0 / (1024d0 * pi * mass_p * n * cLight**5 * tobs**3))**(0.125d0)
      else
         Rshk = (4d0 * cLight * tobs / L)**(1d0 / 7d0) * L
         Gshk = (4d0 * cLight * tobs / L)**(-3d0 / 7d0)
      end if adiab
      Gshk = dmin1(G0, Gshk)
   end subroutine blastwave_approx_SPN98

   !> Synchrotron spectra and light curves as in eqs. (11) of Sari, Piran &
   !! Narayan (1998)
   subroutine syn_afterglow_SPN98(nu_o, t_o, E0, epse, epsB, G0, pind, n, d_lum, adiab, flux)
      implicit none
      real(dp), intent(in) :: E0, epsB, epse, G0, n, d_lum, pind, t_o, nu_o
      logical :: adiab
      real(dp), intent(out) :: flux
      real(dp) :: E52, d28, G2, Fmax, nu_c, nu_m, t0, tdy

      E52 = E0 / 1d52
      d28 = d_lum / 1d28
      G2 = G0 / 100d0

      !   Transition between slow and fast cooling
      if ( adiab ) then
         t0 = 210d0 * (epsB * epse)**2 * E52 * n
      else
         t0 = 4.6d0 * (epsB * epse)**1.4d0 * (E52 / G2)**0.8d0 * n**0.6d0
      end if

      tdy = sec2dy(t_o)

      blastwave: if ( adiab ) then
         nu_c = 2.7d12 * epsB**(-1.5d0) * (E52 * tdy)**(-0.5d0) / n
         nu_m = 5.7d14 * epse**2 * dsqrt(epsB * E52) * tdy**(-1.5d0)
         Fmax = 1.1d5 * E52 * dsqrt(epsB * n) / d28**2
      else
         nu_c = 1.3d13 * (G2 / E52)**(4d0 / 7d0) * epsB**(-1.5d0) * tdy**(-2d0 / 7d0) * n**(-13d0 / 14d0)
         nu_m = 1.2d14 * epse**2 * dsqrt(epsB) * (E52 / G2)**(4d0 / 7d0) * n**(-1d0 / 14d0) * tdy**(-12d0 / 7d0)
         Fmax = 4.5d3 * dsqrt(epsB) * (E52 / G2)**(8d0 / 7d0) * n**(5d0 / 14d0) * tdy**(-3d0 / 7d0) / d28**2
      end if blastwave

      cooling_fast_or_slow: if ( tdy < t0 ) then
         if ( nu_c > nu_o ) then
            flux = (nu_o / nu_c)**(1d0 / 3d0) * Fmax
         else if ( nu_m >= nu_o .and. nu_o >= nu_c ) then
            flux = Fmax / dsqrt(nu_o / nu_c)
         else
            flux = Fmax * (nu_o / nu_m)**(-0.5d0 * pind) / dsqrt(nu_m / nu_c)
         end if
      else
         if ( nu_m > nu_o ) then
            flux = (nu_o / nu_m)**(1d0 / 3d0) * Fmax
         else if ( nu_c >= nu_o .and. nu_o >= nu_m ) then
            flux = Fmax * (nu_o / nu_m)**(-0.5d0 * (pind - 1d0))
         else
            flux = Fmax * (nu_o / nu_c)**(-0.5d0 * pind) * (nu_c / nu_m)**(-0.5d0 * (pind - 1d0))
         end if
      end if cooling_fast_or_slow

   end subroutine syn_afterglow_SPN98
   !############################################################################


   !> Blast-wave solution. Deceleration radius as in eq. (1) of RM92
   subroutine deceleration_radius(Rd1, Rd2, E0, G0, Aw, with_wind, s)
      implicit none
      real(dp), intent(in) :: G0, E0, Aw
      logical, intent(in) :: with_wind
      real(dp), intent(in), optional :: s
      real(dp), intent(out) :: Rd1, Rd2
      if ( with_wind .and. .not. present(s) ) &
            call an_error("deceleration_radius: Wind index s not declared")
      if ( with_wind ) then
         !     Eq. (5) in PK00
         Rd1 = ( (3d0 - s) * E0 / (4d0 * pi * Aw * mass_p * (cLight * G0)**2) )**(1d0 / (3d0 - s))
         Rd2 = Rd1 * 0.25d0 ! see R_B in PVP14, p. 3
      else
         !     Eq. (1) in RM92
         Rd1 = ( 3d0 * E0 / (4d0 * pi * (G0 * cLight)**2 * mass_p * Aw) )**(1d0 / 3d0)
         Rd2 = Rd1 * 2d0**(-2d0 / 3d0) ! see R_B in PVP14, p. 3
      end if
   end subroutine deceleration_radius

   !> Adiabatic blast wave models.
   function adiab_blastwave(Rshk, sol_kind, G0, E0, Aw, s, Npts, filename) result(Gshk)
      implicit none
      real(dp), intent(in) :: Rshk
      integer, intent(in) :: sol_kind
      integer, intent(in), optional :: Npts
      character(len=*), intent(in), optional :: filename
      real(dp), intent(in), optional :: G0, E0, Aw, s
      real(dp) :: M0, x, Gshk, R0
      real(dp), dimension(:), allocatable :: rshknum, gshknum
      integer :: i

      select case( sol_kind )
      case(1)
      !---> Eqs. (9)-(10) in CD99
      M0 = E0 / (G0 * cLight**2)
         x = 4d0 * pi * mass_p * Aw * Rshk**3 / 3d0
      Gshk = (x + G0 * M0) / dsqrt(M0**2 + 2d0 * G0 * M0 * x + x**2)
      case(2)
         if ( .not. present(s) ) call an_error("deceleration_radius: Wind index s not declared")
        !---> Eqs. (4)-(5) in PK00
         R0 = ( (3d0 - s) * E0 / (4d0 * pi * mass_p * cLight**2 * Aw * G0**2))**(1d0 / (3d0 - s) )
         x = Rshk / R0
         Gshk = 0.5d0 * x**(s - 3d0) * G0 * ( dsqrt( 4d0 * x**(3d0 - s) + 1d0 + (2d0 * x**(3d0 - s) / G0)**2 ) - 1d0 )
      case(3)
         ! Numerical interpolation for a blast wave using an external file (ASCII)
         if ( .not. (present(Npts) .or. present(filename)) ) then
            call an_error("deceleration_radius: Arguments Npts or filename not declared")
         else
            call realloc(rshknum, Npts)
            call realloc(gshknum, Npts)
         end if
         open(444, file=trim(filename), status='old')
         do i = 1, Npts
            read(444,*) rshknum(i), gshknum(i)
         end do
         close(444)

         i = 1
         do while ( Rshk > rshknum(i) )
            i = 1 + i
         end do

         ! Interpolate
         call linint(rshknum(i-1), rshknum(i), Rshk, gshknum(i-1), gshknum(i), Gshk)
         !!TODO: Try polint. You can save the do while above, and have a more accurate value than linint.
         ! call polint(rshknum, gshknum, Rshk, Gshk, dR)
      case default
         call an_error("adiab_blastwave: wrong value of sol_kind")
      end select

   end function adiab_blastwave


   !> Depending on the model, the blast wave cross sectional area may be
   !! isotropic or beamed. This subroutine returns the the cross sectional area,
   !! volume, radius/thickness and Omega_j of the emitting region. The emitting
   !! region may be a blob or a slab.
   subroutine bw_crossec_area(beam_kind, blob, Rbw, Gbulk, theta_j0, Rb, volume, csa, Oj)
      implicit none
      integer, intent(in) :: beam_kind
      real(dp), intent(in) :: Rbw, theta_j0, Gbulk
      logical, intent(in) :: blob
      real(dp), intent(out) :: csa, volume, Rb, Oj
      real(dp) :: theta_j
      !---> Uniform isotropic or beamed?
      select case( beam_kind )
         case(0)!> Isotropic blast-wave
            theta_j = pi
         case(1)!> Half blob
            theta_j = halfpi
         case(2)!> Classic beamed jet
            theta_j = 1d0 / Gbulk
         case(3)!> Beamed jet with initial opening angle theta_j0
            theta_j = theta_j0
         case(4)
            theta_j = theta_j0 + 1d0 / (Gbulk * dsqrt(3d0))
         case(5)
            theta_j = theta_j0 + 1d0 / Gbulk
         case default
            call an_error("bw_crossec_area: wrong value of beam_kind")
      end select
      Oj = 2d0 * pi * (1d0 - dcos(theta_j))
      csa = Oj * Rbw**2
      if ( blob ) then
         Rb = Rbw * theta_j
         volume = 4d0 * pi * Rb**3 / 3d0
      else
         Rb = Rbw / (12d0 * Gbulk)
         volume = csa * Rb
      end if
   end subroutine bw_crossec_area

#if 0
   !> Blast-wave following Bianco & Ruffini (2005)
   subroutine bw_solver(Rin, Gin, Ein, Min, Nin, Gout, Eout, Mout, Nout, dR)
      implicit none
      real(dp), intent(in) :: Rin, Gin, Ein, Min, Nin, dR
      real(dp), intent(out) :: Gout, Eout, Mout, Nout
      real(dp), dimension(1, 4) :: yerr, dydx
      real(dp), dimension(2, 4) :: y
      y(1, 1) = Nin ! ISM mass
      y(1, 2) = Ein ! internal energy
      y(1, 3) = Gin ! bulk Lorentz factor
      y(1, 4) = Min ! mass-energy
      call bw_derivs(Rin, y(1, :), dydx)
      call rkck(y(1, :), dydx, Rin, dR, y(2, :), yerr, bw_derivs)
   contains
      !> Drivatives for blast-wave following Bianco & Ruffini (2005)
      subroutine bw_derivs(x, y, dydx)
         implicit none
         real(dp), intent(in) :: x
         real(dp), dimension(:), intent(in) :: y
         real(dp), dimension(:), intent(inout) :: dydx
         dydx(1) = fourpi * mass_p * n * x**2                    ! rate of swept up ISM mass
         dydx(2) = (y(3) - 1d0) * dydx(1) * cLight**2            ! internal energy deriv
         dydx(3) = - ((y(3)**2 - 1d0) / y(4)) * dydx(1)          ! bulk Lorentz factor deriv
         dydx(4) = ((1d0 - eps) * dydx(2) / cLight**2) + dydx(1) ! mass-energy deriv
      end subroutine bw_derivs
   end subroutine bw_solver


   !  #####                                      #
   ! #     # ##### #####  #    #  ####           # ###### #####
   ! #         #   #    # #    # #    #          # #        #
   !  #####    #   #    # #    # #               # #####    #
   !       #   #   #####  #    # #         #     # #        #
   ! #     #   #   #   #  #    # #    #    #     # #        #
   !  #####    #   #    #  ####   ####      #####  ######   #

   !> Read a blast-wave properties from numerical simulation shock profile
   !! @param filename input ASCII file
   !! @param r output radius of the shock
   !! @param theta output direction of the shock
   !! @param Gbulk output bulk Lorentz factor
   !! @param nlines output total number of directions
   subroutine bw_mezcal(filename, th_los, l_los, r, th, Gbulk, rho, mu_obs_r, mu_obs_v)
      implicit none
      real(dp), intent(in) :: th_los, l_los
      character(len=*), intent(in) :: filename
      real(dp), intent(out), dimension(:) :: r, th, Gbulk, rho, mu_obs_r, mu_obs_v
      integer :: i, io, nlines
      real(dp) :: v, vr, vh, th_v, l_v
      nlines = count_lines(filename) - 1
      ! call realloc(r, nlines)
      ! call realloc(v, nlines)
      ! call realloc(th, nlines)
      ! call realloc(Gbulk, nlines)
      ! allocate(r(nlines), v(nlines), theta(nlines), Gbulk(nlines))
      open(77, file=trim(filename), status='old', action='read')
      read(77, *)
      !!!TODO: Set the correct vectors, angles and trigonometry in general for the observer
      do i=1, nlines
         !! Reading columns
         read(77, *) r(i), th(i), vr, vh, rho(i)
         r(i) = r(i) * cLight
         !! Calculating the bulk Lorentz factor
         v = dsqrt(vr**2 + vh**2)
         Gbulk(i) = gofb(v)
         Gsh(i) = Gbulk(i) * sqrt2
         !! Calculating the observing viewing angle
         if ( th(i) < th_los ) then
            th_v = halfpi - dabs(th_los - th(i))
         else if ( th(i) == th_los ) then
            th_v = th_los
         else
            th_v = th(i) - th_los
         end if
         l_v = dsqrt(l_los**2 + r(i)**2 - l_los * r(i) * (dsin(th(i)) * dsin(th_los) + dcos(th(i)) * dcos(th_los)))
         mu_obs_r(i) = dsin(th(i)) * dsin(th_v) + dcos(th(i)) * dcos(th_v)
         mu_obs_v(i) = (vr * l_v) / (v * dsqrt(l_v**2 + th_v**2))
      end do
      close(77)
   end subroutine bw_mezcal
#endif

end module blastwave
