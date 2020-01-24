module Aglow_models
   use data_types
   use constants
   use misc
   use SRtoolkit
   implicit none

contains

   !############################################################################
   !   #####  ######  #     #  #####   #####
   !  #     # #     # ##    # #     # #     #
   !  #       #     # # #   # #     # #     #
   !   #####  ######  #  #  #  ######  #####
   !        # #       #   # #       # #     #
   !  #     # #       #    ## #     # #     #
   !   #####  #       #     #  #####   #####
   !
   !----- Description -----
   !   Here are the analytical expressions to reproduce the results in
   !   Sari, Piran & Narayan, 1998, ApJ, 497, L17 (SPN98).
   !############################################################################

   !
   !     Evolution model of a blast wave as in eqs. (9) and (10) of SPN98
   !
   subroutine blastwave_approx(tobs, z, G0, E0, n, Gshk, Rshk, adiabatic)
      implicit none
      real(dp), intent(in) :: tobs, G0, E0, n, z
      logical, intent(in) :: adiabatic
      real(dp), intent(out) :: Rshk
      real(dp) :: M, L, Gshk, t
      t = tobs / (1d0 + z) ! <-- transforming to the comoving frame
      M = E0 / (G0 * cLight**2)
      L = (17d0 * M / (16d0 * pi * mass_p * n))**(1d0 / 3d0)
      if ( adiabatic ) then
         Rshk = (17d0 * E0 * t / (4d0 * pi * mass_p * n * cLight))**(0.25d0)
         Gshk = (17d0 * E0 / (1024d0 * pi * mass_p * n * cLight**5 * t**3))**(0.125d0)
      else
         Rshk = (4d0 * cLight * t / L)**(1d0 / 7d0) * L
         Gshk = (4d0 * cLight * t / L)**(-3d0 / 7d0)
      end if
      Gshk = dmin1(G0, Gshk)
   end subroutine blastwave_approx


   !
   !     Synchrotron spectra and light curves as in eqs. (11) SPN98
   !
   subroutine syn_aglow_SPN98(nuo, to, z, E0, epse, epsB, G0, pind, n, d_lum, adiab, flux)
      implicit none
      real(dp), intent(in) :: E0, epsB, epse, G0, n, d_lum, pind, z
      real(dp), intent(in), dimension(:) :: nuo, to
      logical :: adiab
      real(dp), intent(out), dimension(:,:) :: flux
      integer :: i, j
      real(dp) :: E52, d28, G2, Fmax, nu_c, nu_m, t0, tdy
      real(dp), dimension(size(nuo, dim=1)) :: nu
      real(dp), dimension(size(to, dim=1)) :: t

      E52 = E0 / 1d52
      d28 = d_lum / 1d28
      G2 = G0 / 100d0
      nu = nuo * (1d0 + z)
      t = to / (1d0 + z)
      ! -----[ Transition between slow and fast cooling ]-----
      if ( adiab ) then
         t0 = 210d0 * (epsB * epse)**2 * E52 * n
      else
         t0 = 4.6d0 * (epsB * epse)**1.4d0 * (E52 / G2)**0.8d0 * n**0.6d0
      end if

      do i=1, size(t)

         tdy = t(i) / 86400d0 ! <-- convert seconds to days

         ! -----[ Is the blast wave adiabatic or radiavtive? ]-----
         if ( adiab ) then
            nu_c = 2.7d12 * epsB**(-1.5d0) * (E52 * tdy)**(-0.5d0) / n
            nu_m = 5.7d14 * epse**2 * dsqrt(epsB * E52) * tdy**(-1.5d0)
            Fmax = 1.1d5 * E52 * dsqrt(epsB * n) / d28**2
         else
            nu_c = 1.3d13 * (G2 / E52)**(4d0 / 7d0) * epsB**(-1.5d0) * tdy**(-2d0 / 7d0) * n**(-13d0 / 14d0)
            nu_m = 1.2d14 * epse**2 * dsqrt(epsB) * (E52 / G2)**(4d0 / 7d0) * n**(-1d0 / 14d0) * tdy**(-12d0 / 7d0)
            Fmax = 4.5d3 * dsqrt(epsB) * (E52 / G2)**(8d0 / 7d0) * n**(5d0 / 14d0) * tdy**(-3d0 / 7d0) / d28**2
         end if

         do j = 1, size(nu)
            ! -----[ Are we fast or slow cooling? ]-----
            if ( tdy < t0 ) then
               if ( nu_c > nu(j) ) then
                  flux(j, i) = (nu(j) / nu_c)**(1d0 / 3d0) * Fmax
               else if ( nu_m >= nu(j) .and. nu(j) >= nu_c ) then
                  flux(j, i) = Fmax / dsqrt(nu(j) / nu_c)
               else
                  flux(j, i) = Fmax * (nu(j) / nu_m)**(-0.5d0 * pind) / dsqrt(nu_m / nu_c)
               end if
            else
               if ( nu_m > nu(j) ) then
                  flux(j, i) = (nu(j) / nu_m)**(1d0 / 3d0) * Fmax
               else if ( nu_c >= nu(j) .and. nu(j) >= nu_m ) then
                  flux(j, i) = Fmax * (nu(j) / nu_m)**(-0.5d0 * (pind - 1d0))
               else
                  flux(j, i) = Fmax * (nu(j) / nu_c)**(-0.5d0 * pind) * (nu_c / nu_m)**(-0.5d0 * (pind - 1d0))
               end if
            end if
         end do

      end do

   end subroutine syn_aglow_SPN98


   !############################################################################
   !  ######  #     #
   !  #     # #  #  #     ####   ####  #      #    # ##### #  ####  #    #
   !  #     # #  #  #    #      #    # #      #    #   #   # #    # ##   #
   !  ######  #  #  #     ####  #    # #      #    #   #   # #    # # #  #
   !  #     # #  #  #         # #    # #      #    #   #   # #    # #  # #
   !  #     # #  #  #    #    # #    # #      #    #   #   # #    # #   ##
   !  ######   ## ##      ####   ####  ######  ####    #   #  ####  #    #
   !
   !----- Description -----
   !   Here are the solutions and properties of the blast wave. The solutions
   !   are both analytic and numerical. Properties like the geometry, beaming,
   !   etc.
   !############################################################################

   !
   !     Deceleration radius as in eq. (1) of RM92
   !
   subroutine deceleration_radius(Rd1, Rd2, E0, G0, Aw, with_wind, s)
      implicit none
      real(dp), intent(in) :: G0, E0, Aw, s
      logical, intent(in) :: with_wind
      real(dp), intent(out) :: Rd1, Rd2
      if ( with_wind ) then
         !     Eq. (5) in PK00
         Rd1 = ((3d0 - s) * E0 / (4d0 * pi * Aw * mass_p * (cLight * G0)**2))**(1d0 / (3d0 - s))
         Rd2 = Rd1 * 0.25d0 ! see R_B in PVP14, p. 3
      else
         !     Eq. (1) in RM92
         Rd1 = (3d0 * E0 / (4d0 * pi * (G0 * cLight)**2 * mass_p * Aw))**(1d0 / 3d0)
         Rd2 = Rd1 * 2d0**(-2d0 / 3d0) ! see R_B in PVP14, p. 3
      end if
   end subroutine deceleration_radius


   !
   !     Analytic solution for the adiabatic blast wave
   !
   function adiab_blast_wave(Rshk, G0, E0, Aw, with_wind, s) result(Gshk)
      implicit none
      real(dp), intent(in) :: Rshk, G0, E0, Aw, s
      logical, intent(in) :: with_wind
      real(dp) :: M0, x, Gshk, R0
      if ( with_wind ) then
         !-----> Eqs. (4)-(5) in PK00
         R0 = ((3d0 - s) * E0 / (4d0 * pi * mass_p * cLight**2 * Aw * G0**2))**(1d0 / (3d0 - s))
         x = Rshk / R0
         Gshk = 0.5d0 * x**(s - 3d0) * G0 * ( dsqrt( 4d0 * x**(3d0 - s) + 1d0 + (2d0 * x**(3d0 - s) / G0)**2 ) - 1d0 )
      else
         !-----> Eqs. (9)-(10) in CD99
         M0 = E0 / (G0 * cLight**2)
         x = 4d0 * pi * mass_p * Aw * Rshk**3 / 3d0
         Gshk = (x + G0 * M0) / dsqrt(M0**2 + 2d0 * G0 * M0 * x + x**2)
      end if
   end function adiab_blast_wave


   !
   !     Blas wave cross sectional area
   !
   subroutine bw_crossec_area(G0, Rbw, Gbulk, theta_j0, beam_kind, blob, Rb, volume, csa, Oj)
      implicit none
      integer, intent(in) :: beam_kind
      real(dp), intent(in) :: Rbw, theta_j0, Gbulk, G0
      logical, intent(in) :: blob
      real(dp), intent(out) :: csa, volume, Rb, Oj
      real(dp) :: theta_j

      ! -----[ Uniform isotropic or beamed? ]-----
      if ( beam_kind >= 0 ) then

         select case( beam_kind )
         case(0)
            theta_j = 1d0 / G0
         case(1)
            theta_j = theta_j0
         case(2)
            theta_j = theta_j0 + 1d0 / (Gbulk * dsqrt(3d0))
         case(3)
            theta_j = theta_j0 + 1d0 / Gbulk
         case default
            call an_error("bw_crossec_area: wrong value of beam_kind")
         end select

         Oj = 2d0 * pi * (1d0 - dcos(theta_j))

         if ( blob ) then
            Rb = Rbw * theta_j
            volume = 4d0 * pi * Rb**3 / 3d0
            if ( beam_kind == 0 ) then
               csa = 2d0 * pi * Rb**2
            else
               csa = Oj * Rbw**2
            end if
         else
            Rb = Rbw / (Gbulk * 12d0)
            csa = Oj * Rbw**2
            volume = csa * Rb
         end if

      else

         Oj = 4d0 * pi
         ! Rb = Rbw / Gbulk
         Rb = Rbw / (Gbulk * 12d0)
         volume = 4d0 * pi * Rbw**2 * Rb
         csa = 4d0 * pi * Rbw**2

      end if

   end subroutine bw_crossec_area


   !############################################################################
   !    ##   #####  #   ##   #####      ####   ####   ####  #
   !   #  #  #    # #  #  #  #    #    #    # #    # #    # #
   !  #    # #    # # #    # #####     #      #    # #    # #
   !  ###### #    # # ###### #    #    #      #    # #    # #
   !  #    # #    # # #    # #    #    #    # #    # #    # #
   !  #    # #####  # #    # #####      ####   ####   ####  ######
   !
   !----- Description -----
   !   Here are the ways to estimate the adiabatic cooling.
   !############################################################################

   !
   !     Numerical form of the adiabatic term
   !
   function adiab_cool_num(vol1, vol2, dt) result(dotg_ad)
      implicit none
      real(dp), intent(in) :: vol1, vol2, dt
      real(dp) :: dotg_ad
      dotg_ad = (dlog(vol2) - dlog(vol1)) / (3d0 * dt)
   end function adiab_cool_num


end module Aglow_models
