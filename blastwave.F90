module blastwave
   use data_types
   use constants
   use transformers
   use misc
   use SRtoolkit
   implicit none

   !> blast-wave type
   !! @param der course [derrotero]
   !! @param t_com time (comoving)
   !! @param r_lab position (lab)
   !! @param Gbulk bulk Lorentz factor
   !! @param cs cross-sectional area
   !! @param vol volume of the shocked region
   type blast_wave
      real(dp) :: der, vol, cs
      real(dp), allocatable, dimension(:) :: t, r, Gbulk
   end type blast_wave

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
         Rshk = (17d0 * E0 * tobs / (4d0 * pi * mass_p * n * cLight))**(0.25d0)
         Gshk = (17d0 * E0 / (1024d0 * pi * mass_p * n * cLight**5 * tobs**3))**(0.125d0)
      else
         Rshk = (4d0 * cLight * tobs / L)**(1d0 / 7d0) * L
         Gshk = (4d0 * cLight * tobs / L)**(-3d0 / 7d0)
      end if adiab
      Gshk = dmin1(G0, Gshk)
   end subroutine blastwave_approx_SPN98

   !> Synchrotron spectra and light curves as in eqs. (11) of Sari, Piran &
   !! Narayan (1998)
   subroutine syn_afterglow_SPN98(nuo, to, z, E0, epse, epsB, G0, pind, n, d_lum, adiab, flux)
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

      !   Transition between slow and fast cooling
      if ( adiab ) then
         t0 = 210d0 * (epsB * epse)**2 * E52 * n
      else
         t0 = 4.6d0 * (epsB * epse)**1.4d0 * (E52 / G2)**0.8d0 * n**0.6d0
      end if

      evolution: do i = 1, size(t)

         tdy = sec2dy(t(i))

         blastwave: if ( adiab ) then
            nu_c = 2.7d12 * epsB**(-1.5d0) * (E52 * tdy)**(-0.5d0) / n
            nu_m = 5.7d14 * epse**2 * dsqrt(epsB * E52) * tdy**(-1.5d0)
            Fmax = 1.1d5 * E52 * dsqrt(epsB * n) / d28**2
         else
            nu_c = 1.3d13 * (G2 / E52)**(4d0 / 7d0) * epsB**(-1.5d0) * tdy**(-2d0 / 7d0) * n**(-13d0 / 14d0)
            nu_m = 1.2d14 * epse**2 * dsqrt(epsB) * (E52 / G2)**(4d0 / 7d0) * n**(-1d0 / 14d0) * tdy**(-12d0 / 7d0)
            Fmax = 4.5d3 * dsqrt(epsB) * (E52 / G2)**(8d0 / 7d0) * n**(5d0 / 14d0) * tdy**(-3d0 / 7d0) / d28**2
         end if blastwave

         spectrum: do j = 1, size(nu)
            cooling_fast_or_slow: if ( tdy < t0 ) then
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
            end if cooling_fast_or_slow
         end do spectrum

      end do evolution

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


   !> Analytic solution for the adiabatic blast wave.
   function adiab_blast_wave(Rshk, G0, E0, Aw, with_wind, s) result(Gshk)
      implicit none
      real(dp), intent(in) :: Rshk, G0, E0, Aw, s
      logical, intent(in) :: with_wind
      real(dp) :: M0, x, Gshk, R0
      if ( with_wind ) then
         !---> Eqs. (4)-(5) in PK00
         R0 = ( (3d0 - s) * E0 / (4d0 * pi * mass_p * cLight**2 * Aw * G0**2))**(1d0 / (3d0 - s) )
         x = Rshk / R0
         Gshk = 0.5d0 * x**(s - 3d0) * G0 * ( dsqrt( 4d0 * x**(3d0 - s) + 1d0 + (2d0 * x**(3d0 - s) / G0)**2 ) - 1d0 )
      else
         !---> Eqs. (9)-(10) in CD99
         M0 = E0 / (G0 * cLight**2)
         x = 4d0 * pi * mass_p * Aw * Rshk**3 / 3d0
         Gshk = (x + G0 * M0) / dsqrt(M0**2 + 2d0 * G0 * M0 * x + x**2)
      end if
   end function adiab_blast_wave


   !> Depending on the model, the blast wave cross sectional area may be 
   !! isotropic or beamed. This subroutine returns the the cross sectional area,
   !! volume, radius/thickness and Omega_j of the emitting region. The emitting
   !! region may be a blob or a slab.
   subroutine bw_crossec_area(beam_kind, blob, G0, Rbw, Gbulk, theta_j0, Rb, volume, csa, Oj)
      implicit none
      integer, intent(in)   :: beam_kind
      real(dp), intent(in)  :: Rbw, theta_j0, Gbulk, G0
      logical, intent(in)   :: blob
      real(dp), intent(out) :: csa, volume, Rb, Oj
      real(dp)              :: theta_j

      !---> Uniform isotropic or beamed?
      iso_or_beamed: if ( beam_kind >= 0 ) then

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
            ! Rb = Rbw / (Gbulk * 12d0)
            Rb = Rbw / (12d0 * (Gbulk + 0.75d0))
            csa = Oj * Rbw**2
            volume = csa * Rb
         end if

      else

         !--->  Isotropic spherical blast-wave
         Oj = 4d0 * pi
         ! Rb = Rbw / Gbulk
         ! Rb = Rbw / (Gbulk * 12d0)
         Rb = Rbw / (12d0 * (Gbulk + 0.75d0))
         volume = 4d0 * pi * Rbw**2 * Rb
         csa = 4d0 * pi * Rbw**2

      end if iso_or_beamed

   end subroutine bw_crossec_area


   !  #####                                      #              
   ! #     # ##### #####  #    #  ####           # ###### ##### 
   ! #         #   #    # #    # #    #          # #        #   
   !  #####    #   #    # #    # #               # #####    #   
   !       #   #   #####  #    # #         #     # #        #   
   ! #     #   #   #   #  #    # #    #    #     # #        #   
   !  #####    #   #    #  ####   ####      #####  ######   #
   !
   !> Read a blast-wave properties from numerical simulation shock profile
   !! @param filename input ASCII file
   !! @param r output radius of the shock
   !! @param theta output direction of the shock
   !! @param Gbulk output bulk Lorentz factor
   !! @param nlines output total number of directions
   subroutine bw_mezcal(filename, nlines, r, theta, Gbulk)
      implicit none
      integer, intent(in) :: nlines
      character(len=*), intent(in) :: filename
      real(dp), intent(out), allocatable, dimension(:) :: r, theta, Gbulk
      integer :: i, io
      real(dp) :: x, y, vx,vy
      ! real(dp), allocatable, dimension(:) :: v
      if ( nlines /= count_lines(filename) ) call an_error("bw_mezcal: nlines and number of lines in "//trim(filename)//"are not the same")
      call realloc(r, nlines)
      ! call realloc(v, nlines)
      call realloc(theta, nlines)
      call realloc(Gbulk, nlines)
      ! allocate(r(nlines), v(nlines), theta(nlines), Gbulk(nlines))
      open(77, file=trim(filename), iostat=io, status='old', action='read')
      if (io /= 0) stop "Cannot open file!"
      do i=1, nlines
         read(77, *) x, y, vx, vy, theta(i), Gbulk(i)
         r(i) = dsqrt(x**2 + y**2)
         ! v(i) = dsqrt(vx**2 + vy**2)
      end do
      close(77)
   end subroutine bw_mezcal

end module blastwave
