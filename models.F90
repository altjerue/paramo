module models
   use data_types
   use constants
   use SRtoolkit
   implicit none

contains

   function deceleration_radius(E0, G0, n) result(Rd)
      implicit none
      real(dp), intent(in) :: G0, E0, n
      real(dp) :: Rd
      Rd = (3d0 * E0 / (4d0 * pi * (G0 * cLight)**2 * mass_p * n))**(1d0 / 3d0)
   end function deceleration_radius

   subroutine adiab_blast_wave(Rshk, R0, G0, E0, n, Gshk)
      ! ************************************************************************
      !  Description:
      !     This is the setup for the blast wave in Petropoulou & Mastichiadis,
      !     2009, A&A, 507, 599
      ! ************************************************************************
      implicit none
      real(dp), intent(in) :: R0, Rshk, G0, E0, n
      real(dp), intent(out) :: Gshk
      real(dp) :: M0, l, x
      M0 = E0 / (G0 * cLight**2)
      x = Rshk / R0
      l = 4d0 * pi * mass_p * n * R0**3 / (3d0 * M0)
      Gshk = (l * (x**3 - 1d0) + G0) / dsqrt( 1d0 + 2d0 * G0 * l * (x**3 - 1d0) + (l * (x**3 - 1d0))**2 )
   end subroutine adiab_blast_wave

   subroutine blastwave_approx(t, G0, E0, n, Rshk, Gshk, adiabatic)
      ! ************************************************************************
      !  Description:
      !     This is the setup for the model in Sari, Piran & Narayan, 1998,
      !     ApJ, 497, L17.
      ! ************************************************************************
      implicit none
      real(dp), intent(in) :: t, G0, E0, n
      logical, intent(in) :: adiabatic
      real(dp), intent(out) :: Rshk, Gshk
      real(dp) :: M, L
      M = E0 / (G0 * cLight**2)
      L = (17d0 * M / (16d0 * pi * mass_p * n))**(1d0 / 3d0)
      if ( adiabatic ) then
         Rshk = (17d0 * E0 * t / (4d0 * pi * mass_p * n * cLight))**(0.25d0)
         Gshk = (17d0 * E0 / (1024d0 * pi * mass_p * n * cLight**5 * t**3))**(0.125d0)
      else
         Rshk = (4d0 * cLight * t / L)**(1d0 / 7d0) * L
         Gshk = (4d0 * cLight * t / L)**(-3d0 / 7d0)
      end if
   end subroutine blastwave_approx

end module models
