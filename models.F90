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
      !     This is the setup for the model in Sari, Piran & Narayan, 1998,
      !     ApJ, 497, L17.
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

   subroutine adiab_blast_wave_approx(t, Rd, G0, E0, n, z, Rshk, Gshk)
      ! ************************************************************************
      !  Description:
      !     This is the setup for the model in Sari, Piran & Narayan, 1998,
      !     ApJ, 497, L17.
      ! ************************************************************************
      implicit none
      real(dp), intent(in) :: t,  z, G0, E0, n, Rd
      real(dp), intent(out) :: Rshk, Gshk
      Rshk = (17d0 * t * E0 / (4d0 * pi * mass_p * cLight * n * (1d0 + z)))**(1d0 / 4d0)
      Gshk = dmin1(G0, (17d0 * E0 * (1d0 + z)**3 / (1024d0 * pi * mass_p * cLight**5 * n * t**3))**(1d0 / 8d0))
      Gshk = dmax1(1.0001d0, Gshk)
   end subroutine adiab_blast_wave_approx

end module models
