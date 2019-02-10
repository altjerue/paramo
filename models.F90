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

   subroutine adiab_blast_wave(t, z, G0, E0, n, Rd, Gshk, Rshk)
      ! ************************************************************************
      !  Description:
      !     This is the setup for the model in Sari, Piran & Narayan, 1998,
      !     ApJ, 497, L17.
      ! ************************************************************************
      implicit none
      real(dp), intent(in) :: t, z, G0, E0, n, Rd
      real(dp), intent(out) :: Gshk, Rshk
      real(dp) :: P0, B0, Pshk, td!, M0, LSedov
      ! M0 = E0 / (G0 * cLight**2)
      B0 = bofg(G0)
      P0 = pofg(G0)
      ! lSedov = Rd * G0**(2d0 / 3d0)
      td = (1d0 + z) * Rd / (B0 * cLight * G0**2)
      if ( t <= td ) then
         Rshk = Rd * t / td
         Pshk = P0
      else if ( t > td .and. t <= td * G0**(8d0 / 3d0) ) then
         Rshk = Rd * (2d0 * t / td)**(1d0 / 4d0)
         Pshk = P0 * (2d0 * t / td)**(-3d0 / 8d0) / sqrt2
      else
         Rshk = Rd * (5d0 * t / (sqrt2**3 * G0 * td))**(2d0 / 5d0)
         Pshk = P0 * (5d0 * t / (sqrt2**3 * G0 * td))**(-3d0 / 5d0) / sqrt2
      end if
      Gshk = dsqrt(Pshk**2 + 1d0)
   end subroutine adiab_blast_wave

end module models
