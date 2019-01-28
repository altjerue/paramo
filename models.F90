module models
   use data_types
   use constants
   use SRtoolkit
   implicit none

contains

   subroutine adiab_blast_wave(t, G0, E0, n, Gshock, Rshock)
      ! ************************************************************************
      !  Description:
      !     This is the setup for the model in Sari, Piran & Narayan, 1998,
      !     ApJ, 497, L17.
      ! ************************************************************************
      implicit none
      real(dp), intent(in) :: t, G0, E0, n
      real(dp), intent(out) :: Gshock, Rshock
      real(dp) :: M0, Rd, P0, B0, Pshock, td
      ! l = 4d0 * pi * rho0 * R0**3 / (3d0 * M0)
      M0 = E0 / (G0 * cLight**2)
      B0 = bofg(G0)
      P0 = G0 * B0
      Rd = ( 3d0 * E0 / (4d0 * pi * (G0 * cLight)**2 * mass_p * n) )**(1d0 / 3d0)
      td = Rd / (B0 * cLight * G0**2)
      Rshock = shock_radius(t, td, Rd, G0)
      Pshock = shock_momentum(t, td, P0, G0)
      Gshock = dsqrt(Pshock**2 + 1d0)

   contains

      function shock_radius(tt, ttd, xxd, gg) result(Rsh)
         implicit none
         real(dp), intent(in) :: tt, ttd, xxd, gg
         real(dp) :: Rsh
         if ( tt <= ttd ) then
            Rsh = xxd * tt / ttd
         else if ( tt > ttd .and. tt <= ttd * gg**(8d0 / 3d0) ) then
            Rsh = xxd * (2d0 * tt / td)**0.25d0
         else
            Rsh = xxd * (5d0 * tt / (dsqrt(8d0) * gg * ttd))**0.4d0
         end if
      end function shock_radius

      function shock_momentum(tt, ttd, pp, gg) result(Psh)
         implicit none
         real(dp), intent(in) :: tt, ttd, pp, gg
         real(dp) :: Psh
         if ( tt <= ttd ) then
            Psh = pp
         else if ( tt > ttd .and. tt <= ttd * gg**(8d0 / 3d0) ) then
            Psh = pp * (2d0 * tt / td)**(-0.375d0) / dsqrt(2d0)
         else
            Psh = pp * (5d0 * tt / (dsqrt(8d0) * gg * ttd))**(-0.6d0) / dsqrt(2d0)
         end if
      end function shock_momentum

   end subroutine adiab_blast_wave

end module models
