module models
   use data_types
   use constants
   implicit none
contains
   !  #####  ######  #     #  #####   #####
   ! #     # #     # ##    # #     # #     #
   ! #       #     # # #   # #     # #     #
   !  #####  ######  #  #  #  ######  #####
   !       # #       #   # #       # #     #
   ! #     # #       #    ## #     # #     #
   !  #####  #       #     #  #####   #####
   subroutine shock_afterglow(t, G0, E0, eps_e, eps_B, n, B, Gshock, Rshock, n_bs, ue_bs, adiab)
      ! ************************************************************************
      !  Description:
      !     This is the setup for the model in Sari, Piran & Narayan, 1998,
      !     ApJ, 497, L17.
      ! ************************************************************************
      implicit none
      real(dp), intent(in) :: eps_e, eps_B, t, G0, E0, n
      real(dp), intent(out) :: B, Gshock, Rshock, n_bs, ue_bs
      logical, optional, value :: adiab
      real(dp) :: M0,L

      if ( .not. present(adiab) ) adiab = .true.

      M0 = E0 / (G0 * cLight**2)
      L = ( 17d0 * M0 / (16d0 * pi * mass_p * n) )**(1d0 / 3d0)
   
      Gshock = shock_Lorentz(t, E0, G0, n, L, adiab)
      Rshock = shock_radius(t, E0, n, L, adiab)
      B = dsqrt(32d0 * pi * mass_p * eps_B * n) * Gshock * cLight
      
      n_bs = 4d0 * Gshock * n
      ue_bs = 4d0 * Gshock**2 * n * mass_p * cLight**2

   contains

      ! function shock_energy(Gsh, L, n, ad) result(E)
      !    implicit none
      !    real(dp), intent(in) :: Gsh, L, n
      !    logical, intent(in) :: ad
      !    if ( ad ) then
      !       E =  16d0 * pi * Gsh**2 * R**3 * n * mp * cspeed**2 / 17d0
      !    else
      !       E = 16d0 * pi * Gsh * L**3 * n * mp * cspeed**2 / 17d0
      !    end if
      ! end function shock_energy

      function shock_radius(tt, E, nn, LL, ad) result(R)
         implicit none
         real(dp), intent(in) :: E, nn, tt, LL
         logical, intent(in) :: ad
         real(dp) :: R
         if ( ad ) then
            R = ( 17d0 * E * tt / (4d0 * pi * mass_p * nn * cLight) )**0.25
         else
            R = ( 4d0 * cLight * tt / LL)**(1d0 / 7d0) * LL
         end if
      end function shock_radius


      function shock_Lorentz(tt, E, G, nn, LL, ad) result(Gsh)
         implicit none
         real(dp), intent(in) :: E, G, nn, tt, LL
         logical, intent(in) :: ad
         real(dp) :: Gsh

         if ( tt <= 0d0 ) then
            Gsh = G
            return
         end if

         if ( ad ) then
            Gsh = dmax1(1d0, ( 17d0 * E / (1024d0 * pi * mass_p * nn * cLight**5 * tt**3) )**0.125)
         else
            Gsh = dmax1(1d0, ( 4d0 * cLight * tt / LL)**(-3d0 / 7d0))
         end if
      end function shock_Lorentz

   end subroutine shock_afterglow

end module models
