module params
   use data_types
   use misc
   implicit none

   integer :: par_numbins, par_numdt, par_numdf, par_time_grid
   real(dp) :: par_R, par_R0, par_d_lum, par_z, par_gamma_bulk, par_theta_obs, &
      par_sigma, par_b_index, par_theta_e, par_zetae, par_L_j, par_eps_e, &
      par_tstep, par_tmax, par_tmin, par_eps_B, par_tvar, par_frec, par_qind, &
      par_g1, par_g2, par_gmin, par_gmax, par_nu_ext, par_uext, par_numin, &
      par_numax, par_E0, par_n_ext, par_B, par_mu_mag, par_eta_j
   public :: par_numbins, par_numdt, par_numdf, par_time_grid
   public :: par_R, par_R0, par_d_lum, par_z, par_gamma_bulk, par_theta_obs, &
      par_sigma, par_b_index, par_theta_e, par_zetae, par_L_j, par_eps_e, &
      par_tstep, par_tmax, par_tmin, par_eps_B, par_tvar, par_frec, par_qind, &
      par_g1, par_g2, par_gmin, par_gmax, par_nu_ext, par_uext, par_numin, &
      par_numax, par_E0, par_B, par_mu_mag, par_eta_j

contains

   subroutine read_params(par_file)
      implicit none
      character(*), intent(in) :: par_file
      integer :: ios
      open(unit=77, file=par_file, iostat=ios)
      if ( ios /= 0 ) call an_error("Paramo: Parameter file "//par_file//" could not be opened")
      read(77, *) par_R
      read(77, *) par_R0
      read(77, *) par_d_lum
      read(77, *) par_z
      read(77, *) par_theta_obs
      read(77, *) par_gamma_bulk
      read(77, *) par_mu_mag
      read(77, *) par_sigma
      read(77, *) par_frec
      read(77, *) par_b_index
      read(77, *) par_B
      read(77, *) par_theta_e
      read(77, *) par_zetae
      read(77, *) par_tstep
      read(77, *) par_tmax
      read(77, *) par_tmin
      read(77, *) par_tvar
      read(77, *) par_L_j
      read(77, *) par_eta_j
      read(77, *) par_E0
      read(77, *) par_n_ext
      read(77, *) par_eps_e
      read(77, *) par_eps_B
      read(77, *) par_g1
      read(77, *) par_g2
      read(77, *) par_gmin
      read(77, *) par_gmax
      read(77, *) par_qind
      read(77, *) par_nu_ext
      read(77, *) par_uext
      read(77, *) par_numin
      read(77, *) par_numax
      read(77, *) par_numbins
      read(77, *) par_numdt
      read(77, *) par_numdf
      read(77, *) par_time_grid
      close(77)
   end subroutine read_params

end module params
