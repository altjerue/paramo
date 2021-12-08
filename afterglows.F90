#define ADIABATIC_NONE     (0)
#define ADIABATIC_VOL_EVOL (1)
#define ADIABATIC_MSB00    (2)
#define ADIABATIC_COOLING_TYPE (ADIABATIC_MSB00)

!!!NOTE
!!! To keep indices clean the order is:
!!!

! #####  #        ##    ####  #####       #    #   ##   #    # ######
! #    # #       #  #  #        #         #    #  #  #  #    # #
! #####  #      #    #  ####    #   ##### #    # #    # #    # #####
! #    # #      ######      #   #         # ## # ###### #    # #
! #    # #      #    # #    #   #         ##  ## #    #  #  #  #
! #####  ###### #    #  ####    #         #    # #    #   ##   ######
!
!> 1D blast-wave from the self-similar solution
subroutine bw1D_afterglow(params_file, output_file, with_wind, cool_withKN, blob)
   use data_types
   use constants
   use params
   use misc
   use pwl_integ
#ifdef HDF5
   use hdf5
   use h5_inout
#endif
   use SRtoolkit
   use distribs
   use radiation
   ! use pairs
   use blastwave
   use specialf
   implicit none

   character(len=*), intent(in) :: params_file
   logical, intent(in) :: cool_withKN, with_wind, blob
   character(len=*),intent(inout) :: output_file
   integer, parameter :: nmod = 50
   character(len=*), parameter :: screan_head = &
      '---------------------------------------------------------------------'&
      //new_line('A')//&
      ' | Iteration | Obser. time |   BW radius |  gamma_bulk |      Bfield |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------', &
      on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
#ifdef HDF5
   integer(HID_T) :: file_id, group_id
   integer :: herror
#endif
   integer :: i, j, k, numbins, numdf, numt, time_grid, flow_kind, bw_sol_type, Npts
   real(dp) :: uB, uext, L_j, gmin, gmax, numin, numax, pind, B, R0, Rmax, &
         tinj, g1, g2, tstep, Q0, tmax, d_lum, z, n_ext, urad_const, Aw, sind, &
         theta_obs, mu_obs, nu_ext, tesc_e, uext0, eps_e, tlc, g1_const, Rd2, &
         Ejet, eps_B, E0, gamma_bulk0, L_e, nu_ext0, tmin, td, Rd, dr, &
         b_const , eps_g2, theta_j0, cs_area, n_ext0, g2_const, &
         Omega_j, dt, f_esc
   real(dp), allocatable, dimension(:) :: freqs, t, Inu, gg, urad, beta_bulk, &
         nu_obs, t_obs, gamma_bulk, R, Dopp, tcool, Rb, volume, dotg_tmp
   real(dp), allocatable, dimension(:,:) :: dotg, n_e, jnut, jmbs, jssc, jeic, &
         ambs, anut, Qinj, tau_gg, pow_syn, Ddiff
   logical :: bw_approx, radius_evol, pwl_over_trpzd_integ, &
         with_ic, ssa_boiler, numerical


   !  #####    ##   #####    ##   #    #  ####
   !  #    #  #  #  #    #  #  #  ##  ## #
   !  #    # #    # #    # #    # # ## #  ####
   !  #####  ###### #####  ###### #    #      #
   !  #      #    # #   #  #    # #    # #    #
   !  #      #    # #    # #    # #    #  ####
   call read_params(params_file)
   eps_e = par_eps_e
   eps_B = par_eps_B
   eps_g2 = par_eps_acc
   gamma_bulk0 = par_gamma_bulk
   d_lum = par_d_lum
   z = par_z
   tstep = par_tstep
   tmin = par_tmin
   gmin = par_gmin
   gmax = par_gmax
   pind = par_pind
   numin = par_numin
   numax = par_numax
   nu_ext0 = par_nu_ext
   uext0 = par_uext
   numbins = par_NG
   numt = par_NT
   numdf = par_NF
   n_ext0 = par_n_ext
   g1 = par_g1
   g2 = par_g2
   time_grid = par_time_grid
   R0 = par_R0
   Rmax = par_R
   tmax = par_tmax
   E0 = par_E0
   f_esc = par_fesc


   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   write(*, "('--> Simulation setup')")

   allocate(t(0:numt), freqs(numdf), Inu(numdf), gg(numbins), urad(numbins), &
         nu_obs(numdf), t_obs(0:numt), R(0:numt), tcool(numbins), Rb(0:numt), &
         gamma_bulk(0:numt), beta_bulk(0:numt), Dopp(0:numt), volume(0:numt), &
         dotg_tmp(numbins))
   allocate(n_e(numbins, 0:numt), Ddiff(numbins, 0:numt), &
         dotg(numbins, 0:numt), ambs(numdf, numt), jmbs(numdf, numt), &
         jnut(numdf, numt), jssc(numdf, numt), anut(numdf, numt), &
         jeic(numdf, numt), Qinj(numbins, 0:numt), pow_syn(numdf, numbins), &
         tau_gg(numdf, numt))

   call K1_init
   call K2_init

   !!!!!TODO: Reduce hardcodedness. Transform all these into arguments of the subroutine
   with_ic = .true.
   bw_approx = .false.
   radius_evol = .false.
   pwl_over_trpzd_integ = .false.
   flow_kind = 0 !1
   ssa_boiler = .false.
   numerical = .false.
   Npts = 200

   !-----> About the observer
   theta_obs = par_theta_obs * pi / 180d0
   mu_obs = dcos(theta_obs)

   !-----> Initializing blast wave
   theta_j0 = 0.2d0
   beta_bulk(0) = bofg(gamma_bulk0)
   call bw_crossec_area(flow_kind, blob, R0, gamma_bulk0, theta_j0, Rb(0), volume(0), cs_area, Omega_j)

   !-----> External medioum
   if ( with_wind ) then
      !!!NOTE: n_ext0 is A_{*}
      sind = 2d0
      Aw = n_ext0 * 3d35
      n_ext = Aw * R0**(-sind)
      bw_sol_type = 2
   else
      sind = 0d0
      Aw = n_ext0
      n_ext = n_ext0
      bw_sol_type = 1
   end if
   if ( numerical ) bw_sol_type = 3
   call deceleration_radius(Rd, Rd2, E0, gamma_bulk0, Aw, with_wind, sind)
   td = (1d0 + z) * Rd / (4d0 * gamma_bulk0**2 * cLight)

   !---> True outflow energy and decceleration radius
   Ejet = E0 * Omega_j / (4d0 * pi)

   !---> Locating the emission region
   if ( bw_approx ) then
      t_obs(0) = tstep
      call blastwave_approx_SPN98(gamma_bulk0, E0, Aw, t_obs(0), gamma_bulk(0), R(0), .true.)
      beta_bulk(0) = bofg(gamma_bulk(0))
      Dopp(0) = Doppler(gamma_bulk(0), mu_obs)
      t(0) = R(0) / (beta_bulk(0) * gamma_bulk(0) * cLight)
   else
      R(0) = R0
      if ( numerical ) then
         gamma_bulk(0) = adiab_blastwave(R(0), bw_sol_type, Npts=Npts, &
               filename='/home/lcombi/Dropbox/repository/astrograv_tools/comala_new/paramo/numerical_gshk.dat')
      else
         gamma_bulk(0) = adiab_blastwave(R(0), bw_sol_type, G0=gamma_bulk0, E0=E0, Aw=Aw, s=sind)
      end if
      beta_bulk(0) = bofg(gamma_bulk(0))
      Dopp(0) = Doppler(gamma_bulk(0), mu_obs)
      t_obs(0) = (1d0 + z) * R(0) / (beta_bulk(0) * gamma_bulk(0) * cLight * Dopp(0))
      t(0) = R(0) / (beta_bulk(0) * gamma_bulk(0) * cLight)
   end if

   call bw_crossec_area(flow_kind, blob, R(0), gamma_bulk(0), theta_j0, Rb(0), volume(0), cs_area, Omega_j)

   !---> External medioum
   n_ext = Aw * R(0)**(-sind)

   !---> Magnetic field
   b_const = dsqrt(32d0 * pi * eps_B * mass_p) * cLight
   !B = b_const * dsqrt(n_ext) * gamma_bulk(0)
   !B = b_const * dsqrt(n_ext * (gamma_bulk(0) - 1d0) * gamma_bulk(0))
   B = b_const * dsqrt(n_ext * (gamma_bulk(0) - 1d0) * (gamma_bulk(0) + 0.75d0))
   uB = B**2 / (8d0 * pi)

   !---> Radiation fields
   uext = uext0 * gamma_bulk(0)**2 * (1d0 + beta_bulk(0)**2 / 3d0) ! eq. (5.25) in DM09
   nu_ext = nu_ext0 * gamma_bulk(0)
   urad_const = 4d0 * sigmaT * cLight / (3d0 * energy_e)

   !---> Minimum and maximum Lorentz factors of the particles distribution
   g2_const = dsqrt(6d0 * pi * eCharge * eps_g2 / sigmaT)
   g1_const = eps_e * mass_p * (pind - 2d0) / ((pind - 1d0) * mass_e)
   g2 = g2_const / dsqrt(B)
   g1 = g1_const * (gamma_bulk(0) - 1d0)

   !---> Time-scales
   tlc = Rb(0) / cLight
   tesc_e = f_esc * tlc! * R(0) / (cLight * gamma_bulk(0))!
   tinj = 1d200

   !---> Fraction of accreted kinetic energy injected into non-thermal electrons
   L_e = eps_e * cs_area * n_ext * energy_p * cLight * beta_bulk(0) * gamma_bulk(0) * (gamma_bulk(0) - 1d0)
   Q0 = L_e / ( g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) * volume(0) * energy_e )
   !Q0 = L_e / ((g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) &
   !      - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6)) * volume(0) * energy_e)

   write(*, "('mu_obs  =', ES15.7)") mu_obs
   write(*, "('Omega_j =', ES15.7)") Omega_j
   write(*, "('volume  =', ES15.7)") volume(0)
   write(*, "('Gamma_0 =', ES15.7)") gamma_bulk(0)
   write(*, "('Aw      =', ES15.7)") Aw
   write(*, "('sind    =', ES15.7)") sind
   write(*, "('E_0     =', ES15.7)") E0
   write(*, "('E_jet   =', ES15.7)") Ejet
   write(*, "('L_e     =', ES15.7)") L_e
   write(*, "('Q_0     =', ES15.7)") Q0
   write(*, "('pind    =', ES15.7)") pind
   write(*, "('B0      =', ES15.7)") B
   write(*, "('u_B     =', ES15.7)") uB
   write(*, "('u_ext   =', ES15.7)") uext
   write(*, "('nu_ext  =', ES15.7)") nu_ext
   write(*, "('R_0     =', ES15.7)") R(0)
   write(*, "('t_0     =', ES15.7)") t(0)
   write(*, "('t_obs   =', ES15.7)") t_obs(0)
   write(*, "('R_d     =', ES15.7)") Rd
   write(*, "('Doppler =', ES15.7)") Dopp(0)
   write(*, "('gamma_1 =', ES15.7)") g1
   write(*, "('gamma_2 =', ES15.7)") g2

   build_f: do j = 1, numdf
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numdf - 1)) )
   end do build_f

   build_g: do k = 1, numbins
      gg(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numbins - 1))
      ! gg(k) = (gmin - 1d0) * ((gmax - 1d0) / (gmin - 1d0))**(dble(k - 1) / dble(numbins - 1)) + 1d0
   end do build_g

   dotg(:, 0) = urad_const * (uB + uext) * pofg(gg)**2
   Inu = 0d0
   Ddiff(:, 0) = 1d-200
   Qinj(:, 0) = injection_pwl(t(0), tinj, gg, g1, g2, pind, Q0)
   n_e(:, 0) = Qinj(:, 0)

   write(*, *) ''
   write(*, "(' Using tstep = ', F7.3)") tstep
   write(*, *) "Wrting data in: ", trim(output_file)
   write(*, *) ''
   write(*, *) screan_head


   ! ###### #    #  ####  #      #    # ##### #  ####  #    #
   ! #      #    # #    # #      #    #   #   # #    # ##   #
   ! #####  #    # #    # #      #    #   #   # #    # # #  #
   ! #      #    # #    # #      #    #   #   # #    # #  # #
   ! #       #  #  #    # #      #    #   #   # #    # #   ##
   ! ######   ##    ####  ######  ####    #   #  ####  #    #
   time_loop: do i = 1, numt

      !-----> Localizing the emission region
      if ( bw_approx ) then
         t_obs(i) = t_obs(0) * (tmax / t_obs(0))**(dble(i) / dble(numt))
         call blastwave_approx_SPN98(gamma_bulk0, E0, Aw, t_obs(i), gamma_bulk(i), R(i), .true.)
         dr = R(i) - R(i - 1)
         Dopp(i) = Doppler(gamma_bulk(i), mu_obs)
         beta_bulk(i) = bofg(gamma_bulk(i))
         call rk2_arr(t(i - 1), 1d0 / (beta_bulk(i - 1:i) * gamma_bulk(i - 1:i) * cLight), dr, t(i))
         dt = t(i) - t(i - 1)
      else
         R(i) = R0 * (Rmax / R0)**(dble(i) / dble(numt))
         dr = R(i) - R(i - 1)
         if ( numerical ) then
            gamma_bulk(i) = adiab_blastwave(R(i), bw_sol_type, Npts=Npts, &
                  filename='/home/lcombi/Dropbox/repository/astrograv_tools/comala_new/paramo/numerical_gshk.dat')
         else
            gamma_bulk(i) = adiab_blastwave(R(i), bw_sol_type, G0=gamma_bulk0, E0=E0, Aw=Aw, s=sind)
         end if
         beta_bulk(i) = bofg(gamma_bulk(i))
         Dopp(i) = Doppler(gamma_bulk(i), mu_obs)
         call rk2_arr(t(i - 1), 1d0 / (beta_bulk(i - 1:i) * gamma_bulk(i - 1:i) * cLight), dr, t(i))
         dt = t(i) - t(i - 1)
      end if
      call rk2_arr(t_obs(i - 1), (1d0 + z) / Dopp(i - 1:i), dt, t_obs(i))
      Rb(i) = R(i) / gamma_bulk(i)

      !  ###### ###### #####
      !  #      #      #    #
      !  #####  #####  #    #
      !  #      #      #    #
      !  #      #      #    #
      !  ###### ###### #####
      call FP_FinDif_difu(dt, &
            &             gg, &
            &             n_e(:, i - 1), &
            &             n_e(:, i), &
            &             dotg(:, i - 1), &
            &             Ddiff(:, i - 1), &
            &             Qinj(:, i - 1), &
            &             1d200, &
            &             tlc)

      call bw_crossec_area(flow_kind, blob, R(i), gamma_bulk(i), theta_j0, Rb(i), volume(i), cs_area, Omega_j)

      !-----> External medioum
      n_ext = Aw * R(i)**(-sind)

      !-----> Magnetic field assuming equipartition
      !B = b_const * dsqrt(n_ext) * gamma_bulk(i)
      !B = b_const * dsqrt(n_ext * (gamma_bulk(i) - 1d0) * gamma_bulk(i))
      B = b_const * dsqrt(n_ext * (gamma_bulk(i) - 1d0) * (gamma_bulk(i) + 0.75d0))
      uB = B**2 / (8d0 * pi)

      !-----> Radiation fields
      uext = uext0 * gamma_bulk(i)**2 * (1d0 + beta_bulk(i)**2 / 3d0) ! eq. (5.25) in DM09
      nu_ext = nu_ext0 * gamma_bulk(i)

      !-----> Minimum and maximum Lorentz factors of the particles distribution
      g2 = g2_const / dsqrt(B)
      g1 = g1_const * (gamma_bulk(i) - 1d0)

      !-----> Time-scales
      tlc = Rb(i) / cLight
      !!!!!COMBAK: here may be implemented escape time-scale in PVP14, eq. (27)
      ! call bolometric_integ(freqs, opt_depth(anut(:, i - 1), 2d0 * Rb(i)), tau)
      ! tesc_e = tlc * 0.75d0 * Rb(i) * (1d0 + (1d0 - dexp(-tau)) / (1d0 + dexp(-tau))) / cLight
      !!!!!!!!!!!!!!!
      tesc_e = f_esc * tlc! * R(i) / (cLight * gamma_bulk(i))!
      tinj = 1d200

      !-----> Fraction of accreted kinetic energy injected into non-thermal electrons
      L_e = eps_e * cs_area * n_ext * energy_p * cLight * beta_bulk(i) * gamma_bulk(i) * (gamma_bulk(i) - 1d0)
      Q0 = L_e / ( g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) * volume(i) * energy_e )
      !!!!!NOTE: The expression below corresponds to the normalization in Eq. (13)
      !!!!!      in PM09
      !Q0 = L_e / ((g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) - &
      !      g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6)) * volume(i) * energy_e)
      Qinj(:, i) = injection_pwl(t(i), tinj, gg, g1, g2, pind, Q0)


      !  #####    ##   #####  #   ##   ##### #  ####  #    #
      !  #    #  #  #  #    # #  #  #    #   # #    # ##   #
      !  #    # #    # #    # # #    #   #   # #    # # #  #
      !  #####  ###### #    # # ######   #   # #    # #  # #
      !  #   #  #    # #    # # #    #   #   # #    # #   ##
      !  #    # #    # #####  # #    #   #   #  ####  #    #
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j, k)
      do j = 1, numdf
         freqs(j) = nu_com_f(nu_obs(j), z, Dopp(i))
         call syn_emissivity(jmbs(j, i), freqs(j), gg, n_e(:, i), B)
         call syn_absorption(ambs(j, i), freqs(j), gg, n_e(:, i), B)
      end do
      !$OMP END PARALLEL DO

      !!!!!TODO: pairs optical depth
      ! if ( hPlanck * freqs(j) > 5d11 * (0.01d0 / 0.02d0) * (100d0 / gamma_bulk(i)) ) then
      !    tau_gg(j, i) = 0.16d0 * (R(i) / 1d16) * (100d0 / gamma_bulk(i)) * (uext / 1d-7) * (1d11 / (hPlanck * freqs(j))) * (0.01d0 / 0.02d0)**2
      ! else
      !    tau_gg(j, i) = 0d0
      ! end if

      anut(:, i) = ambs(:, i)! + tau_gg(j, i) / (2d0 * Rb(i))

      call RadTrans(Inu, Rb(i), jmbs(:, i), ambs(:, i))

      !-----> Synchrotron boiler from GGS88
      if ( ssa_boiler ) then
         !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j, k)
         do j = 1, numdf
            do k = 1, numbins
               ! pow_syn(j, k) = dsqrt(3d0) * eCharge**3 * B * RMA_new(freqs(j) / (nuconst * B), gg(k))
               !-----> Expression below is Eq. (3.63) in my thesis
               pow_syn(j, k) = 1.315d-28 * nuconst * B * volume(i) * RMA_new(freqs(j) / (nuconst * B), gg(k))
            end do
         end do
         !$OMP END PARALLEL DO
         do k=1,numbins
            call bolometric_integ(freqs, Inu * pow_syn(:, k) / freqs**2, Ddiff(k, i))
            Ddiff(k, i) = Ddiff(k, i) * gg(k) * pofg(gg(k)) / (2d0 * mass_e * energy_e)
         end do
      else
         Ddiff(:, i) = 1d-200
      end if

      !-----> Inverse Compton emissivity
      if ( with_ic ) then
         !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j)
         do j = 1, numdf
            ! call IC_emis_full(freqs(j), freqs, gg, n_e(:, i), Inu, jssc(j, i))
            call IC_iso_powlaw(jssc(j, i), freqs(j), freqs, Inu, n_e(:, i), gg)
            call IC_iso_monochrom(jeic(j, i), freqs(j), uext, nu_ext, n_e(:, i), gg)
            jnut(j, i) = jmbs(j, i) + jssc(j, i) + jeic(j, i)
         end do
         !$OMP END PARALLEL DO
      else
         jnut(:, i) = jmbs(:, i)
      end if


      !   ####   ####   ####  #      # #    #  ####
      !  #    # #    # #    # #      # ##   # #    #
      !  #      #    # #    # #      # # #  # #
      !  #      #    # #    # #      # #  # # #  ###
      !  #    # #    # #    # #      # #   ## #    #
      !   ####   ####   ####  ###### # #    #  ####

      !> Synchrotron cooling
      dotg(:, i) = urad_const * uB * pofg(gg)**2 &
      !> Adiabatic cooling
#if( ADIABATIC_COOLING_TYPE == ADIABATIC_MSB00 )
         + cLight * beta_bulk(i) * gamma_bulk(i) * gg / R(i)
      !!!!!NOTE: using time-scile in Hao's paper, eq. (11)
      !dotg(:, i) = dotg(:, i) + 1.6d0 * cLight * gamma_bulk(i) * pofg(gg) / R(i)
#elif( ADIABATIC_COOLING_TYPE == ADIABATIC_VOL_EVOL )
         + pofg(gg) * dlog(volume(l, i) / volume(l, i - 1)) / (3d0 * dt)
#else
         + 0d0
#endif
      ! call bolometric_integ(freqs, 4d0 * pi * Inu / cLight, urad)
      ! call RadTrans_blob(Inu, R, jssc(:, i) + jeic(:, i), anut(:, i))
      call RadTrans_blob(Inu, Rb(i), jmbs(:, i), ambs(:, i))
      call rad_cool_pwl(dotg_tmp, gg, freqs, 4d0 * pi * Inu / cLight, cool_withKN)
      dotg(:, i) = dotg(:, i) + dotg_tmp
      call rad_cool_mono(dotg_tmp, gg, nu_ext, uext, cool_withKN)
      dotg(:, i) = dotg(:, i) + dotg_tmp

      !   ####  #    #     ####   ####  #####  ###### ###### #    #
      !  #    # ##   #    #      #    # #    # #      #      ##   #
      !  #    # # #  #     ####  #      #    # #####  #####  # #  #
      !  #    # #  # #         # #      #####  #      #      #  # #
      !  #    # #   ##    #    # #    # #   #  #      #      #   ##
      !   ####  #    #     ####   ####  #    # ###### ###### #    #
      if ( mod(i, nmod) == 0 .or. i == 1 ) &
            write(*, on_screen) i, t_obs(i), R(i), gamma_bulk(i), B

   end do time_loop
   write(*, "(' ---------------------------------------------------------------------')")

   !  ####    ##   #    # # #    #  ####
   ! #       #  #  #    # # ##   # #    #
   !  ####  #    # #    # # # #  # #
   !      # ###### #    # # #  # # #  ###
   ! #    # #    #  #  #  # #   ## #    #
   !  ####  #    #   ##   # #    #  ####
   write(*, "('--> Saving')")
#ifdef HDF5
   write(*, "('--> Opening HDF5')")
   ! ------  Opening output file  ------
   call h5open_f(herror)
   call h5io_createf(output_file, file_id, herror)

   ! ------  Saving initial parameters  ------
   call h5io_createg(file_id, "Parameters", group_id, herror)
   call h5io_wint0(group_id, 'numt', numt, herror)
   call h5io_wint0(group_id, 'numdf', numdf, herror)
   call h5io_wint0(group_id, 'numbins', numbins, herror)
   call h5io_wint0(group_id, 'time-grid', time_grid, herror)
   call h5io_wdble0(group_id, 't_max', tmax, herror)
   call h5io_wdble0(group_id, 't_min', tmin, herror)
   call h5io_wdble0(group_id, 'tstep', tstep, herror)
   call h5io_wdble0(group_id, 'd_lum', d_lum, herror)
   call h5io_wdble0(group_id, 'redshift', z, herror)
   call h5io_wdble0(group_id, 'Gamma_bulk0', gamma_bulk0, herror)
   call h5io_wdble0(group_id, 'view-angle', par_theta_obs, herror)
   call h5io_wdble0(group_id, 'gamma_min', gmin, herror)
   call h5io_wdble0(group_id, 'gamma_max', gmax, herror)
   call h5io_wdble0(group_id, 'gamma_1', g1, herror)
   call h5io_wdble0(group_id, 'gamma_2', g2, herror)
   call h5io_wdble0(group_id, 'pwl-index', pind, herror)
   call h5io_wdble0(group_id, 'L_j', L_j, herror)
   call h5io_wdble0(group_id, 'epsilon_e', eps_e, herror)
   call h5io_wdble0(group_id, 'epsilon_B', eps_B, herror)
   call h5io_wdble0(group_id, 'nu_min', numin, herror)
   call h5io_wdble0(group_id, 'nu_max', numax, herror)
   call h5io_wdble0(group_id, 'nu_ext0', nu_ext0, herror)
   call h5io_wdble0(group_id, 'u_ext0', uext0, herror)
   call h5io_wdble0(group_id, 'E0', E0, herror)
   call h5io_wdble0(group_id, 'Ejet', Ejet, herror)
   call h5io_wdble0(group_id, 'R0', R0, herror)
   call h5io_wdble0(group_id, 'n_ext', n_ext0, herror)
   call h5io_closeg(group_id, herror)

   ! ------  Saving data  ------
   call h5io_wdble0(file_id, 'Rd', Rd, herror)
   call h5io_wdble0(file_id, 'td', td, herror)
   call h5io_wdble1(file_id, 'time', t(1:), herror)
   call h5io_wdble1(file_id, 't_obs', t_obs(1:), herror)
   call h5io_wdble1(file_id, 'Rb', Rb(1:), herror)
   call h5io_wdble1(file_id, 'Rbw', R(1:), herror)
   call h5io_wdble1(file_id, 'volume', volume(1:), herror)
   call h5io_wdble1(file_id, 'Gamma_bulk', gamma_bulk(1:), herror)
   call h5io_wdble1(file_id, 'Doppler', Dopp(1:), herror)
   call h5io_wdble1(file_id, 'nu', freqs, herror)
   call h5io_wdble1(file_id, 'nu_obs', nu_obs, herror)
   call h5io_wdble1(file_id, 'gamma', gg, herror)
   call h5io_wdble2(file_id, 'jnut', jnut, herror)
   call h5io_wdble2(file_id, 'jmbs', jmbs, herror)
   call h5io_wdble2(file_id, 'jssc', jssc, herror)
   call h5io_wdble2(file_id, 'jeic', jeic, herror)
   call h5io_wdble2(file_id, 'anut', anut, herror)
   call h5io_wdble2(file_id, 'ambs', ambs, herror)
   call h5io_wdble2(file_id, 'Qinj', Qinj, herror)
   call h5io_wdble2(file_id, 'n_e', n_e(:, 1:), herror)
   call h5io_wdble2(file_id, 'cool-coef', dotg(:, 1:), herror)
   call h5io_wdble2(file_id, 'diffusion', Ddiff(:, 1:), herror)

   ! ------  Closing output file  ------
   call h5io_closef(file_id, herror)
   call h5close_f(herror)
#endif
   write(*, "('==========  FINISHED  ==========')")
   write(*,*) ''

end subroutine bw1D_afterglow


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if 0
! #    # ###### ######  ####    ##   #
! ##  ## #          #  #    #  #  #  #
! # ## # #####     #   #      #    # #
! #    # #        #    #      ###### #
! #    # #       #     #    # #    # #
! #    # ###### ######  ####  #    # ######
!
!> Blast-wave from mezcal simulations.
!! @param params_file input file with list of parameters
!! @output_file output name of HDF5 file
!! @param KNcool input boolean radiative cooling with or without Klein-Nishina
!! @param assume_blob input boolean geometry of the emission region: .true. for blob, .false. for slab
subroutine mezcal(params_file, output_file, with_ic, KNcool, assume_blob)
   use data_types
   use constants
   use params
   use misc
   use pwl_integ
#ifdef HDF5
   use hdf5
   use h5_inout
#endif
   use SRtoolkit
   use distribs
   use radiation
   use blastwave
   use specialf
   implicit none
   character(len=*), intent(in) :: params_file
   logical, intent(in) :: assume_blob, KNcool, with_ic
   character(len=*),intent(inout) :: output_file
   integer, parameter :: nmod=50
   character(len=*), parameter :: screan_head = &
      '| Iteration | Obser. time |   BW radius |  gamma_bulk |      Bfield |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------',&
      on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
#ifdef HDF5
   integer(HID_T) :: file_id, group_id
   integer :: herror
#endif
   integer :: i, j, k, l, numt, numf, numg, numd
   real(dp) :: uB, uext, L_j, gmin, gmax, numin, numax, pind, B, R0, Rmax, &
         tinj, g1, g2, tstep, Q0, tmax, d_lum, z, n_ext, urad_const, Aw, sind, &
         theta_obs, mu_obs, nu_ext, tesc_e, uext0, eps_e, tlc, g1_const, Rd2, &
         Ejet, eps_B, E0, gamma_bulk0, L_e, nu_ext0, tmin, td, Rd, dr, &
         b_const, beta_bulk, eps_g2, theta_j0, cs_area, n_ext0, g2_const, &
         Omega_j, dt, f_esc
   real(dp), allocatable, dimension(:) :: freqs, t, Inu, gamma_e, urad, &
         nu_obs, t_obs, gamma_bulk, R, tcool, Rb, volume, dotg_tmp
   real(dp), allocatable, dimension(:,:) :: dotg, n_e, jnut, jmbs, jssc, jeic, &
         ambs, anut, Qinj, tau_gg, pow_syn, Ddiff, Dopp
   logical :: full_rad_cool, bw_approx, radius_evol, pwl_over_trpzd_integ
   ! type(blast_wave), allocatable, dimension(:) :: bw

   !--->  Params
   call read_params(params_file)
   eps_e = par_eps_e
   eps_B = par_eps_B
   eps_g2 = par_eps_acc
   d_lum = par_d_lum
   z = par_z
   gmin = par_gmin
   gmax = par_gmax
   pind = par_pind
   numin = par_numin
   numax = par_numax
   numg = par_NG
   R0 = par_R0
   E0 = par_E0
   write(ifile, fmt) 0
   mezcal_file = "rads/iso_shock."//trim(ifile)//".out"!"gaus_shock."//trim(ifile)//".out"
   numh = count_lines(mezcal_file)
   dt = 0.1d0

   !--->  Setup
   write(*, "('--> Simulation setup')")

   allocate(nu_com(numf), gamma_e(numg), theta(numh), t_lab(0:numt), Inu(numf))
   allocate(r(numh, 0:numt), gamma_bulk(numh, 0:numt), Dopp(numh, 0:numt), &
         t_com(numh, 0:numt), volume(numh, 0:numt), dotg(numg, numh), &
         rho(numh, 0:numt), beta_bulk(numh, 0:numt), Bfield(numh, 0:numt), &
         r_mu_obs(numh, 0:numt), v_mu_obs(numh, 0:numt), Qinj(numg, numh))
   allocate(n_e(numg, numh, 0:numt), jnut(numf, numh, numt), &
         jssc(numf, numh, numt), anut(numf, numh, numt))
   zeros_arr = zeros1D(numg, .true.)

   call K1_init
   call K2_init

   !!!TODO: Transform all these into arguments of the subroutine
!   with_ic = .true.
!   full_rad_cool = .true.
!   bw_approx = .true.
!   radius_evol = .false.
!   pwl_over_trpzd_integ = .false.
   b_const = dsqrt(32d0 * pi * eps_B) * cLight
   g1_const = eps_e * mass_p * (pind - 2d0) / ((pind - 1d0) * mass_e)
   g2_const = dsqrt(6d0 * pi * eCharge * eps_g2 / sigmaT)
   t_lab(0) = 0d0

   build_f: do i = 1, numf
      nu_com(i) = numin * (numax / numin)**(dble(i - 1) / dble(numf - 1))
   end do build_f

   build_g: do i = 1, numg
      gamma_e(i) = gmin * (gmax / gmin)**(dble(i - 1) / dble(numg - 1))
   end do build_g

   call bw_mezcal(mezcal_file, theta_los, d_lum, r(:, 0), theta, gamma_bulk(:, 0), rho(:, 0), r_mu_obs(:, 0), v_mu_obs(:, 0))
   beta_bulk(:, 0) = bofg(gamma_bulk(:, 0))
   t_lab(0) = r(1, 0) / beta_bulk(1, 0) * cLight
   ! TODO: setup time in the lab frame

   do l=1,numh

      !> Element of arc
      if ( l == 1 ) then
         dth = 0.5d0 * (theta(2) + theta(1))
      else if ( l == numh ) then
         dth = halfpi - 0.5d0 * (theta(numh) + theta(numh - 1))
      else
         dth = 0.5d0 * (theta(l + 1) - theta(l - 1))
      end if

      !> Crossectional area
      call bw_crossec_area(1, assume_blob, r(l, 0), gamma_bulk(l, 0), dth, rb, volume(l, 0), cs_area, Omega_j)

      !> Magnetic field
      Bfield(l, 0) = b_const * dsqrt(rho(l, 0) * (gamma_bulk(l, 0) - 1d0) * gamma_bulk(l, 0))
      uB = Bfield(l, 0)**2 / (8d0 * pi)

      !> Minimum and maximum Lorentz factors of the particles distribution
      g2 = g2_const / dsqrt(Bfield(l, 0))
      g1 = g1_const * (gamma_bulk(l, 0) - 1d0)

      !> Fraction of accreted kinetic energy injected into non-thermal electrons
      L_e = eps_e * cs_area * rho(l, 0) * cLight**3 * beta_bulk(l, 0) * gamma_bulk(l, 0) * (gamma_bulk(l, 0) - 1d0)
      Q0 = L_e * pwl_norm(volume(l, 0) * energy_e, pind - 1d0, g1, g2)
      ! Q0 = L_e / ( ( g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) &
      !       - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6) ) &
      !       * (4d0 * pi * r(l, 0)**3 / 3d0) * energy_e )
      ! L_e = 4d0 * gamma_bulk(l, 0) * rho(l, 0) * cLight**2
      ! Q0 = L_e * pwl_norm(volume(l, 0), pind, g1, g2)
      n_e(:, l, 0) = injection_pwl(t_com(l, 0), 1d200, gamma_e, g1, g2, pind, Q0)

   end do

   write(*,"('--> Calculating the emission')")
   if ( KNcool ) then
      write(*, "('--> Radiative cooling: Klein-Nishina')")
      output_file="KNcool-"//trim(output_file)
   else
      write(*, "('--> Radiative cooling: Thomson')")
      output_file="Thcool-"//trim(output_file)
   end if
   write(*,*) ''
   write(*, "('--> Calculating the emission')")
   write(*, *) ''
   write(*, *) "Wrting data in: ", trim(output_file)
   write(*, *) ''
   write(*, *) screan_head


   !>  Evolution
   time_loop: do i=1, numt

      write(ifile, fmt) i
      mezcal_file = "rads/iso_shock."//trim(ifile)//".out"!"gaus_shock."//trim(ifile)//".out"
      call bw_mezcal(mezcal_file, theta_los, d_lum, r(:, i), theta, gamma_bulk(:, i), rho(:, i), r_mu_obs(:, i), v_mu_obs(:, i))

      derroteros_loop: do l=1, numh

         dr = r(l, i) - r(l, i - 1)
         !> Element of arc
         if ( l == 1 ) then
            dth = 0.5d0 * (theta(2) + theta(1))
         else if ( l == numh ) then
            dth = halfpi - 0.5d0 * (theta(numh) + theta(numh - 1))
         else
            dth = 0.5d0 * (theta(l + 1) - theta(l - 1))
         end if

         !> Crossectional area
         call bw_crossec_area(1, assume_blob, r(l, i), gamma_bulk(l, i), dth, rb, volume(l, i), cs_area, Omega_j)

         Dopp(l, i) = Doppler(gamma_bulk(l, i), v_mu_obs(l, i))
         beta_bulk(l, i) = bofg(gamma_bulk(l, i))
         call rk2_arr(t_com(l, i - 1), 1d0 / (beta_bulk(l, i - 1:i) * gamma_bulk(l, i-1:i) * cLight), dr, t_com(l, i))
         dt = t_com(l, i) - t_com(l, i - 1)
         rb = r(l, i) / (12d0 * gamma_bulk(l, i))
         tlc = rb / cLight

         !> Solving the Fokker-Planck eq.
         call FP_FinDif_difu(dt, &
               &             gamma_e, &
               &             n_e(:, l, i - 1), &
               &             n_e(:, l, i), &
               &             dotg(:, i - 1), &
               &             zeros_arr, &
               &             Qinj(:, l), &
               &             1d200, &
               &             tlc)

         !> Magnetic field
         Bfield(l, i) = b_const * dsqrt(rho(l, i) * (gamma_bulk(l, i) - 1d0) * gamma_bulk(l, i))
         uB = Bfield(l, i)**2 / (8d0 * pi)

         !> Minimum and maximum Lorentz factors of the particles distribution
         g2 = g2_const / dsqrt(Bfield(l, i))
         g1 = g1_const * (gamma_bulk(l, i) - 1d0)

         !> Fraction of accreted kinetic energy injected into non-thermal electrons
         L_e = eps_e * 4d0 * pi * r(l, i)**2 * rho(l, i) * cLight**3 * beta_bulk(l, i) * gamma_bulk(l, i) * (gamma_bulk(l, i) - 1d0)
         Q0 = L_e * pwl_norm(4d0 * pi * r(l, i)**3 * energy_e / (12d0 * gamma_bulk(l, i)), pind - 1d0, g1, g2)
         ! Q0 = L_e / ( ( g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) &
         !       - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6) ) &
         !       * (4d0 * pi * r(l, i)**3 / 3d0) * energy_e )
         ! L_e = 4d0 * gamma_bulk(l, i) * n_ext * eps_e
         ! Q0 = L_e * pwl_norm(4d0 * pi * r(l, i)**3 / 3d0, pind, g1, g2)
         Qinj(:, l) = injection_pwl(t_com(l, i), 1d200, gamma_e, g1, g2, pind, Q0)
         n_e(:, l, i) = Qinj(:, l)
         !> Synchrotron cooling
         dotg(:, l) = b_const * Bfield(l, i)**2 * gamma_e**2 &
         !> Adiabatic cooling
#if( ADIABATIC_COOLING_TYPE == ADIABATIC_MSB00 )
               + cLight * beta_bulk(l, i) * gamma_bulk(l, i) * gamma_e / r(l, i)
#elif( ADIABATIC_COOLING_TYPE == ADIABATIC_VOL_EVOL )
               + pofg(gamma_e) * dlog(volume(l, i) / volume(l, i - 1)) / (3d0 * dt)
#else
               + 0d0
#endif

         !> Radiation
         !$OMP PARALLEL DO SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j)
         do j = 1, numf
            call syn_emissivity(jsyn(l, j, i), nu_com(j), gamma_e, n_e(:, l, i), Bfield(l, i))
            call syn_absorption(anut(l, j, i), nu_com(j), gamma_e, n_e(:, l, i), Bfield(l, i))
         end do
         !$OMP END PARALLEL DO

         if ( with_ic ) then
            call RadTrans(Inu, rb, jsyn(l, :, i), anut(l, :, i))
            !$OMP PARALLEL DO SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j)
            do j = 1, numf
               call IC_iso_powlaw(jssc(l, j, i), nu_com(j), nu_com, Inu, n_e(:, l, i), gamma_e)
               ! call IC_iso_monochrom(jeic(j, i), freqs(j), uext, nu_ext, n_e(:, l, i), gamma_e)
               jnut(l, j, i) = jsyn(l, j, i) + jssc(l, j, i)! + jeic(j, i)
            end do
            !$OMP END PARALLEL DO
         else
            jnut(l, :, i) = jsyn(l, :, i)
         end if

      end do derroteros_loop

      if ( mod(i, nmod) == 0 .or. i == 1 ) &
            write(*, on_screen) i, t_lab(i), r(i, 0)

   end do time_loop

   !  ####    ##   #    # # #    #  ####
   ! #       #  #  #    # # ##   # #    #
   !  ####  #    # #    # # # #  # #
   !      # ###### #    # # #  # # #  ###
   ! #    # #    #  #  #  # #   ## #    #
   !  ####  #    #   ##   # #    #  ####
   write(*, "('--> Saving')")
#ifdef HDF5
   ! ------  Opening output file  ------
   call h5open_f(herror)
   call h5io_createf(output_file, file_id, herror)
   ! ------  Saving initial parameters  ------
   call h5io_createg(file_id, "Parameters", group_id, herror)
   call h5io_wint0 (group_id, 'NT',        numt, herror)
   call h5io_wint0 (group_id, 'NF',        numf, herror)
   call h5io_wint0 (group_id, 'NG',        numg, herror)
   call h5io_wint0 (group_id, 'NH',        numh, herror)
   call h5io_wdble0(group_id, 'd_lum',     d_lum, herror)
   call h5io_wdble0(group_id, 'redshift',  z, herror)
   call h5io_wdble0(group_id, 'gamma_min', gmin, herror)
   call h5io_wdble0(group_id, 'gamma_max', gmax, herror)
   call h5io_wdble0(group_id, 'gamma_1',   g1, herror)
   call h5io_wdble0(group_id, 'gamma_2',   g2, herror)
   call h5io_wdble0(group_id, 'pwl-index', pind, herror)
   call h5io_wdble0(group_id, 'epsilon_e', eps_e, herror)
   call h5io_wdble0(group_id, 'epsilon_B', eps_B, herror)
   call h5io_wdble0(group_id, 'nu_min',    numin, herror)
   call h5io_wdble0(group_id, 'nu_max',    numax, herror)
   call h5io_wdble0(group_id, 'E0',        E0, herror)
   call h5io_closeg(group_id, herror)
   ! ------  Saving data  ------
   call h5io_wdble1(file_id, 't_lab',      t_lab, herror)
   call h5io_wdble1(file_id, 'theta',      theta, herror)
   call h5io_wdble1(file_id, 'nu_com',     nu_com, herror)
   call h5io_wdble1(file_id, 'gamma_e',    gamma_e, herror)
   call h5io_wdble2(file_id, 't_com',      t_com, herror)
   call h5io_wdble2(file_id, 'radius',     r, herror)
   call h5io_wdble2(file_id, 'volume',     volume, herror)
   call h5io_wdble2(file_id, 'Gamma_bulk', gamma_bulk, herror)
   call h5io_wdble2(file_id, 'Doppler',    Dopp, herror)
   call h5io_wdble2(file_id, 'r_mu_obs',   r_mu_obs, herror)
   call h5io_wdble2(file_id, 'v_mu_obs',   v_mu_obs, herror)
   call h5io_wdble2(file_id, 'density',    rho, herror)
   call h5io_wdble3(file_id, 'jnut',       jnut, herror)
   call h5io_wdble3(file_id, 'jsyn',       jsyn, herror)
   call h5io_wdble3(file_id, 'jssc',       jssc, herror)
   ! call h5io_wdble3(file_id, 'jeic',       jeic, herror)
   call h5io_wdble3(file_id, 'anut',       anut, herror)
   ! call h5io_wdble3(file_id, 'asyn',       ambs, herror)
   call h5io_wdble3(file_id, 'n_e',        n_e, herror)
   ! ------  Closing output file  ------
   call h5io_closef(file_id, herror)
   call h5close_f(herror)
#endif
   write(*, "('==========  FINISHED  ==========')")
   write(*,*) ''

end subroutine mezcal
#endif
