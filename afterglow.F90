subroutine afterglow(params_file, output_file, with_abs, with_cool, with_ic)
   use data_types
   use constants
   use params
   use misc
   use pwl_integ
   use hdf5
   use h5_inout
   use SRtoolkit
   use anaFormulae
   use dist_evol
   use radiation
   use pairs
   use Aglow_models
   use K1
   use K2
   implicit none

   character(len=*), intent(in) :: output_file, params_file
   logical, intent(in) :: with_cool, with_ic, with_abs
   integer, parameter :: nmod = 50
   character(len=*), parameter :: screan_head = &
      '| Iteration | Obser. time |   BW radius |  gamma_bulk |       N_tot |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------', &
      on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
   integer(HID_T) :: file_id, group_id
   integer :: i, j, k, numbins, numdf, numdt, time_grid, herror, beam_kind
   real(dp) :: uB, uext, L_j, gmin, gmax, numin, numax, pind, B, R0, &
      tinj, g1, g2, tstep, Q0, tmax, d_lum, z, n_ext, urad, &
      theta_obs, mu_obs, nu_ext, tesc_e, uext0, volume, eps_e, tlc, &
      eps_B, E0, gamma_bulk0, L_e, nu_ext0, tmin, td, Rd, dr, &
      b_const, beta_bulk, eps_g2, theta_j0, Rb, cs_area, volume_p, &
      g2_const, g1_const, energy_e, energy_p, Rd2, Omega_j, dt
   real(dp), allocatable, dimension(:) :: freqs, t, Ntot, Inu, gg, &
      nu_obs, t_obs, gamma_bulk, R, D, tcool, gc, nu_c, Ddiff
   real(dp), allocatable, dimension(:, :) :: dotg, n_e, jnut, jmbs, jssc, jeic, &
      ambs, anut, Qinj, Fmbs, Feic, Fssc, Fnut, tau_gg
   logical :: blob, adiab_coef

   !!!!!!!! TODO:
   !!  [ ] Cases for each model to compare:
   !!       - SPN98
   !!       - PM09
   !!       - PVP14
   !!       - Hao19
   !!  [ ] Winds
   !!  [ ] Pair production
   !!  [ ] Solve photons kinetic equation or Integrate emissivities
   !!  [ ] Time dependent adiabatic cooling
   !!  [ ] Dilution term
   !!  [X] Save cooling break lorentz factor
   !!  [X] Geometry of the emission region: blob (size of the jet), slab (DM09), a region from eq. (30) in PVP14
   !!  [X] Jet vs isotropic expansion
   !!  [X] Swept up material

   !  #####    ##   #####    ##   #    #  ####
   !  #    #  #  #  #    #  #  #  ##  ## #
   !  #    # #    # #    # #    # # ## #  ####
   !  #####  ###### #####  ###### #    #      #
   !  #      #    # #   #  #    # #    # #    #
   !  #      #    # #    # #    # #    #  ####
   call read_params(params_file)
   eps_e = par_eps_e
   eps_B = par_eps_B
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
   numbins = par_numbins
   numdt = par_numdt
   numdf = par_numdf
   n_ext = par_n_ext
   g1 = par_g1
   g2 = par_g2
   time_grid = par_time_grid
   R0 = par_R0
   tmax = par_tmax
   E0 = par_E0

   !-----> Misc.
   energy_e = mass_e * cLight**2
   energy_p = mass_p * cLight**2

   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   write(*, "('--> Simulation setup')")

   allocate(t(0:numdt), freqs(numdf), Ntot(numdt), Inu(numdf), gg(numbins), &
      nu_obs(numdf), t_obs(0:numdt), R(0:numdt), tcool(numbins), &
      gamma_bulk(0:numdt), D(0:numdt), gc(0:numdt), nu_c(0:numdt), &
      Ddiff(numbins))
   allocate(n_e(numbins, 0:numdt), dotg(numbins, 0:numdt), &
      ambs(numdf, numdt), jmbs(numdf, numdt), jnut(numdf, numdt), &
      jssc(numdf, numdt), anut(numdf, numdt), jeic(numdf, numdt), &
      Qinj(numbins, 0:numdt), Fnut(numdf, numdt), &
      Fmbs(numdf, numdt), Fssc(numdf, numdt), Feic(numdf, numdt), &
      tau_gg(numdf, numdt))

   call K1_init
   call K2_init

   !-----> About the observer
   theta_obs = par_theta_obs * pi / 180d0
   mu_obs = dcos(theta_obs)
   t_obs(0) = 0d0

   !-----> Initializing blast wave
   theta_j0 = 0.2d0
   beam_kind = -1
   blob = .false.
   adiab_coef = .false.
   beta_bulk = bofg(gamma_bulk0)
   call bw_crossec_area(gamma_bulk0, R0, gamma_bulk0, theta_j0, beam_kind, blob, Rb, volume, cs_area, Omega_j)

   !-----> True outflow energy and decceleration radius
   !!!!!NOTE: multiply by 2d0 if a doulbe-sided jet is considered
   E0 = E0 * Omega_j / (4d0 * pi)
   Rd = deceleration_radius(E0, gamma_bulk0, n_ext) ! eq. (1) in RM92
   td = (1d0 + z) * Rd / (beta_bulk * gamma_bulk0**2 * cLight) ! eq. (11.4) of DM09
   Rd2 = Rd * 2d0**(-2d0 / 3d0) ! see R_B in PVP14, p. 3

   !-----> Locating the emission region
   R(0) = par_R
   gamma_bulk(0) = adiab_blast_wave(R(0), R0, gamma_bulk0, E0, n_ext)
   beta_bulk = bofg(gamma_bulk(0))
   D(0) = Doppler(gamma_bulk(0), mu_obs)
   call bw_crossec_area(gamma_bulk0, R(0), gamma_bulk(0), theta_j0, beam_kind, blob, Rb, volume, cs_area, Omega_j)

   !-----> External medioum
   !!!!!!COMBAK: implement radius dependent external medium
   !!!!!!WARNING: a radius dependent external medium may change the blast wave solution
   !-----> Magnetic field
   b_const = dsqrt(32d0 * pi * eps_B * mass_p) * cLight
   B = b_const * dsqrt(n_ext) * gamma_bulk(0)
   uB = B**2 / (8d0 * pi)

   !-----> Transforming monochromatic ext. rad. field to comoving frame
   uext = uext0 * gamma_bulk(0)**2 * (1d0 + beta_bulk**2 / 3d0) ! eq. (5.25) in DM09
   nu_ext = nu_ext0 * gamma_bulk(0)

   !-----> Characteristic Lorentz factor (comoving) and frequency (observer)
   gc(0) = 6d0 * gamma_bulk(0) * energy_e / (5d0 * sigmaT * R(0) * uB)
   nu_c(0) = nu_obs_f(nuconst * B * gc(0)**2, z, D(0))

   !-----> Minimum and maximum Lorentz factors of the particles distribution
   eps_g2 = 0.35
   g2_const = dsqrt(6d0 * pi * eCharge * eps_g2 / sigmaT)
   g1_const = eps_e * mass_p * (pind - 2d0) / ((pind - 1d0) * mass_e)
   g2 = g2_const / dsqrt(B)
   g1 = g1_const * (gamma_bulk(0) - 1d0)

   !-----> Fraction of accreted kinetic energy injected into non-thermal electrons
   L_e = eps_e * cs_area * n_ext * energy_p * cLight * beta_bulk * gamma_bulk(0) * (gamma_bulk(0) - 1d0)
   Q0 = L_e * pwl_norm(volume * energy_e, pind - 1d0, g1, g2)

   write(*, "('mu_obs  =', ES15.7)") mu_obs
   write(*, "('Doppler =', ES15.7)") D(0)
   write(*, "('gamma_1 =', ES15.7)") g1
   write(*, "('gamma_2 =', ES15.7)") g2
   write(*, "('L_e     =', ES15.7)") L_e
   write(*, "('Q_inj   =', ES15.7)") Q0
   write(*, "('B0      =', ES15.7)") B
   write(*, "('u_B     =', ES15.7)") uB
   write(*, "('u_ext   =', ES15.7)") uext
   write(*, "('nu_ext  =', ES15.7)") nu_ext
   write(*, "('Gamma_0 =', ES15.7)") gamma_bulk(0)
   write(*, "('R_0     =', ES15.7)") R(0)
   write(*, "('R_dec   =', ES15.7)") Rd
   write(*, "('t_dec   =', ES15.7)") td
   write(*, "('t_dyn   =', ES15.7)") tlc
   write(*, "('tinj    =', ES15.7)") tinj
   write(*, "('tesc    =', ES15.7)") tesc_e

   build_f: do j = 1, numdf
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numdf - 1)) )
   end do build_f

   build_g: do k = 1, numbins
      gg(k) = (gmin - 1d0) * ((gmax - 1d0) / (gmin - 1d0))**(dble(k - 1) / dble(numbins - 1)) + 1d0
   end do build_g

   dotg(:, 0) = 4d0 * sigmaT * uB * pofg(gg)**2 / (3d0 * mass_e * cLight)
   Inu = 0d0
   Ddiff = 1d-200!0.5d0 * gg**2 / tacc !
   n_e(:, 0) = injection_pwl(0d0, 1d200, gg, g1, g2, pind, Q0)
   if ( with_cool ) then
      urad = uext
   else
      urad = 0d0
   end if

   write(*, "('---> Calculating the emission')")
   write(*, *) ''
   write(*, "(' Using tstep = ', F7.3)") tstep
   write(*, "(' Initial observer time = ', F7.3)") tmin
   write(*, *) "Wrting data in: ", trim(output_file)
   write(*, *) ''
   write(*, *) screan_head

   ! ###### #    #  ####  #      #    # ##### #  ####  #    #
   ! #      #    # #    # #      #    #   #   # #    # ##   #
   ! #####  #    # #    # #      #    #   #   # #    # # #  #
   ! #      #    # #    # #      #    #   #   # #    # #  # #
   ! #       #  #  #    # #      #    #   #   # #    # #   ##
   ! ######   ##    ####  ######  ####    #   #  ####  #    #
   time_loop: do i = 1, numdt
      t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numdt - 1)) )
      dt = t(i) - t(i - 1)

      !-----> Locating the emission region
      dr = dt * beta_bulk * gamma_bulk(i - 1) * cLight
      R(i) = R(i - 1) + dr
      gamma_bulk(i) = adiab_blast_wave(R(i), R0, gamma_bulk0, E0, n_ext)
      D(i) = Doppler(gamma_bulk(i), mu_obs)
      beta_bulk = bofg(gamma_bulk(i))
      t_obs(i) = t_obs(i - 1) + 0.5d0 * dt * (1d0 / D(i) + 1d0 / D(i - 1))
      call bw_crossec_area(gamma_bulk0, R(i), gamma_bulk(i), theta_j0, beam_kind, blob, Rb, volume, cs_area, Omega_j)

      !-----> External medioum
      !!!!!COMBAK
      !-----> Magnetic field
      B = b_const * dsqrt(n_ext) * gamma_bulk(i)
      uB = B**2 / (8d0 * pi)

      !-----> Transforming monochromatic ext. rad. field to comoving frame
      uext = uext0 * gamma_bulk(i)**2 * (1d0 + beta_bulk**2 / 3d0) ! eq. (5.25) in DM09
      nu_ext = nu_ext0 * gamma_bulk(i)

      !-----> Minimum and maximum Lorentz factors of the particles distribution
      g2 = g2_const / dsqrt(B)
      g1 = g1_const * (gamma_bulk(i) - 1d0)

      !-----> Time-scales
      tlc = Rb / cLight
      tesc_e = 1d200!1.01d0 * tlc!
      tinj = 1d200
      ! tcool = 1d0 / (nu0(:, i) * pofg(gg) + Aadi)
      ! tadiab = 1d0 / (Aadi * gg(ig) * bofg(gg(ig)))

      !-----> Fraction of accreted kinetic energy injected into non-thermal electrons
      L_e = eps_e * cs_area * n_ext * energy_p * cLight * beta_bulk * gamma_bulk(0) * (gamma_bulk(0) - 1d0)
      Q0 = L_e * pwl_norm(volume * energy_e, pind - 1d0, g1, g2)
      Qinj(:, i) = injection_pwl(t(i), tinj, gg, g1, g2, pind, Q0)

      !  ###### ###### #####
      !  #      #      #    #
      !  #####  #####  #    #
      !  #      #      #    #
      !  #      #      #    #
      !  ###### ###### #####
      ! Ddiff(:, i) = 1d-200
      call FP_FinDif_difu(dt, &
            &             pofg(gg), &
            &             n_e(:, i - 1), &
            &             n_e(:, i), &
            &             dotg(:, i - 1), &
            &             Ddiff, &
            &             Qinj(:, i), &
            &             tesc_e, Rb)

      !!!!!TEMP: this is just to test deleting low energy particles
      ! do k = 1, numbins
      !    if ( g(k) < dmin1(g1, gc(i - 1)) ) n(k, i) = 0d0
      ! end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !  #####    ##   #####  #   ##   ##### #  ####  #    #
      !  #    #  #  #  #    # #  #  #    #   # #    # ##   #
      !  #    # #    # #    # # #    #   #   # #    # # #  #
      !  #####  ###### #    # # ######   #   # #    # #  # #
      !  #   #  #    # #    # # #    #   #   # #    # #   ##
      !  #    # #    # #####  # #    #   #   #  ####  #    #
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numdf
         freqs(j) = nu_com_f(nu_obs(j), z, D(i))
         call mbs_emissivity(jmbs(j, i), freqs(j), gg, n_e(:, i), B)
         ambs(j, i) = 0d0
         if ( with_abs ) call mbs_absorption(ambs(j, i), freqs(j), gg, n_e(:, i), B)
      end do
      !$OMP END PARALLEL DO

      if ( blob ) then
         call RadTrans_blob(Inu, Rb, jmbs(:, i), ambs(:, i))
      else
         call RadTrans(Inu, Rb, jmbs(:, i), ambs(:, i))
      end if

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numdf
         if ( with_ic ) then
            call IC_iso_powlaw(jssc(j, i), freqs(j), freqs, Inu, n_e(:, i - 1) / volume, gg)
            call IC_iso_monochrom(jeic(j, i), freqs(j), uext, nu_ext, n_e(:, i - 1) / volume, gg)
            ! call IC_emis_full(freqs(j), freqs, gg, n_e(:, i - 1) / volume, Inu, jssc(j, i))
            ! call IC_emis_full(freqs(j), nu_ext, gg, n_e(:, i - 1) / volume, uext * cLight / (4d0 * pi), jeic(j, i))
         else
            jssc(j, i) = 0d0
            jeic(j, i) = 0d0
         end if

         !!!!!COMBAK: pairs optical depth
         ! if ( hPlanck * freqs(j) > 5d11 * (0.01d0 / 0.02d0) * (100d0 / gamma_bulk(i)) ) then
         !    tau_gg(j, i) = 0.16d0 * (R(i) / 1d16) * (100d0 / gamma_bulk(i)) * (uext / 1d-7) * (1d11 / (hPlanck * freqs(j))) * (0.01d0 / 0.02d0)**2
         ! else
         ! tau_gg(j, i) = 0d0
         ! end if

         anut(j, i) = ambs(j, i)! + tau_gg(j, i) / (2d0 * Rb)

         if ( blob ) then
            Fmbs(j, i) = D(i)**4 * volume * freqs(j) * jmbs(j, i) * opt_depth_blob(anut(j, i), Rb) / (4d0 * pi * d_lum**2)
            Fssc(j, i) = D(i)**4 * volume * freqs(j) * jssc(j, i) * opt_depth_blob(anut(j, i), Rb) / (4d0 * pi * d_lum**2)
            Feic(j, i) = D(i)**4 * volume * freqs(j) * jeic(j, i) * opt_depth_blob(anut(j, i), Rb) / (4d0 * pi * d_lum**2)
         else

            ! Fmbs(j, i) = D(i)**4 * volume * freqs(j) * jmbs(j, i) / (4d0 * pi * d_lum**2)
            ! Fssc(j, i) = D(i)**4 * volume * freqs(j) * jssc(j, i) / (4d0 * pi * d_lum**2)
            ! Feic(j, i) = D(i)**4 * volume * freqs(j) * jeic(j, i) / (4d0 * pi * d_lum**2)
            Fmbs(j, i) = D(i)**4 * volume * freqs(j) * jmbs(j, i) * opt_depth_slab(anut(j, i), Rb) / d_lum**2
            Fssc(j, i) = D(i)**4 * volume * freqs(j) * jssc(j, i) * opt_depth_slab(anut(j, i), Rb) / d_lum**2
            Feic(j, i) = D(i)**4 * volume * freqs(j) * jeic(j, i) * opt_depth_slab(anut(j, i), Rb) /d_lum**2

         end if

         jnut(j, i) = jmbs(j, i) + jssc(j, i) + jeic(j, i)
         Fnut(j, i) = Fmbs(j, i) + Fssc(j, i) + Feic(j, i)

      end do
      !$OMP END PARALLEL DO

      !   ####   ####   ####  #      # #    #  ####
      !  #    # #    # #    # #      # ##   # #    #
      !  #      #    # #    # #      # # #  # #
      !  #      #    # #    # #      # #  # # #  ###
      !  #    # #    # #    # #      # #   ## #    #
      !   ####   ####   ####  ###### # #    #  ####
      dotg(:, i) = 4d0 * sigmaT * uB * pofg(gg)**2 / (3d0 * mass_e * cLight)
      !!!!!COMBAK: setup radiative and adiabatic cooling
      !-----> Radiative cooling
      !-----> Adiabatic cooling
      ! if ( i == 1 .or. adiab_coef_an ) then
      !    if ( R(i) < Rd2 ) then
      !       gindex = 1d0
      !    else
      !       gindex = 2.5d0
      !    end if
      !    nu0(:, i) = (4d0 * sigmaT * uB * pofg(gg)**2 / (3d0 * mass_e * cLight)) + (pofg(gg) * bofg(gg) * ((2d0 / gindex) + 1d0) / (3d0 * t(i)))
      ! else
      !    nu0(:, i) = (4d0 * sigmaT * uB * pofg(g)**2 / (3d0 * mass_e * cLight)) + pofg(gg) * (dlog(volume) - dlog(volume_p)) / (3d0 * dt)
      ! end if
      ! Aadi = (8d0 / 5d0) * cLight * beta_bulk * gamma_bulk(i) / R(i)
      ! nu0(:, i) = 4d0 * sigmaT * (uB + urad) / (3d0 * mass_e * cLight)

      !-----> Characteristic Lorentz factor (comoving) and frequency (observer)
      gc(i) = 6d0 * gamma_bulk(i) * energy_e / (5d0 * sigmaT * R(i) * uB)
      nu_c(i) = nu_obs_f(nuconst * B * gc(i)**2, z, D(i))

      !-----> Saving volume of this time-step
      volume_p = volume

      !   ####  #    #     ####   ####  #####  ###### ###### #    #
      !  #    # ##   #    #      #    # #    # #      #      ##   #
      !  #    # # #  #     ####  #      #    # #####  #####  # #  #
      !  #    # #  # #         # #      #####  #      #      #  # #
      !  #    # #   ##    #    # #    # #   #  #      #      #   ##
      !   ####  #    #     ####   ####  #    # ###### ###### #    #
      Ntot(i) = 0.5d0 * sum((n_e(:numbins - 1, i - 1) + n_e(2:, i - 1)) * (gg(2:) - gg(:numbins - 1)))

      if ( mod(i, nmod) == 0 .or. i == 1 ) &
         write(*, on_screen) i, t_obs(i - 1), R(i - 1), gamma_bulk(i - 1), Ntot(i)

   end do time_loop


   !  ####    ##   #    # # #    #  ####
   ! #       #  #  #    # # ##   # #    #
   !  ####  #    # #    # # # #  # #
   !      # ###### #    # # #  # # #  ###
   ! #    # #    #  #  #  # #   ## #    #
   !  ####  #    #   ##   # #    #  ####
   write(*, "('--> Saving')")

   ! ------  Opening output file  ------
   call h5open_f(herror)
   call h5io_createf(output_file, file_id, herror)

   ! ------  Saving initial parameters  ------
   call h5io_createg(file_id, "Parameters", group_id, herror)

   call h5io_wint0(group_id, 'numdt', numdt, herror)
   call h5io_wint0(group_id, 'numdf', numdf, herror)
   call h5io_wint0(group_id, 'numbins', numbins, herror)
   call h5io_wint0(group_id, 'time-grid', time_grid, herror)

   call h5io_wdble0(group_id, 't_max', tmax, herror)
   call h5io_wdble0(group_id, 't_min', tmin, herror)
   call h5io_wdble0(group_id, 'tstep', tstep, herror)
   call h5io_wdble0(group_id, 'd_lum', d_lum, herror)
   call h5io_wdble0(group_id, 'redshift', z, herror)
   call h5io_wdble0(group_id, 'gamma_bulk0', gamma_bulk0, herror)
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
   call h5io_wdble0(group_id, 'R0', R0, herror)
   call h5io_wdble0(group_id, 'n_ext', n_ext, herror)

   call h5io_closeg(group_id, herror)

   ! ------  Saving data  ------
   call h5io_wdble0(file_id, 'Rd', Rd, herror)
   call h5io_wdble0(file_id, 'td', td, herror)

   call h5io_wdble1(file_id, 'time', t(1:), herror)
   call h5io_wdble1(file_id, 't_obs', t_obs(1:), herror)
   call h5io_wdble1(file_id, 'Ntot', Ntot, herror)
   call h5io_wdble1(file_id, 'Rbw', R(1:), herror)
   call h5io_wdble1(file_id, 'Gamma_bulk', gamma_bulk(1:), herror)
   call h5io_wdble1(file_id, 'nu', freqs, herror)
   call h5io_wdble1(file_id, 'nu_obs', nu_obs, herror)
   call h5io_wdble1(file_id, 'gamma', gg, herror)
   call h5io_wdble1(file_id, 'gamma_c', gc(1:), herror)
   call h5io_wdble1(file_id, 'nu_c', nu_c(1:), herror)

   call h5io_wdble2(file_id, 'jnut', jnut, herror)
   call h5io_wdble2(file_id, 'jmbs', jmbs, herror)
   call h5io_wdble2(file_id, 'jssc', jssc, herror)
   call h5io_wdble2(file_id, 'jeic', jeic, herror)
   call h5io_wdble2(file_id, 'anut', anut, herror)
   call h5io_wdble2(file_id, 'ambs', ambs, herror)
   call h5io_wdble2(file_id, 'Fnut', Fnut, herror)
   call h5io_wdble2(file_id, 'Fmbs', Fmbs, herror)
   call h5io_wdble2(file_id, 'Fssc', Fssc, herror)
   call h5io_wdble2(file_id, 'Feic', Feic, herror)
   call h5io_wdble2(file_id, 'Qinj', Qinj, herror)
   call h5io_wdble2(file_id, 'n_e', n_e(:, 1:), herror)
   call h5io_wdble2(file_id, 'cool-coef', dotg(:, 1:), herror)

   ! ------  Closing output file  ------
   call h5io_closef(file_id, herror)
   call h5close_f(herror)

   write(*, "('==========  FINISHED  ==========')")
   write(*,*) ''

end subroutine afterglow
