subroutine afterglow(params_file, output_file, with_cool)
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
   use models
   use K1
   use K2
   implicit none

   character(len=*), intent(in) :: output_file, params_file
   logical, intent(in) :: with_cool
   integer, parameter :: nmod = 10
   character(len=*), parameter :: screan_head = &
      '| Iteration |        Time |   Time step |    nu_0(g2) |       N_tot |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------', &
      on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
   integer(HID_T) :: file_id, group_id
   integer :: i, j, k, numbins, numdf, numdt, time_grid, herror
   real(dp) :: uB, uext, urad, ue, R, L_j, gmin, gmax, numin, numax, qind, B, &
      tacc, g1, g2, tstep, zetae, Q0, theta_e, tmax, d_lum, z, D, ne, n_ext, &
      gamma_bulk, theta_obs, R0, mu_obs, nu_ext, tesc, x, y, dt_var, &
      volume, sigma, beta_bulk, eps_e, eps_B, E0, gamma_bulk0, Iind, L_e, tau, &
      eta_delta, bw_radius, nu_ext0
   real(dp), allocatable, dimension(:) :: freqs, t, Ntot, Inu, gg, sen_lum, &
      dt, nu_obs, t_obs
   real(dp), allocatable, dimension(:, :) :: nu0, nn, jnut, jmbs, jssc, jeic, &
      ambs, anut, Qinj, Ddif, Fmbs, Feic, Fssc, Fnut

   call read_params(params_file)
   eps_e = par_eps_e
   eps_B = par_eps_B
   gamma_bulk0 = par_gamma_bulk
   d_lum = par_d_lum
   ! theta_obs = par_theta_obs
   z = par_z
   bw_radius = par_R0
   tstep = par_tstep
   tmax = par_tmax
   dt_var = par_tvar
   g1 = par_g1
   g2 = par_g2
   gmin = par_gmin
   gmax = par_gmax
   qind = par_qind
   numin = par_numin
   numax = par_numax
   nu_ext0 = par_nu_ext
   numbins = par_numbins
   numdt = par_numdt
   numdf = par_numdf

   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   write(*, "('--> Simulation setup')")

   allocate(t(0:numdt), freqs(numdf), Ntot(numdt), Inu(numdf), sen_lum(numdt), &
      dt(numdt), nu_obs(numdf), t_obs(numdt))
   allocate(nn(0:numdt, numbins), nu0(0:numdt, numbins), gg(numbins), &
      ambs(numdf, numdt), jmbs(numdf, numdt), jnut(numdf, numdt), &
      jssc(numdf, numdt), anut(numdf, numdt), jeic(numdf, numdt), &
      Qinj(numdt, numbins), Ddif(numdt, numbins), Fnut(numdf, numdt), &
      Fmbs(numdf, numdt), Fssc(numdf, numdt), Feic(numdf, numdt))

   call K1_init
   call K2_init

   E0 = 1d55
   n_ext = 0.1d0
   theta_obs = par_theta_obs * pi / 180d0
   mu_obs = dcos(theta_obs)
   gamma_bulk = gamma_bulk0
   eta_delta = 1d0
   ! R = R0

   tacc = tmax !0.5d0 * R / cLight

   write(*, "('theta_obs =', ES15.7)") theta_obs
   write(*, "('nu_ext0   =', ES15.7)") nu_ext0


   build_f: do j=1,numdf
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numdf - 1)) )
   end do build_f

   build_g: do k = 1, numbins
      gg(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numbins - 1))
   end do build_g

   t(0) = 0d0
   nn = 0d0

   write(*, "('---> Calculating the emission')")
   write(*, *) ''
   write(*, "(' Using tstep = ', F5.3)") tstep
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

      t_obs(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numdt - 1)) )

      call shock_afterglow(t_obs(i) / (1d0 + z), gamma_bulk0, E0, eps_e, eps_B, n_ext, B, gamma_bulk, bw_radius, ne, ue, .true.)
      R = eta_delta * gamma_bulk * cLight * dt_var / (1d0 + z)
      volume = 4d0 * pi * R**3 / 3d0 !(12d0 * gamma_bulk)
      beta_bulk = bofg(gamma_bulk)
      D = Doppler(gamma_bulk, mu_obs)
      t(i) = D * t_obs(i) / (1d0 + z) !t_com_f(t_obs(i), z, gamma_bulk, bw_radius, mu_obs)!
      dt(i) = t(i) - t(i - 1)
      tesc = 1d200!1.5d0 * R / cLight

      L_e = eps_e * 4d0 * pi * R**2 * n_ext * mass_p * cLight**3 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0)
      Q0 = L_e * pwl_norm(mass_e * cLight**2 * volume, qind - 1d0, g1, g2)

      B = dsqrt(32d0 * pi * mass_p * eps_B * n_ext) * gamma_bulk * cLight
      uB = 0.125d0 * B**2 / pi
      y = 125d0
      ! u_ph = 7.5d-7
      uext = y * eps_B * 4d0 * n_ext * mass_p * cLight**2 * gamma_bulk**2 * (1d0 + beta_bulk**2 / 3d0)
      nu_ext = nu_ext0 * gamma_bulk

      if ( i == 1 )  nu0(i - 1, :) = 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)

      Qinj(i, :) = injection(t(i), tacc, gg, g1, g2, qind, theta_e, 0d0, Q0)
      Ddif(i, :) = 1d-200 !0.5d0 * gg**2 / tacc !
      call FP_FinDif_difu(dt(i), gg, nn(i - 1, :), nn(i, :), nu0(i - 1, :), Ddif(i, :), Qinj(i, :), tesc)
      ! call FP_FinDif_cool(dt(i), gg, nn(i - 1, :), nn(i, :), nu0(i - 1, :), Qinj(i, :), tesc)

      do j=1,numdf
         freqs(j) = nu_com_f(nu_obs(j), z, gamma_bulk, mu_obs)
      end do
      call mbs_emissivity(jmbs(:, i), freqs, gg, nn(i, :), B)
      call mbs_absorption(ambs(:, i), freqs, gg, nn(i, :), B)
      call RadTrans_blob(Inu, R, jmbs(:, i), ambs(:, i))
      call SSC_pwlEED(jssc(:, i), freqs, Inu, nn(i, :), gg)
      call EIC_pwlEED(jeic(:, i), freqs, uext, nu_ext, nn(i, :), gg)

      jnut(:, i) = jmbs(:, i) + jssc(:, i) + jeic(:, i)
      anut(:, i) = ambs(:, i)

      if ( with_cool ) then
         do j = 2, numdf
            if ( Inu(j - 1) > 1d-200 .and. Inu(j) > 1d-200 ) then
               Iind = -dlog(Inu(j) / Inu(j - 1)) / dlog(freqs(j) / freqs(j - 1))
               if ( Iind > 8d0 ) Iind = 8d0
               if ( Iind < -8d0 ) Iind = -8d0
               urad = urad + Inu(j - 1) * freqs(j - 1) * Pinteg(freqs(j) / freqs(j - 1), Iind, 1d-9)
            end if
         end do
         urad = 4d0 * pi * urad / cLight
         urad = urad + IC_cool(gg, freqs, R * (jssc(:, i) + jeic(:, i)))
      else
         urad = 0d0
      end if

      nu0(i, :) = 4d0 * sigmaT * (uB + urad) / (3d0 * mass_e * cLight)

      do j = 1, numdf
         tau = dmax1(1d-200, 2d0 * R * ambs(j, i))
         Fmbs(j, i) = 3d0 * D**4 * volume * freqs(j) * opt_depth_blob(tau) * jmbs(j, i) / (4d0 * pi * d_lum**2 * tau)
         Fssc(j, i) = D**4 * volume * freqs(j) * jssc(j, i) / (4d0 * pi * d_lum**2)
         Feic(j, i) = D**4 * volume * freqs(j) * jeic(j, i) / (4d0 * pi * d_lum**2)
         Fnut(j, i) = Fmbs(j, i) + Fssc(j, i) + Feic(j, i)
      end do

      Ntot(i) = sum(nn(i, 2:) * (gg(2:) - gg(:numbins - 1)))
      if ( mod(i, nmod) == 0 .or. i == 1 ) &
         write(*, on_screen) i, t(i), nu_ext, nu0(i, numbins), Ntot(i)

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
   call h5io_wdble0(group_id, 'tstep', tstep, herror)
   call h5io_wdble0(group_id, 'sigma', sigma, herror)
   call h5io_wdble0(group_id, 'R', R, herror)
   call h5io_wdble0(group_id, 'R0', R0, herror)
   call h5io_wdble0(group_id, 'd_lum', d_lum, herror)
   call h5io_wdble0(group_id, 'redshift', z, herror)
   call h5io_wdble0(group_id, 'gamma_bulk', gamma_bulk, herror)
   call h5io_wdble0(group_id, 'view-angle', theta_obs, herror)
   call h5io_wdble0(group_id, 'gamma_min', gmin, herror)
   call h5io_wdble0(group_id, 'gamma_max', gmax, herror)
   call h5io_wdble0(group_id, 'gamma_1', g1, herror)
   call h5io_wdble0(group_id, 'gamma_2', g2, herror)
   call h5io_wdble0(group_id, 'pwl-index', qind, herror)
   call h5io_wdble0(group_id, 'L_j', L_j, herror)
   call h5io_wdble0(group_id, 'Theta_e', theta_e, herror)
   call h5io_wdble0(group_id, 'zeta_e', zetae, herror)
   call h5io_wdble0(group_id, 'nu_min', numin, herror)
   call h5io_wdble0(group_id, 'nu_max', numax, herror)
   call h5io_closeg(group_id, herror)

   ! ------  Saving data  ------
   call h5io_wdble1(file_id, 'time', t(1:), herror)
   call h5io_wdble1(file_id, 't_obs', t_obs, herror)
   call h5io_wdble1(file_id, 'Ntot', Ntot, herror)
   call h5io_wdble1(file_id, 'sen_lum', sen_lum, herror)
   call h5io_wdble1(file_id, 'frequency', freqs, herror)
   call h5io_wdble1(file_id, 'nu_obs', nu_obs, herror)
   call h5io_wdble1(file_id, 'gamma', gg, herror)

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
   call h5io_wdble2(file_id, 'distrib', nn(1:, :), herror)
   call h5io_wdble2(file_id, 'nu0_tot', nu0(1:, :), herror)

   ! ------  Closing output file  ------
   call h5io_closef(file_id, herror)
   call h5close_f(herror)

   write(*, "('==========  FINISHED  ==========')")
   write(*,*) ''

end subroutine afterglow