subroutine afterglow(params_file, output_file, with_cool, with_ic)
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
   use models
   use K1
   use K2
   implicit none

   character(len=*), intent(in) :: output_file, params_file
   logical, intent(in) :: with_cool, with_ic
   integer, parameter :: nmod = 50
   character(len=*), parameter :: screan_head = &
      '| Iteration |        Time |   Time step |  gamma_bulk |       N_tot |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------', &
      on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
   integer(HID_T) :: file_id, group_id
   integer :: i, j, k, numbins, numdf, numdt, time_grid, herror
   real(dp) :: uB, uext, L_j, gmin, gmax, numin, numax, qind, B, R0, &
      tacc, tinj, g1, g2, tstep, Q0, tmax, d_lum, z, n_ext, &
      theta_obs, mu_obs, nu_ext, tesc_e, uext0, volume, eps_e, tlc,&
      eps_B, E0, gamma_bulk0, L_e, nu_ext0, Aadi, tmin, Rd, td, R, dr, &
      b_index, beta_bulk, volume_p, R_p, eps_g2
   real(dp), allocatable, dimension(:) :: freqs, t, Ntot, Isyn, gg, sen_lum, &
      dt, nu_obs, t_obs, gamma_bulk, Rbw, D, tcool, urad
   real(dp), allocatable, dimension(:, :) :: nu0, n_e, jnut, jmbs, jssc, jeic, &
      ambs, anut, Qinj, Ddif, Fmbs, Feic, Fssc, Fnut

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
   qind = par_qind
   numin = par_numin
   numax = par_numax
   nu_ext0 = par_nu_ext
   uext0 = par_uext
   numbins = par_numbins
   numdt = par_numdt
   numdf = par_numdf
   E0 = par_E0
   n_ext = par_n_ext
   g1 = par_g1
   g2 = par_g2
   time_grid = par_time_grid
   R0 = par_R0
   tmax = par_tmax
   b_index = par_b_index


   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   write(*, "('--> Simulation setup')")

   allocate(t(0:numdt), freqs(numdf), Ntot(numdt), Isyn(numdf), sen_lum(numdt), &
      dt(numdt), nu_obs(numdf), t_obs(0:numdt), Rbw(0:numdt), tcool(numbins), &
      gamma_bulk(0:numdt), D(0:numdt), urad(numbins))
   allocate(n_e(numbins, 0:numdt), nu0(numbins, 0:numdt), gg(numbins), &
      ambs(numdf, numdt), jmbs(numdf, numdt), jnut(numdf, numdt), &
      jssc(numdf, numdt), anut(numdf, numdt), jeic(numdf, numdt), &
      Qinj(numbins, 0:numdt), Ddif(numbins, 0:numdt), Fnut(numdf, numdt), &
      Fmbs(numdf, numdt), Fssc(numdf, numdt), Feic(numdf, numdt))

   call K1_init
   call K2_init

   !
   !   ---------->    Initial conditions     <----------
   !
   theta_obs = par_theta_obs * pi / 180d0
   mu_obs = dcos(theta_obs)
   beta_bulk = bofg(gamma_bulk0)
   gamma_bulk(0) = gamma_bulk0
   D(0) = Doppler(gamma_bulk0, mu_obs)
   Rd = deceleration_radius(E0, gamma_bulk0, n_ext)
   Rbw(0) = R0
   R = R0 / gamma_bulk0
   R_p = R
   volume = 4d0 * pi * R**3 / 3d0
   volume_p = volume
   B = dsqrt(32d0 * pi * eps_B * mass_p * n_ext) * cLight * gamma_bulk0
   ! B = dsqrt(32d0 * pi * eps_B * mass_p * n_ext) * cLight * 0.5d0 * gamma_bulk0 * Rd**(1.5d0) *  R0**(-b_index)
   uB = 0.125d0 * B**2 / pi
   eps_g2 = 1d-1
   g2 = dsqrt(6d0 * pi * eCharge * eps_g2 / (sigmaT * B))
   g1 = eps_e * (gamma_bulk0 - 1d0) * mass_p * (qind - 2d0) / ((qind - 1d0) * mass_e)
   ! L_e = eps_e * 4d0 * pi * R0**2 * n_ext * mass_p * cLight**3 * beta_bulk * gamma_bulk0 * (gamma_bulk0 - 1d0)
   L_e = eps_e * 2d0 * pi * (1d0 - dcos(1d0 / gamma_bulk0)) * R0**2 * n_ext * mass_p * cLight**3 * beta_bulk * gamma_bulk0 * (gamma_bulk0 - 1d0)
   Q0 = L_e * pwl_norm(mass_e * cLight**2 * volume, qind - 1d0, g1, g2)
   td = (1d0 + z) * Rd / (beta_bulk * cLight * gamma_bulk0**2)
   t(0) = 0d0
   t_obs(0) = 0d0
   tlc = R / cLight
   ! tinj = tlc
   ! tesc_e = 1.5d0 * tlc
   ! tacc = tesc_e

   write(*, "('theta_obs =', ES15.7)") theta_obs
   write(*, "('nu_ext0   =', ES15.7)") nu_ext0
   write(*, "('Doppler   =', ES15.7)") D(0)
   write(*, "('td        =', ES15.7)") td
   write(*, "('Rd        =', ES15.7)") Rd
   write(*, "('Gamma0    =', ES15.7)") gamma_bulk0
   write(*, "('R         =', ES15.7)") R

   build_f: do j = 1, numdf
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numdf - 1)) )
   end do build_f

   build_g: do k = 1, numbins
      gg(k) = (gmin - 1d0) * ((gmax - 1d0) / (gmin - 1d0))**(dble(k - 1) / dble(numbins - 1)) + 1d0
   end do build_g

   nu0 = 0d0
   Aadi = 0d0
   Ddif(:, 0) = 1d-200 ! 0.5d0 * gg**2 / tacc !
   n_e(:, 0) = injection_pwl(0d0, 1d200, gg, g1, g2, qind, Q0)

   write(*, "('---> Calculating the emission')")
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
   time_loop: do i = 1, numdt

      select case(time_grid)
      case(1)
         t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numdt - 1)) )
         dt(i) = t(i) - t(i - 1)
      case(2)
         call time_step(dt(i), gg, Ddif(:, i - 1), nu0(:, i - 1) * gg**2 + Aadi * gg, tstep, tlc)
         t(i) = t(i - 1) + dt(i)
      case(3)
         t(i) = tmax * dble(i) / dble(numdt)
         dt(i) = t(i) - t(i - 1)
      case default
         write(*, *) "Wrong time-grid selection"
      end select
      dr = dt(i) * beta_bulk * gamma_bulk(i - 1) * cLight
      Rbw(i) = Rbw(i - 1) + dr
      call adiab_blast_wave(Rbw(i), R0, gamma_bulk0, E0, n_ext, gamma_bulk(i))
      beta_bulk = bofg(gamma_bulk(i))
      D(i) = Doppler(gamma_bulk(i), mu_obs)
      t_obs(i) = t_obs(i - 1) + 0.5d0 * (1d0 + z) * dt(i) * ( 1d0 / D(i) + 1d0 / D(i - 1) )
      R = Rbw(i) / gamma_bulk(i)
      tlc = R / cLight
      volume = 4d0 * pi * R**3 / 3d0
      B = dsqrt(32d0 * pi * eps_B * mass_p * n_ext) * cLight * gamma_bulk(i)
      ! B = dsqrt(32d0 * pi * eps_B * mass_p * n_ext) * cLight * 0.5d0 * gamma_bulk0 * Rd**(1.5d0) *  Rbw(i)**(-b_index)
      uB = 0.125d0 * B**2 / pi
      g2 = dsqrt(6d0 * pi * eCharge * eps_g2 / (sigmaT * B))
      g1 = eps_e * (gamma_bulk(i) - 1d0) * mass_p * (qind - 2d0) / ((qind - 1d0) * mass_e)
      ! L_e = eps_e * 4d0 * pi * Rbw(i)**2 * n_ext * mass_p * cLight**3 * beta_bulk * gamma_bulk(i) * (gamma_bulk(i) - 1d0)
      L_e = eps_e * 2d0 * pi * (1d0 - dcos(1d0 / gamma_bulk(i))) * Rbw(i)**2 * n_ext * mass_p * cLight**3 * beta_bulk * gamma_bulk(i) * (gamma_bulk(i) - 1d0)
      Q0 = L_e * pwl_norm(mass_e * cLight**2 * volume, qind - 1d0, g1, g2)
      uext = uext0 * gamma_bulk(i)**2 * (1d0 + beta_bulk**2 / 3d0)
      nu_ext = nu_ext0 * gamma_bulk(i)

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numdf
         freqs(j) = nu_com_f(nu_obs(j), z, D(i))
         call mbs_emissivity(jmbs(j, i), freqs(j), gg, n_e(:, i - 1), B)
         call mbs_absorption(ambs(j, i), freqs(j), gg, n_e(:, i - 1), B)
         Isyn(j) = opt_depth_blob(ambs(j, i), R) * jmbs(j, i) * 2d0 * R
         ! Isyn(j) = opt_depth_slab(ambs(j, i), 2d0 * R) * jmbs(j, i) * 2d0 * R
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numdf
         if ( with_ic ) then
            call IC_iso_powlaw(jssc(j, i), freqs(j), freqs, Isyn, n_e(:, i - 1), gg)
            call IC_iso_monochrom(jeic(j, i), freqs(j), uext, nu_ext, n_e(:, i - 1), gg)
         else
            jssc(j, i) = 0d0
            jeic(j, i) = 0d0
         end if
         jnut(j, i) = jmbs(j, i) + jssc(j, i) + jeic(j, i)
         anut(j, i) = ambs(j, i)

         Fmbs(j, i) = D(i)**4 * volume * freqs(j) * jmbs(j, i) * opt_depth_blob(ambs(j, i), R) / d_lum**2
         Fssc(j, i) = D(i)**4 * volume * freqs(j) * jssc(j, i) / (4d0 * pi * d_lum**2)
         Feic(j, i) = D(i)**4 * volume * freqs(j) * jeic(j, i) / (4d0 * pi * d_lum**2)
         Fnut(j, i) = Fmbs(j, i) + Fssc(j, i) + Feic(j, i)
      end do
      !$OMP END PARALLEL DO

      if ( with_cool ) then
         urad = bolometric_integ(freqs, 4d0 * pi * Isyn / cLight)
         urad = urad + uext
         call RadTrans_blob(Isyn, R, jssc(:, i) + jeic(:, i), anut(:, i))
         urad = urad + IC_cool(gg, freqs, 4d0 * pi * Isyn / cLight)
      else
         urad = 0d0
      end if

      Aadi = 3d0 * beta_bulk * gamma_bulk(i) * cLight / (2d0 * Rbw(i))
      nu0(:, i) = 4d0 * sigmaT * (uB + urad) / (3d0 * mass_e * cLight)
      ! tesc_e = 1.5d0 * tlc
      ! tacc = tesc_e
      ! tcool = 1d0 / (nu0(:, i) * pofg(gg) + Aadi)
      ! tadiab = 1d0 / (Aadi * gg(ig) * bofg(gg(ig)))
      Qinj(:, i) = injection_pwl(t(i), 1d200, gg, g1, g2, qind, Q0)
      Ddif(:, i) = 1d-200 !0.5d0 * gg**2 / tacc !
      call FP_FinDif_difu(dt(i), &
            &             pofg(gg), &
            &             n_e(:, i - 1), &
            &             n_e(:, i), &
            &             nu0(:, i) * pofg(gg)**2 + Aadi * pofg(gg), &
            &             Ddif(:, i), &
            &             Qinj(:, i), &
            ! &             tesc_e)
            &             1d200)

      volume_p = volume
      R_p = R

      Ntot(i) = 0.5d0 * sum((n_e(:numbins - 1, i) + n_e(2:, i)) * (gg(2:) - gg(:numbins - 1)))

      if ( mod(i, nmod) == 0 .or. i == 1 ) &
         write(*, on_screen) i, t(i), dt(i), gamma_bulk(i), Ntot(i)

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
   call h5io_wdble0(group_id, 'pwl-index', qind, herror)
   call h5io_wdble0(group_id, 'L_j', L_j, herror)
   call h5io_wdble0(group_id, 'epsilon_e', eps_e, herror)
   call h5io_wdble0(group_id, 'epsilon_B', eps_B, herror)
   call h5io_wdble0(group_id, 'nu_min', numin, herror)
   call h5io_wdble0(group_id, 'nu_max', numax, herror)
   call h5io_wdble0(group_id, 'nu_ext0', nu_ext0, herror)
   call h5io_wdble0(group_id, 'u_ext0', uext0, herror)
   call h5io_wdble0(group_id, 'E0', E0, herror)
   call h5io_wdble0(group_id, 'n_ext', n_ext, herror)

   call h5io_closeg(group_id, herror)

   ! ------  Saving data  ------
   call h5io_wdble0(file_id, 'Rd', Rd, herror)
   call h5io_wdble0(file_id, 'td', td, herror)

   call h5io_wdble1(file_id, 'time', t(1:), herror)
   call h5io_wdble1(file_id, 't_obs', t_obs(1:), herror)
   call h5io_wdble1(file_id, 'Ntot', Ntot, herror)
   call h5io_wdble1(file_id, 'Rbw', Rbw(1:), herror)
   call h5io_wdble1(file_id, 'Gamma_bulk', gamma_bulk(1:), herror)
   call h5io_wdble1(file_id, 'nu', freqs, herror)
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
   call h5io_wdble2(file_id, 'Qinj', Qinj, herror)
   call h5io_wdble2(file_id, 'n_e', n_e(:, 1:), herror)
   call h5io_wdble2(file_id, 'nu0_tot', nu0(:, 1:), herror)

   ! ------  Closing output file  ------
   call h5io_closef(file_id, herror)
   call h5close_f(herror)

   write(*, "('==========  FINISHED  ==========')")
   write(*,*) ''

end subroutine afterglow
