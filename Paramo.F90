subroutine Paramo(params_file, output_file, hyb_dis, with_cool, with_abs, with_ssc)
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
   use K1
   use K2
   implicit none

   character(len=*), intent(in) :: output_file, params_file
   logical, intent(in) :: with_cool, with_abs, with_ssc, hyb_dis
   integer, parameter :: nmod = 10
   character(len=*), parameter :: screan_head = &
      '| Iteration |        Time |   Time step |    nu_0(g2) |       N_tot |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------',&
      on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
   integer(HID_T) :: file_id, group_id
   integer :: i, j, k, numbins, numdf, numdt, time_grid, herror
   real(dp) :: uB, uext, urad, R, L_j, gmin, gmax, numin, numax, qind, B, &
      tacc, g1, g2, tstep, zetae, Qth, Qnth, theta_e, tmax, d_lum, z, D, &
      gamma_bulk, theta_obs, R0, b_index, mu_obs, mu_com, nu_ext, tesc, kappa, &
      volume, sigma, beta_bulk, eps_e, L_B, mu_mag, Iind, eps_B, tau
   real(dp), allocatable, dimension(:) :: freqs, t, Ntot, Inu, gg, sen_lum, &
      dt, nu_obs, t_obs, dg
   real(dp), allocatable, dimension(:, :) :: nu0, nn, jnut, jmbs, jssc, jeic, &
      ambs, anut, Qinj, Ddif, flux

   call read_params(params_file)
   R = par_R
   R0 = par_R0
   d_lum = par_d_lum
   z = par_z
   gamma_bulk = par_gamma_bulk
   theta_obs = par_theta_obs
   zetae = par_zetae
   tstep = par_tstep
   tmax = par_tmax
   L_j = par_L_j
   eps_e = par_eps_e
   eps_B = par_eps_B
   sigma = par_sigma
   g1 = par_g1
   g2 = par_g2
   gmin = par_gmin
   gmax = par_gmax
   qind = par_qind
   nu_ext = par_nu_ext
   uext = par_uext
   numin = par_numin
   numax = par_numax
   numbins = par_numbins
   numdt = par_numdt
   numdf = par_numdf
   time_grid = par_time_grid


   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   allocate(t(0:numdt), freqs(numdf), Ntot(numdt), Inu(numdf), sen_lum(numdt), &
      dt(numdt), nu_obs(numdf), t_obs(numdt), dg(numbins))
   allocate(nn(0:numdt, numbins), nu0(0:numdt, numbins), gg(numbins), &
      ambs(numdf, numdt), jmbs(numdf, numdt), jnut(numdf, numdt), &
      jssc(numdf, numdt), anut(numdf, numdt), jeic(numdf, numdt), &
      Qinj(numdt, numbins), Ddif(numdt, numbins), flux(numdf, numdt))

   call K1_init
   call K2_init

   !   # #    # # #####     ####   ####  #    # #####
   !   # ##   # #   #      #    # #    # ##   # #    #
   !   # # #  # #   #      #      #    # # #  # #    #
   !   # #  # # #   #      #      #    # #  # # #    #
   !   # #   ## #   #      #    # #    # #   ## #    #
   !   # #    # #   #       ####   ####  #    # #####
   beta_bulk = bofg(gamma_bulk)
   mu_obs = dcos(theta_obs * pi / 180d0)
   mu_com = mu_com_f(gamma_bulk, mu_obs)
   D = Doppler(gamma_bulk, mu_obs)

   ! --->    External radiation field
   nu_ext = nu_ext / D ! * gamma_bulk
   uext = uext * gamma_bulk**2 * (1d0 + beta_bulk**2 / 3d0) !  Eq. (5.25) Dermer & Menon (2009)

   ! --->    Magnetic field
   L_B = sigma * L_j / (1d0 + sigma)
   uB = L_B / (pi * cLight * beta_bulk * (gamma_bulk * R)**2) ! B**2 / (8d0 * pi)
   B = dsqrt(uB * 8d0 * pi)
   mu_mag = gamma_bulk * (1d0 + sigma)
   nu0 = 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)

   ! --->   Injection of particles
   volume = 4d0 * pi * R**3 / 3d0
   tesc = 1.5d0 * R / cLight
   tacc = 0.5d0 * R / cLight
   if ( hyb_dis ) then
      kappa = 3d0 * theta_e + dexp(K1_func(-dlog(theta_e))) / dexp(K2_func(-dlog(theta_e)))
      Qth = (1d0 - zetae) * (L_j - L_B) / (kappa * volume * mass_e * cLight**2)
      Qnth = zetae * Qth * pwl_norm(1d0 - zetae, qind, g1, g2)
   else
      kappa = 1d0
      Qth = 0d0
      ! Qnth = eps_e * (L_j - L_B) * pwl_norm(volume * mass_e * cLight**2, qind - 1d0, g1, g2)
      Qnth = eps_e * L_B * pwl_norm(volume * mass_e * cLight**2, qind - 1d0, g1, g2)
   end if

   write(*, *) ' ---> Simulation setup'
   write(*, *) ''
   write(*, "('Doppler =', ES15.7)") D
   write(*, "('Q_nth   =', ES15.7)") Qnth
   write(*, "('Q_th    =', ES15.7)") Qth
   write(*, "('t_dyn   =', ES15.7)") R / cLight
   write(*, "('t_esc   =', ES15.7)") tesc
   write(*, "('t_acc   =', ES15.7)") tacc
   write(*, "('L_B     =', ES15.7)") L_B
   write(*, "('u_B     =', ES15.7)") uB
   write(*, "('B       =', ES15.7)") B
   write(*, "('mu      =', ES15.7)") mu_mag


   build_f: do j=1,numdf
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numdf - 1)) )
      freqs(j) = nu_com_f(nu_obs(j), z, gamma_bulk, mu_obs)
   end do build_f

   build_g: do k = 1, numbins
      gg(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numbins - 1))
      if ( k > 1 ) dg(k) = gg(k) - gg(k - 1)
   end do build_g
   dg(1) = dg(2)

   t(0) = 0d0
   nn = 0d0
   jnut = 0d0

   write(*, *) '---> Calculating the emission'
   write(*, *) ''
   write(*, "(' Using tstep = ', F5.3)") tstep
   write(*, *) 'Wrting data in: ', trim(output_file)
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
         t_obs(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numdt - 1)) )
      case(2)
         t_obs(i) = t_obs(i - 1) + tstep / (nu0(max0(1, i - 1), numbins) * g2)
      case(3)
         t_obs(i) = tmax * dble(i) / dble(numdt)
      case default
         write(*, *) "Wrong time-grid selection"
      end select
      t(i) = D * t_obs(i) / (1d0 + z)
      dt(i) = t(i) - t(i - 1)


      !  ###### ###### #####
      !  #      #      #    #
      !  #####  #####  #    #
      !  #      #      #    #
      !  #      #      #    #
      !  ###### ###### #####
      Qinj(i, :) = injection(t(i), tacc, gg, g1, g2, qind, theta_e, Qth, Qnth)
      Ddif(i, :) = 1d-200 !0.5d0 * gg**2 / tacc !
      call FP_FinDif_difu(dt(i), gg, nn(i - 1, :), nn(i, :), nu0(i - 1, :), Ddif(i, :), Qinj(i, :), tesc)
      ! call FP_FinDif_cool(dt(i), gg, nn(i - 1, :), nn(i, :), nu0(i - 1, :), Qinj(i, :), tesc)
      !   ----->   Then we compute the light path
      sen_lum(i) = sum(dt(:i)) * cLight


      !  #####    ##   #####  #   ##   ##### #  ####  #    #
      !  #    #  #  #  #    # #  #  #    #   # #    # ##   #
      !  #    # #    # #    # # #    #   #   # #    # # #  #
      !  #####  ###### #    # # ######   #   # #    # #  # #
      !  #   #  #    # #    # # #    #   #   # #    # #   ##
      !  #    # #    # #####  # #    #   #   #  ####  #    #
      !
      ! ----->>   magnetobrem   <<-----
      call mbs_emissivity(jmbs(:, i), freqs, gg, nn(i, :), B)

      ! ----->>   Absorption   <<-----
      if ( with_abs ) then
         call mbs_absorption(ambs(:, i), freqs, gg, nn(i, :), B)
      else
         ambs(:, i) = 0d0
      end if

      ! ----->>   Inverse Compton   <<-----
      call RadTrans_blob(Inu, R, jmbs(:, i), ambs(:, i))
      if ( with_ssc ) then
         call SSC_pwlEED(jssc(:, i), freqs, Inu, nn(i, :), gg)
         call EIC_pwlEED(jeic(:, i), freqs, uext, nu_ext, nn(i, :), gg)
      else
         jeic(:, i) = 0d0
         jssc(:, i) = 0d0
      end if

      ! ----->>   Totals   <<-----
      jnut(:, i) = jmbs(:, i) + jssc(:, i) + jeic(:, i)
      anut(:, i) = ambs(:, i)


      !   ####   ####   ####  #      # #    #  ####
      !  #    # #    # #    # #      # ##   # #    #
      !  #      #    # #    # #      # # #  # #
      !  #      #    # #    # #      # #  # # #  ###
      !  #    # #    # #    # #      # #   ## #    #
      !   ####   ####   ####  ###### # #    #  ####
      if ( b_index /= 0d0 ) then
         uB = 0.125d0 * (B * (1d0 + (cLight * gamma_bulk * t(i) / R0))**(-b_index))**2 / pi
      end if

      if ( with_cool ) then
         urad = 0d0
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

      ! ----->   N_tot
      Ntot(i) = sum(nn(i, :) * dg, mask=nn(i, :) > 1d-200)

      ! ----->   nu F_nu
      do j = 1, numdf
         tau = dmax1(1d-200, 2d0 * R * ambs(j, i))
         flux(j, i) = D**4 * volume * ( 3d0 * freqs(j) * opt_depth_blob(tau) * jmbs(j, i) / tau + freqs(j) * (jssc(j, i) + jeic(j, i)) ) / (4d0 * pi * d_lum**2)
      end do

      if ( mod(i, nmod) == 0 .or. i == 1 ) &
         write(*, on_screen) i, t(i), dt(i), nu0(i, numbins), Ntot(i)

   end do time_loop


   !  ####    ##   #    # # #    #  ####
   ! #       #  #  #    # # ##   # #    #
   !  ####  #    # #    # # # #  # #
   !      # ###### #    # # #  # # #  ###
   ! #    # #    #  #  #  # #   ## #    #
   !  ####  #    #   ##   # #    #  ####
   write(*, *) "---> Saving"

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
   call h5io_wdble0(group_id, 'Bfield-index', b_index, herror)
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
   call h5io_wdble0(file_id, 't_acc', tacc, herror)
   call h5io_wdble0(file_id, 't_esc', tesc, herror)
   call h5io_wdble0(file_id, 'Bfield', B, herror)
   call h5io_wdble0(file_id, 'mu', mu_mag, herror)
   call h5io_wdble0(file_id, 'Q_th', Qth, herror)
   call h5io_wdble0(file_id, 'Q_nth', Qnth, herror)

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
   call h5io_wdble2(file_id, 'flux', flux, herror)
   call h5io_wdble2(file_id, 'distrib', nn(1:, :), herror)
   call h5io_wdble2(file_id, 'nu0_tot', nu0(1:, :), herror)

   ! ------  Closing output file  ------
   call h5io_closef(file_id, herror)
   call h5close_f(herror)

   write(*,*) '=======  FINISHED  ======='
   write(*,*) ''

end subroutine Paramo
