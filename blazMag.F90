subroutine blazMag(params_file, output_file, with_cool, with_abs, with_ssc)
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
   implicit none

   character(len=*), intent(in) :: output_file, params_file
   logical, intent(in) :: with_cool, with_abs, with_ssc
   integer, parameter :: nmod = 50
   character(len=*), parameter :: screan_head = &
      '| Iteration |        Time |   Time step |    nu_0(g2) |       N_tot |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------',&
      on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
   integer(HID_T) :: file_id, group_id
   integer :: i, j, k, numbins, numdf, numdt, time_grid, herror
   real(dp) :: uB, uext, R, L_j, gmin, gmax, numin, numax, qind, B, D, &
      tacc, g1, g2, tstep, Qnth, tmax, d_lum, z, tvar, tinj,&
      gamma_bulk, theta_obs, R0, mu_obs, nu_ext, tesc, tlc, &
      volume, sigma, beta_bulk, eps_e, L_B, eps_B, f_rec
   real(dp), allocatable, dimension(:) :: freqs, t, Ntot, Inu, gg, &
      dt, nu_obs, t_obs, dg, urad
   real(dp), allocatable, dimension(:, :) :: nu0, nn, jnut, jmbs, jssc, jeic, &
      ambs, anut, Qinj, Ddif, Fmbs, Feic, Fssc, Fnut

   call read_params(params_file)
   ! R = par_R
   ! R0 = par_R0
   d_lum = par_d_lum
   z = par_z
   tstep = par_tstep
   tmax = par_tmax
   L_j = par_L_j
   eps_B = par_eps_B
   f_rec = par_frec
   sigma = par_sigma
   if ( par_mu_mag > 1d0 ) then
      gamma_bulk = par_mu_mag / (1d0 + sigma)
   else
      gamma_bulk = par_gamma_bulk
   end if
   gmin = par_gmin
   gmax = par_gmax
   qind = par_qind
   uext = par_uext
   numin = par_numin
   numax = par_numax
   numbins = par_numbins
   numdt = par_numdt
   numdf = par_numdf
   time_grid = par_time_grid
   tvar =  par_tvar



   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   allocate(t(0:numdt), freqs(numdf), Ntot(numdt), Inu(numdf), &
      dt(numdt), nu_obs(numdf), t_obs(numdt), dg(numbins), urad(numbins))
   allocate(nn(numbins, 0:numdt), nu0(numbins, 0:numdt), gg(numbins), &
      ambs(numdf, numdt), jmbs(numdf, numdt), jnut(numdf, numdt), &
      jssc(numdf, numdt), anut(numdf, numdt), jeic(numdf, numdt), &
      Qinj(numbins, numdt), Ddif(numbins, numdt), Fnut(numdf, numdt), &
      Fmbs(numdf, numdt), Fssc(numdf, numdt), Feic(numdf, numdt))


   !   # #    # # #####     ####   ####  #    # #####
   !   # ##   # #   #      #    # #    # ##   # #    #
   !   # # #  # #   #      #      #    # # #  # #    #
   !   # #  # # #   #      #      #    # #  # # #    #
   !   # #   ## #   #      #    # #    # #   ## #    #
   !   # #    # #   #       ####   ####  #    # #####
   theta_obs = 1d0 / gamma_bulk ! par_theta_obs * pi / 180d0
   beta_bulk = bofg(gamma_bulk)
   mu_obs = dcos(theta_obs)
   D = Doppler(gamma_bulk, mu_obs)
   R = 0.95d0 * D * tvar * cLight / (1d0 + z)

   ! ----->    External radiation field
   nu_ext = par_nu_ext * gamma_bulk !nu_com_f(par_nu_ext, z, D)
   uext = uext * gamma_bulk**2 * (1d0 + beta_bulk**2 / 3d0) !  Eq. (5.25) Dermer & Menon (2009)

   ! ----->    Magnetic field
   L_B = sigma * L_j / (1d0 + sigma)
   uB = L_B / (pi * cLight * beta_bulk * (gamma_bulk * R)**2) ! B**2 / (8d0 * pi)
   B = dsqrt(uB * 8d0 * pi)
   nu0 = 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)

   ! ----->   Injection of particles
   if ( qind > 2d0 ) then
      g1 = (qind - 2d0) * f_rec * sigma * mass_p / ((qind - 1d0) * mass_e)
      ! g2 = par_g2
      g2 = dsqrt(6d0 * pi * eCharge * 1d-6 / (sigmaT * B))
   else if ( qind > 1d0 .and. qind < 2d0 ) then
      g1 = par_g1
      ! g2 = dsqrt(6d0 * pi * eCharge * 1d-3 / (sigmaT * B))
      g2 = ((sigma + 1d0) * (mass_p / mass_e) * ((2d0 - qind) / (qind - 1d0)))**(1d0 / (2d0 - qind)) * g1**((1d0 - qind) / (2d0 - qind))
   else
      g1 = par_g1
      g2 = par_g2
      ! g2 = dsqrt(6d0 * pi * eCharge * 1d-3 / (sigmaT * B))
   end if

   volume = 4d0 * pi * R**3 / 3d0
   tlc = R / cLight
   tesc = 1.5d0 * tlc
   tacc = 0.5d0 * tlc
   tinj = tesc


   Qnth = f_rec * L_B * pwl_norm(volume * mass_e * cLight**2, qind - 1d0, g1, g2)
   ! Qnth = eps_e * (L_j - L_B) * pwl_norm(volume * mass_e * cLight**2, qind - 1d0, g1, g2)

   write(*, "('--> Simulation setup')")
   ! write(*, *) ''
   write(*, "('theta_obs =', ES15.7)") theta_obs
   write(*, "('Doppler   =', ES15.7)") D
   write(*, "('gamma_1   =', ES15.7)") g1
   write(*, "('gamma_2   =', ES15.7)") g2
   write(*, "('Q_nth     =', ES15.7)") Qnth
   write(*, "('t_dyn     =', ES15.7)") tlc
   write(*, "('t_esc     =', ES15.7)") tesc
   write(*, "('t_acc     =', ES15.7)") tacc
   write(*, "('L_B       =', ES15.7)") L_B
   write(*, "('u_B       =', ES15.7)") uB
   write(*, "('B         =', ES15.7)") B
   write(*, "('u_ext     =', ES15.7)") uext
   write(*, "('nu_ext    =', ES15.7)") nu_ext
   write(*, "('mu        =', ES15.7)") par_mu_mag
   write(*, "('Gamma     =', ES15.7)") gamma_bulk

   build_f: do j=1,numdf
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numdf - 1)) )
      freqs(j) = nu_com_f(nu_obs(j), z, D)
   end do build_f

   build_g: do k = 1, numbins
      ! gg(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numbins - 1))
      gg(k) = (gmin - 1d0) * ((gmax - 1d0) / (gmin - 1d0))**(dble(k - 1) / dble(numbins - 1)) + 1d0
      if ( k > 1 ) dg(k) = gg(k) - gg(k - 1)
   end do build_g
   dg(1) = dg(2)

   t(0) = 0d0
   nn(:, 0) = injection_pwl(0d0, 1d200, gg, g1, g2, qind, Qnth)
   jnut = 0d0

   write(*, "('--> Calculating the emission')")
   write(*, *) ''
   write(*, "('Using tstep = ', F5.3)") tstep
   write(*, "('Wrting data in: ', A)") trim(output_file)
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
      case(2)
         t(i) = t(i - 1) + 1d0 / (nu0(i - 1, numbins) * g2)
      case(3)
         t(i) = tmax * dble(i) / dble(numdt)
      case default
         write(*, *) "Wrong time-grid selection"
      end select
      dt(i) = t(i) - t(i - 1)
      t_obs(i) = (1d0 + z) * t(i) / D

      !  #####    ##   #####  #   ##   ##### #  ####  #    #
      !  #    #  #  #  #    # #  #  #    #   # #    # ##   #
      !  #    # #    # #    # # #    #   #   # #    # # #  #
      !  #####  ###### #    # # ######   #   # #    # #  # #
      !  #   #  #    # #    # # #    #   #   # #    # #   ##
      !  #    # #    # #####  # #    #   #   #  ####  #    #
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numdf
         call mbs_emissivity(jmbs(j, i), freqs(j), gg, nn(:, i - 1), B)
         if ( with_abs ) then
            call mbs_absorption(ambs(j, i), freqs(j), gg, nn(:, i - 1), B)
         else
            ambs(j, i) = 0d0
         end if
         call RadTrans_blob(Inu, R, jmbs(:, i), ambs(:, i))
         ! Inu(j) = opt_depth_blob(ambs(j, i), R) * jmbs(j, i) * 2d0 * R
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numdf
         if ( with_ssc ) then
            call IC_iso_powlaw(jssc(j, i), freqs(j), freqs, Inu, nn(:, i - 1), gg)
            call IC_iso_monochrom(jeic(j, i), freqs(j), uext, nu_ext, nn(:, i - 1), gg)
         else
            jeic(j, i) = 0d0
            jssc(j, i) = 0d0
         end if
         jnut(j, i) = jmbs(j, i) + jssc(j, i) + jeic(j, i)
         anut(j, i) = ambs(j, i)
         Fmbs(j, i) = D**4 * volume * freqs(j) * jmbs(j, i) * opt_depth_blob(anut(j, i), R) / d_lum**2
         Fssc(j, i) = D**4 * volume * freqs(j) * jssc(j, i) * opt_depth_blob(anut(j, i), R) / d_lum**2
         Feic(j, i) = D**4 * volume * freqs(j) * jeic(j, i) * opt_depth_blob(anut(j, i), R) / d_lum**2
         ! Fssc(j, i) = 0.25d0 * D**4 * volume * freqs(j) * jssc(j, i) / (pi * d_lum**2)
         ! Feic(j, i) = 0.25d0 * D**4 * volume * freqs(j) * jeic(j, i) / (pi * d_lum**2)
         Fnut(j, i) = Fmbs(j, i) + Fssc(j, i) + Feic(j, i)
         ! Fnut(j, i) = D**4 * volume * freqs(j) * jnut(j, i) * opt_depth_blob(anut(j, i), R) / d_lum**2
      end do
      !$OMP END PARALLEL DO

      !   ####   ####   ####  #      # #    #  ####
      !  #    # #    # #    # #      # ##   # #    #
      !  #      #    # #    # #      # # #  # #
      !  #      #    # #    # #      # #  # # #  ###
      !  #    # #    # #    # #      # #   ## #    #
      !   ####   ####   ####  ###### # #    #  ####
!      if ( b_index /= 0d0 ) then
!         uB = 0.125d0 * (B * (1d0 + (cLight * gamma_bulk * t(i) / R0))**(-b_index))**2 / pi
!      end if

      if ( with_cool ) then
         urad = bolometric_integ(freqs, 4d0 * pi * Inu / cLight)
         call RadTrans_blob(Inu, R, jssc(:, i) + jeic(:, i), anut(:, i))
         urad = urad + IC_cool(gg, freqs, 4d0 * pi * Inu / cLight)
      else
         urad = 0d0
      end if

      !  ###### ###### #####
      !  #      #      #    #
      !  #####  #####  #    #
      !  #      #      #    #
      !  #      #      #    #
      !  ###### ###### #####
      Qinj(:, i) = injection_pwl(t(i), tacc, gg, g1, g2, qind, Qnth)
      Ddif(:, i) = 1d-200 !0.5d0 * gg**2 / tacc !
      nu0(:, i) = 4d0 * sigmaT * (uB + urad) / (3d0 * mass_e * cLight)
      call FP_FinDif_difu(dt(i), &
            &             gg, &
            &             nn(:, i - 1), &
            &             nn(:, i), &
            &             nu0(:, i) * (gg**2 - 1d0), &
            &             Ddif(:, i), &
            &             Qinj(:, i), &
            &             tesc)
      ! call FP_FinDif_cool(dt(i), gg, nn(i - 1, :), nn(i, :), nu0(i - 1, :), Qinj(i, :), tesc)


      ! ----->   N_tot
      Ntot(i) = sum(nn(:, i) * dg, mask=nn(:, i) > 1d-200)

      if ( mod(i, nmod) == 0 .or. i == 1 ) &
         write(*, on_screen) i, t(i), dt(i), nu0(numbins, i), Ntot(i)

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
   call h5io_wdble0(group_id, 'Gamma_bulk', gamma_bulk, herror)
   call h5io_wdble0(group_id, 'view-angle', theta_obs, herror)
   call h5io_wdble0(group_id, 'gamma_min', gmin, herror)
   call h5io_wdble0(group_id, 'gamma_max', gmax, herror)
   call h5io_wdble0(group_id, 'gamma_1', g1, herror)
   call h5io_wdble0(group_id, 'gamma_2', g2, herror)
   call h5io_wdble0(group_id, 'pwl-index', qind, herror)
   call h5io_wdble0(group_id, 'L_j', L_j, herror)
   call h5io_wdble0(group_id, 'nu_min', numin, herror)
   call h5io_wdble0(group_id, 'nu_max', numax, herror)
   call h5io_wdble0(group_id, 'mu_mag', par_mu_mag, herror)
   call h5io_closeg(group_id, herror)

   ! ------  Saving data  ------
   call h5io_wdble0(file_id, 't_acc', tacc, herror)
   call h5io_wdble0(file_id, 't_esc', tesc, herror)
   call h5io_wdble0(file_id, 'Bfield', B, herror)
   call h5io_wdble0(file_id, 'Q_nth', Qnth, herror)

   call h5io_wdble1(file_id, 'time', t(1:), herror)
   call h5io_wdble1(file_id, 't_obs', t_obs, herror)
   call h5io_wdble1(file_id, 'Ntot', Ntot, herror)
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
   call h5io_wdble2(file_id, 'n_e', nn(:, 1:), herror)
   call h5io_wdble2(file_id, 'nu0_tot', nu0(:, 1:), herror)

   ! ------  Closing output file  ------
   call h5io_closef(file_id, herror)
   call h5close_f(herror)

   write(*,*) '=======  FINISHED  ======='
   write(*,*) ''

end subroutine blazMag
