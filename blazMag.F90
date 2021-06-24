subroutine blazMag(params_file, output_file, cool_withKN, with_abs)
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
   implicit none

   character(len=*), intent(in) :: output_file, params_file
   logical, intent(in) :: cool_withKN, with_abs
   integer, parameter :: nmod = 50
   character(len=*), parameter :: screan_head = &
      '| Iteration |        Time |   Time step |    nu_0(g2) |       N_tot |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------',&
      on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
#ifdef HDF5
   integer(HID_T) :: file_id, group_id
   integer :: herror
#endif
   integer :: i, j, k, numbins, numdf, numdt, time_grid
   real(dp) :: uB, uext, R, gmin, gmax, numin, numax, pind, B, D, g1, g2, &
         tstep, Qnth, tmax, d_lum, z, tinj, gamma_bulk, theta_obs, Rdis, &
         mu_obs, nu_ext, tesc, tlc, mu_mag, L_jet, volume, sigma, beta_bulk, L_B, &
         eps_B, f_rec, urad_const, f_esc, eps_acc
   real(dp) :: L_e, Qnth2, Qnth3
   real(dp), allocatable, dimension(:) :: freqs, t, Ntot, Inu, gg, dt, nu_obs, &
         t_obs, dg, urad
   real(dp), allocatable, dimension(:,:) :: dotg, nn, jnut, jmbs, jssc, jeic, &
         ambs, anut, Qinj, Ddiff
   logical :: with_cool


   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   call read_params(params_file)
   d_lum = par_d_lum
   z = par_z
   tstep = par_tstep
   tmax = par_tmax
   eps_B = par_eps_B
   eps_acc = par_eps_acc
   f_rec = par_frec
   L_jet = par_L_j
   gamma_bulk = par_gamma_bulk
   sigma = par_sigma
   ! mu_mag = par_mu_mag
   gmin = par_gmin
   gmax = par_gmax
   pind = par_pind
   numin = par_numin
   numax = par_numax
   numbins = par_numbins
   numdt = par_numdt
   numdf = par_numdf
   time_grid = par_time_grid
   f_esc = par_fesc


   allocate(t(0:numdt), freqs(numdf), Ntot(numdt), Inu(numdf), &
      dt(numdt), nu_obs(numdf), t_obs(numdt), dg(numbins), urad(numbins))
   allocate(nn(numbins, 0:numdt), dotg(numbins, 0:numdt), gg(numbins), &
      ambs(numdf, numdt), jmbs(numdf, numdt), jnut(numdf, numdt), &
      jssc(numdf, numdt), anut(numdf, numdt), jeic(numdf, numdt), &
      Qinj(numbins, numdt), Ddiff(numbins, 0:numdt))


   !   # #    # # #####     ####   ####  #    # #####
   !   # ##   # #   #      #    # #    # ##   # #    #
   !   # # #  # #   #      #      #    # # #  # #    #
   !   # #  # # #   #      #      #    # #  # # #    #
   !   # #   ## #   #      #    # #    # #   ## #    #
   !   # #    # #   #       ####   ####  #    # #####
   with_cool = .false.

   ! TODO: try mu_mag as input
   ! sigma = (mu_mag / gamma_bulk) - 1d0
   mu_mag = gamma_bulk * (sigma + 1d0)    ! baryon loading

   beta_bulk = bofg(gamma_bulk)           ! Beta bulk
   theta_obs = par_theta_obs * pi / 180d0 ! Observer viewing angle
   mu_obs = dcos(theta_obs)
   D = Doppler(gamma_bulk, mu_obs)        ! Doppler factor
   R = par_R                              ! Radius of the blob

   ! ----->    External radiation field
   nu_ext = par_nu_ext * gamma_bulk
   uext = par_uext * gamma_bulk**2 * (1d0 + beta_bulk**2 / 3d0) ! Eq. (5.25) Dermer & Menon (2009)
   urad_const = 4d0 * sigmaT * cLight / (3d0 * energy_e)

   ! ----->    Magnetic field
   L_B = sigma * L_jet / (1d0 + sigma)
   uB = L_B / (2d0 * pi * R**2 * cLight * beta_bulk * gamma_bulk * (gamma_bulk - 1d0))
   B = dsqrt(uB * 8d0 * pi)

   ! ----->   Injection of particles
   if ( pind > 2d0 ) then
      g1 = f_rec * sigma * mass_p * (pind - 2d0) / ((pind - 1d0) * mass_e)
      ! g2 = par_g2
      g2 = dsqrt(6d0 * pi * eCharge * eps_acc / (sigmaT * B))
   else if ( pind > 1d0 .and. pind < 2d0 ) then
      g1 = par_g1
      ! g2 = dsqrt(6d0 * pi * eCharge * 1d-3 / (sigmaT * B))
      g2 = ( f_rec * (sigma + 1d0) * mass_p * (2d0 - pind) * g1**(1d0 - pind) / (mass_e * (pind - 1d0)) )**(1d0 / (2d0 - pind))
   else
      g1 = par_g1
      g2 = par_g2
!       g2 = dsqrt(6d0 * pi * eCharge * 1d-3 / (sigmaT * B))
   end if

   volume = 4d0 * pi * R**3 / 3d0
   tlc = R / cLight
   tesc = f_esc * tlc
   tinj = tlc

   L_e = f_rec * L_B / (1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0))
!   Qnth = f_rec * uB * pwl_norm(tlc * energy_e, pind - 1d0, g1, g2)
   Qnth = f_rec * L_B * pwl_norm(1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0) * volume * energy_e, pind - 1d0, g1, g2)
   Qnth2 = f_rec * L_B * pwl_norm(volume * energy_e, pind - 1d0, g1, g2)
   Qnth3 = f_rec * L_B / ( (g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6)) * 1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0) * volume * energy_e )

   ! ----->   Output with initial setup
   write(*, "('--> Simulation setup')")
   write(*, "('theta_obs =', ES15.7)") par_theta_obs
   write(*, "('Doppler   =', ES15.7)") D
   write(*, "('gamma_1   =', ES15.7)") g1
   write(*, "('gamma_2   =', ES15.7)") g2
   write(*, "('L_jet     =', ES15.7)") L_jet
   write(*, "('L_B       =', ES15.7)") L_B
   write(*, "('L_e       =', ES15.7)") L_e
   write(*, "('L_e,2     =', ES15.7)") Qnth * volume * energy_e
   write(*, "('L_e,3     =', ES15.7)") Qnth2 * volume * energy_e
   write(*, "('Q_nth     =', ES15.7)") Qnth
   write(*, "('Q_nth_old =', ES15.7)") Qnth2
   write(*, "('Q_nth_alt =', ES15.7)") Qnth3
   write(*, "('u_B       =', ES15.7)") uB
   write(*, "('B         =', ES15.7)") B
   write(*, "('u_ext     =', ES15.7)") par_uext
   write(*, "('nu_ext    =', ES15.7)") par_nu_ext
   write(*, "('sigma     =', ES15.7)") sigma
   write(*, "('mu        =', ES15.7)") mu_mag
   write(*, "('Gamma     =', ES15.7)") gamma_bulk
   write(*, "('t_dyn     =', ES15.7)") tlc
   write(*, "('t_esc     =', ES15.7)") tesc
   write(*, "('t_inj     =', ES15.7)") tinj
   write(*, "('R_b       =', ES15.7)") R

   build_f: do j=1, numdf
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numdf - 1)) )
      freqs(j) = nu_com_f(nu_obs(j), z, D)
   end do build_f

   build_g: do k = 1, numbins
      gg(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numbins - 1))
      ! gg(k) = (gmin - 1d0) * ((gmax - 1d0) / (gmin - 1d0))**(dble(k - 1) / dble(numbins - 1)) + 1d0
      if ( k > 1 ) dg(k) = gg(k) - gg(k - 1)
   end do build_g
   dg(1) = dg(2)

   t(0) = 0d0
   dotg(:, 0) = urad_const * (uB + uext) * pofg(gg)**2
   Ddiff(:, 0) = 1d-200

   Qinj(:, 0) = injection_pwl(0d0, tinj, gg, g1, g2, pind, Qnth3)
   nn(:, 0) = Qinj(:, 0)

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

      t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numdt - 1)) )
      dt(i) = t(i) - t(i - 1)
      t_obs(i) = (1d0 + z) * t(i) / D


      !  ###### ###### #####
      !  #      #      #    #
      !  #####  #####  #    #
      !  #      #      #    #
      !  #      #      #    #
      !  ###### ###### #####
      call FP_FinDif_difu(dt(i), &
            &             gg, &
            &             nn(:, i - 1), &
            &             nn(:, i), &
            &             dotg(:, i - 1), &
            &             Ddiff(:, i - 1), &
            &             Qinj(:, i - 1), &
            &             tesc, &
            &             tlc)

      Qinj(:, i) = injection_pwl(t(i), tinj, gg, g1, g2, pind, Qnth)
      Ddiff(:, i) = 1d-200


      !  #####    ##   #####  #   ##   ##### #  ####  #    #
      !  #    #  #  #  #    # #  #  #    #   # #    # ##   #
      !  #    # #    # #    # # #    #   #   # #    # # #  #
      !  #####  ###### #    # # ######   #   # #    # #  # #
      !  #   #  #    # #    # # #    #   #   # #    # #   ##
      !  #    # #    # #####  # #    #   #   #  ####  #    #
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numdf
         call syn_emissivity(jmbs(j, i), freqs(j), gg, nn(:, i), B)
         call syn_absorption(ambs(j, i), freqs(j), gg, nn(:, i), B)
      end do
      !$OMP END PARALLEL DO

      call RadTrans_blob(Inu, R, jmbs(:, i), ambs(:, i))

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numdf
         call IC_iso_powlaw(jssc(j, i), freqs(j), freqs, Inu, nn(:, i), gg)
         call IC_iso_monochrom(jeic(j, i), freqs(j), uext, nu_ext, nn(:, i), gg)
         jnut(j, i) = jmbs(j, i) + jssc(j, i) + jeic(j, i)
         anut(j, i) = ambs(j, i)
         ! Fmbs(j, i) = D**4 * volume * freqs(j) * jmbs(j, i) * opt_depth_blob(anut(j, i), R) / (4d0 * pi * d_lum**2)
         ! Fssc(j, i) = D**4 * volume * freqs(j) * jssc(j, i) * opt_depth_blob(anut(j, i), R) / (4d0 * pi * d_lum**2)
         ! Feic(j, i) = D**4 * volume * freqs(j) * jeic(j, i) * opt_depth_blob(anut(j, i), R) / (4d0 * pi * d_lum**2)
         ! Fnut(j, i) = Fmbs(j, i) + Fssc(j, i) + Feic(j, i)
      end do
      !$OMP END PARALLEL DO


      !   ####   ####   ####  #      # #    #  ####
      !  #    # #    # #    # #      # ##   # #    #
      !  #      #    # #    # #      # # #  # #
      !  #      #    # #    # #      # #  # # #  ###
      !  #    # #    # #    # #      # #   ## #    #
      !   ####   ####   ####  ###### # #    #  ####
      dotg(:, i) = 0d0
      if ( with_cool ) then
         ! call bolometric_integ(freqs, 4d0 * pi * Inu / cLight, urad)
         ! call RadTrans_blob(Inu, R, jssc(:, i) + jeic(:, i), anut(:, i))
         call RadTrans_blob(Inu, R, jmbs(:, i), ambs(:, i))
         call rad_cool(dotg(:, i), gg, freqs, 4d0 * pi * Inu / cLight, cool_withKN)
      end if
      dotg(:, i) = dotg(:, i) + urad_const * (uB + uext) * pofg(gg)**2


      ! ----->   N_tot
      Ntot(i) = sum(nn(:, i) * dg, mask=nn(:, i) > 1d-200)

      if ( mod(i, nmod) == 0 .or. i == 1 ) &
         write(*, on_screen) i, t(i), dt(i), dotg(numbins, i), Ntot(i)

   end do time_loop


   !  ####    ##   #    # # #    #  ####
   ! #       #  #  #    # # ##   # #    #
   !  ####  #    # #    # # # #  # #
   !      # ###### #    # # #  # # #  ###
   ! #    # #    #  #  #  # #   ## #    #
   !  ####  #    #   ##   # #    #  ####
#ifdef HDF5
   write(*, *) "---> Saving"
   !----->   Opening output file
   call h5open_f(herror)
   call h5io_createf(output_file, file_id, herror)
   !----->   Saving initial parameters
   call h5io_createg(file_id, "Parameters", group_id, herror)
   call h5io_wint0(group_id, 'numdt', numdt, herror)
   call h5io_wint0(group_id, 'numdf', numdf, herror)
   call h5io_wint0(group_id, 'numbins', numbins, herror)
   call h5io_wdble0(group_id, 't_max', tmax, herror)
   call h5io_wdble0(group_id, 'tstep', tstep, herror)
   call h5io_wdble0(group_id, 'R_b', R, herror)
   call h5io_wdble0(group_id, 'R_em', Rdis, herror)
   call h5io_wdble0(group_id, 'd_lum', d_lum, herror)
   call h5io_wdble0(group_id, 'redshift', z, herror)
   call h5io_wdble0(group_id, 'Gamma_bulk', gamma_bulk, herror)
   call h5io_wdble0(group_id, 'sigma', sigma, herror)
   call h5io_wdble0(group_id, 'theta_obs_deg', par_theta_obs, herror)
   call h5io_wdble0(group_id, 'gamma_min', gmin, herror)
   call h5io_wdble0(group_id, 'gamma_max', gmax, herror)
   call h5io_wdble0(group_id, 'gamma_1', g1, herror)
   call h5io_wdble0(group_id, 'gamma_2', g2, herror)
   call h5io_wdble0(group_id, 'pwl-index', pind, herror)
   call h5io_wdble0(group_id, 'u_ext', par_uext, herror)
   call h5io_wdble0(group_id, 'nu_ext', par_nu_ext, herror)
   call h5io_wdble0(group_id, 'L_jet', L_jet, herror)
   call h5io_wdble0(group_id, 'nu_min', numin, herror)
   call h5io_wdble0(group_id, 'nu_max', numax, herror)
   call h5io_wdble0(group_id, 'mu_mag', mu_mag, herror)
   call h5io_closeg(group_id, herror)
   !----->   Saving data
   call h5io_wdble0(file_id, 't_inj', tinj, herror)
   call h5io_wdble0(file_id, 't_esc', tesc, herror)

   call h5io_createg(file_id, "electrons-energy", group_id, herror)
   call h5io_wdble0(group_id, 'Bfield', B, herror)
   call h5io_wdble0(group_id, 'L_B', L_B, herror)          ! Eq. (8)
   call h5io_wdble0(group_id, 'uB', uB, herror)            ! Eq. (9)
   call h5io_wdble0(group_id, 'L_e', L_e, herror)          ! Eq. (10)
   call h5io_wdble0(group_id, 'L_e2', f_rec * L_B, herror) ! old, wrong L_e
   call h5io_wdble0(group_id, 'Q_nth', Qnth, herror)       ! Eq. (13)
   call h5io_wdble0(group_id, 'Q_nth2', Qnth2, herror)     ! old, wrong Q0
   call h5io_wdble0(group_id, 'Q_nth3', Qnth3, herror)     ! from full expresion of Eq. (12)
   call h5io_closeg(group_id, herror)
   
   call h5io_wdble1(file_id, 'time', t(1:), herror)
   call h5io_wdble1(file_id, 't_obs', t_obs, herror)
   call h5io_wdble1(file_id, 'nu', freqs, herror)
   call h5io_wdble1(file_id, 'nu_obs', nu_obs, herror)
   call h5io_wdble1(file_id, 'gamma', gg, herror)
   call h5io_wdble2(file_id, 'jnut', jnut, herror)
   call h5io_wdble2(file_id, 'jmbs', jmbs, herror)
   call h5io_wdble2(file_id, 'jssc', jssc, herror)
   call h5io_wdble2(file_id, 'jeic', jeic, herror)
   call h5io_wdble2(file_id, 'anut', anut, herror)
   call h5io_wdble2(file_id, 'ambs', ambs, herror)
   call h5io_wdble2(file_id, 'n_e', nn(:, 1:), herror)
   call h5io_wdble2(file_id, 'dgdt', dotg(:, 1:), herror)
   !----->   Closing output file
   call h5io_closef(file_id, herror)
   call h5close_f(herror)
#endif
   write(*,*) '=======  FINISHED  ======='
   write(*,*) ''

end subroutine blazMag
