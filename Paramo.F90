subroutine Paramo(params_file, output_file, hyb_dis, with_cool, with_abs, with_ssc)
   use data_types
   use constants
   use misc
   use pwl_integ
   use hdf5
   use h5_inout
   use SRtoolkit
   use anaFormulae
   use dist_evol, only: FP_FinDif_difu, FP_FinDif_cool, injection
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
      ' ---------------------------------------------------------------------', &
      on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
   integer(HID_T) :: file_id, group_id
   integer :: i, j, k, numbins, numdf, numdt, ios, time_grid, herror, tstop
   real(dp) :: uB, uext, urad, R, L_j, gmin, gmax, numin, numax, qind, B, &
      tacc, g1, g2, tstep, zetae, Qth, Qnth, theta_e, tmax, d_lum, z, D, &
      gamma_bulk, theta_obs, R0, b_index, mu_obs, mu_com, nu_ext, tesc, kappa, &
      volume, sigma, beta_bulk, eps_e, L_B, mu_mag, Iind
   real(dp), allocatable, dimension(:) :: freqs, t, dg, Ntot, Inu, gg, &
      sen_lum, dfreqs, dtimes, dt, nu_obs, t_obs, to_com
   real(dp), allocatable, dimension(:, :) :: nu0, nn, jnut, jmbs, &
      jssc, jeic, ambs, anut, Qinj, Ddif

   open(unit=77, file=params_file, iostat=ios)
   if ( ios /= 0 ) call an_error("Paramo: Parameter file "//params_file//" could not be opened")
   read(77, *) R
   read(77, *) R0
   read(77, *) d_lum
   read(77, *) z
   read(77, *) gamma_bulk
   read(77, *) theta_obs
   read(77, *) sigma
   read(77, *) b_index
   read(77, *) theta_e
   read(77, *) zetae
   read(77, *) tstep
   read(77, *) tmax
   read(77, *) L_j
   read(77, *) eps_e
   read(77, *) g1
   read(77, *) g2
   read(77, *) gmin
   read(77, *) gmax
   read(77, *) qind
   read(77, *) nu_ext
   read(77, *) uext
   read(77, *) numin
   read(77, *) numax
   read(77, *) numbins
   read(77, *) numdt
   read(77, *) numdf
   read(77, *) time_grid
   close(77)

   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   allocate(t(0:numdt), freqs(numdf), dg(numbins), Ntot(numdt),&
      dfreqs(numdf), dtimes(numdt), Inu(numdf), &
      sen_lum(numdt), dt(numdt), nu_obs(numdf), t_obs(numdt), to_com(numdt))
   allocate(nn(0:numdt, numbins), nu0(0:numdt, numbins), gg(numbins), &
      ambs(numdf, numdt), jmbs(numdf, numdt), jnut(numdf, numdt), &
      jssc(numdf, numdt), anut(numdf, numdt), jeic(numdf, numdt), &
      Qinj(numdt, numbins), Ddif(numdt, numbins))

   dfreqs = 1d0
   dtimes = 1d0

   call K1_init
   call K2_init

   !   # #    # # #####     ####   ####  #    # #####
   !   # ##   # #   #      #    # #    # ##   # #    #
   !   # # #  # #   #      #      #    # # #  # #    #
   !   # #  # # #   #      #      #    # #  # # #    #
   !   # #   ## #   #      #    # #    # #   ## #    #
   !   # #    # #   #       ####   ####  #    # #####
   beta_bulk = bofg(gamma_bulk)

   ! --->    External radiation field
   nu_ext = nu_ext * gamma_bulk
   uext = uext * gamma_bulk**2

   ! --->    Magnetic field
   L_B = sigma * L_j / (1d0 + sigma)
   uB = L_B / (pi * cLight * beta_bulk * (gamma_bulk * R)**2) ! B**2 / (8d0 * pi)
   B = dsqrt(uB * 8d0 * pi)
   mu_mag = gamma_bulk * (1d0 + sigma)

   volume = 4d0 * pi * R**3 / 3d0
   tesc = 1.5 * R / cLight
   tacc = tesc
   if ( hyb_dis ) then
      kappa = 3d0 * theta_e + dexp(K1_func(-dlog(theta_e))) / dexp(K2_func(-dlog(theta_e)))
      Qth = (1d0 - zetae) * (L_j - L_B) / (kappa * volume * mass_e * cLight**2)
      Qnth = zetae * Qth * pwl_norm(1d0 - zetae, qind, g1, g2)
   else
      kappa = 1d0
      Qth = 0d0
      Qnth = eps_e * (L_j - L_B) * pwl_norm(volume * mass_e * cLight**2, qind, g1, g2)
   end if

   write(*, *) ' ---> Simulation setup'
   write(*, *) ''
   write(*, "('g1    =', ES15.7)") g1
   write(*, "('Q_nth =', ES15.7)") Qnth
   write(*, "('Q_th  =', ES15.7)") Qth
   write(*, "('t_esc =', ES15.7)") tesc
   write(*, "('t_acc =', ES15.7)") tacc
   write(*, "('B     =', ES15.7)") B
   write(*, "('mu    =', ES15.7)") mu_mag


   ! ----->>   Viewing angle
   mu_obs = dcos(theta_obs * pi / 180d0)
   mu_com = mu_com_f(gamma_bulk, mu_obs)
   D = Doppler(gamma_bulk, mu_obs)

   build_f: do j=1,numdf
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numdf - 1)) )
      freqs(j) = nu_com_f(nu_obs(j), z, gamma_bulk, mu_com)
   end do build_f
   dfreqs(1:numdf - 1) = 1d0 / (freqs(2:numdf) - freqs(1:numdf - 1))

   build_g: do k=1, numbins
      gg(k) = gmin * ( (gmax / gmin)**(dble(k - 1) / dble(numbins - 1)) )
      if ( k > 1 ) dg(k - 1) = gg(k) - gg(k - 1)
   end do build_g
   dg(numbins) = dg(numbins - 1)

   t(0) = 0d0
   nn = 0d0
   jnut = 0d0

   write(*, *) '---> Calculating the emission'
   write(*, *) ''
   write(*, "(' Using tstep = ', F5.3)") tstep
   write(*, *) 'Wrting data in: ', trim(output_file)
   write(*, *) ''
   write(*, *) screan_head

   !
   ! ###### #    #  ####  #      #    # ##### #  ####  #    #
   ! #      #    # #    # #      #    #   #   # #    # ##   #
   ! #####  #    # #    # #      #    #   #   # #    # # #  #
   ! #      #    # #    # #      #    #   #   # #    # #  # #
   ! #       #  #  #    # #      #    #   #   # #    # #   ##
   ! ######   ##    ####  ######  ####    #   #  ####  #    #
   !
   tstop = numdt
   time_loop: do i = 1, numdt

      select case(time_grid)
      case(1)
         t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numdt - 1)) )
      case(2)
         if ( i == 1 ) then
            t(i) = t(i - 1) + dmin1(tstep / (nu0(max0(1, i - 1), numbins) * g2), 1e-2 * tstep)
         else
            t(i) = t(i - 1) + tstep / (nu0(max0(1, i - 1), numbins) * g2)
         end if
         if ( t(i - 1) >= tmax ) then
            tstop = i - 1
            exit time_loop
         end if
      case(3)
         t(i) = (tmax - tstep) * dble(i) / dble(numdt - 1)
      case default
         write(*, *) "Wrong time-grid selection"
      end select
      dt(i) = t(i) - t(i - 1)
      dtimes(i) = 1d0 / dt(i)


      !  ###### ###### #####
      !  #      #      #    #
      !  #####  #####  #    #
      !  #      #      #    #
      !  #      #      #    #
      !  ###### ###### #####
      Qinj(i, :) = injection(t(i), tacc, gg, g1, g2, qind, theta_e, Qth, Qnth)
      Ddif(i, :) = 1d-200 ! gg**2 / tesc !
      call FP_FinDif_difu(dt(i), gg, nn(i - 1, :), nn(i, :), nu0(i - 1, :), Ddif(i, :), Qinj(i, :), tesc)
      ! call FP_FinDif_cool(dt(i), gg, nn(i - 1, :), nn(i, :), nu0(i - 1, :), Qinj(i, :), tmax)
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
         uB = 0.125 * (B * (1d0 + (cLight * gamma_bulk * t(i) / R0))**(-b_index))**2 / pi
      end if

      if ( with_cool ) then
         urad = 0d0
         do j = 2, numdf
            if ( Inu(j - 1) > 1d-200 .and. Inu(j) > 1d-200) then
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


      ! #     #
      ! ##    #         #####  ####  #####
      ! # #   #           #   #    #   #
      ! #  #  #           #   #    #   #
      ! #   # #           #   #    #   #
      ! #    ##           #   #    #   #
      ! #     #           #    ####    #
      !         #######
      Ntot(i) = sum(nn(i, :) * dg, mask=(nn(i, :) > 1d-200))

      t_obs(i) = t(i)

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
   call h5io_wint0(file_id, 't_stop', tstop, herror)

   call h5io_wdble0(file_id, 't_acc', tacc, herror)
   call h5io_wdble0(file_id, 't_esc', tesc, herror)
   call h5io_wdble0(file_id, 'Bfield', B, herror)
   call h5io_wdble0(file_id, 'mu', mu_mag, herror)
   call h5io_wdble0(file_id, 'Q_th', Qth, herror)
   call h5io_wdble0(file_id, 'Q_nth', Qnth, herror)

   call h5io_wdble1(file_id, 'time', t(1:tstop), herror)
   call h5io_wdble1(file_id, 't_obs', t_obs(:tstop), herror)
   call h5io_wdble1(file_id, 'Ntot', Ntot(:tstop), herror)
   call h5io_wdble1(file_id, 'sen_lum', sen_lum(:tstop), herror)
   call h5io_wdble1(file_id, 'frequency', freqs, herror)
   call h5io_wdble1(file_id, 'nu_obs', nu_obs, herror)
   call h5io_wdble1(file_id, 'gamma', gg, herror)

   call h5io_wdble2(file_id, 'jnut', jnut(:, :tstop), herror)
   call h5io_wdble2(file_id, 'jmbs', jmbs(:, :tstop), herror)
   call h5io_wdble2(file_id, 'jssc', jssc(:, :tstop), herror)
   call h5io_wdble2(file_id, 'jeic', jeic(:, :tstop), herror)
   call h5io_wdble2(file_id, 'anut', anut(:, :tstop), herror)
   call h5io_wdble2(file_id, 'ambs', ambs(:, :tstop), herror)
   call h5io_wdble2(file_id, 'distrib', nn(1:tstop, :), herror)
   call h5io_wdble2(file_id, 'nu0_tot', nu0(1:tstop, :), herror)

   ! ------  Closing output file  ------
   call h5io_closef(file_id, herror)
   call h5close_f(herror)

   write(*,*) '=======  FINISHED  ======='
   write(*,*) ''

end subroutine Paramo
