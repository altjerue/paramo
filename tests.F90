program tests
   use data_types
   use constants
   use misc
   use pwl_integ
   use hdf5
   use h5_inout
   use SRtoolkit
   use anaFormulae
   use dist_evol
   use K1
   use K2
   implicit none

   call steady_state

   write(*,*) '=======  FINISHED  ======='
   write(*,*) ''
   
contains
   
   subroutine steady_state
      implicit none
      integer(HID_T) :: file_id
      integer :: i, k, numg, numt, herror
      real(dp) :: g1, g2, gmin, gmax, tmax, tstep, qind, tacc, tesc, R, &
         theta_e, zeta_e
      real(dp), allocatable, dimension(:) :: t, g, Q0, D0, Ntot, dt, dg, &
         C0
      real(dp), allocatable, dimension(:, :) :: n
      character(*), parameter :: dir = "/Users/jesus/lab/2018/blazMag/output/"

      numg = 192
      numt = 200
      g1 = 1e1
      g2 = 1e6
      gmin = 1.0001d0
      gmax = 1.1d0 * g2
      tstep = 1e0
      tmax = 5e7
      qind = 0d0
      R = 1e16
      theta_e = 2e1
      zeta_e = 1d0

      allocate(g(numg), Q0(numg), D0(numg), t(0:numt), dt(numt), dg(numg), &
         Ntot(numt), C0(numg), n(0:numt, numg))

      build_g: do k = 1, numg
         g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numg - 1))
         if ( k > 1 ) dg(k) = g(k) - g(k - 1)
         ! if ( g(k) >= g1 .and. g(k) <= g2 ) then
         !    n(0, k) = 1d0
         ! else
         !    n(0, k) = 0d0
         ! end if
      end do build_g
      dg(1) = dg(2)

      t(0) = 0d0
      C0 = 3.48d-11 !4d0 * sigmaT * uB / (3d0 * mass_e * cLight)
      tacc = 1d0 / (C0(1) * 10d0**4.5d0) !tesc
      tesc = tacc ! 1d200! 1.5d0 * R / cLight!
      D0 = 0.5d0 * g**2 / tacc
      write(*, *) C0(1), D0(1), tacc
      n(0, :) = injection(1d0, tacc, g, g1, g2, 0d0, theta_e, 0d0, 1d0)
      Q0 = injection(1d0, tacc, g, g1, g2, 0d0, theta_e, 0d0, 1d0)

      time_loop: do i = 1, numt
         
         t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numt - 1)) )
         dt(i) = t(i) - t(i - 1)
         
         call FP_FinDif_difu(dt(i), g, n(i - 1, :), n(i, :), C0, D0, Q0, tesc)
         ! call FP_FinDif_cool(dt(i), g, n(i - 1, :), n(i, :), C0, Q0, tesc)

         Ntot(i) = sum(n(i, :) * dg, mask=n(i, :) > 1d-200)

      end do time_loop

      call h5open_f(herror)
      call h5io_createf(dir//"SSsol.h5", file_id, herror)
      call h5io_wint0(file_id, 'numdt', numt, herror)
      call h5io_wint0(file_id, 'numbins', numg, herror)
      call h5io_wdble0(file_id, 't_max', tmax, herror)
      call h5io_wdble0(file_id, 'tstep', tstep, herror)
      call h5io_wdble0(file_id, 'gamma_min', gmin, herror)
      call h5io_wdble0(file_id, 'gamma_max', gmax, herror)
      call h5io_wdble0(file_id, 'gamma_1', g1, herror)
      call h5io_wdble0(file_id, 'gamma_2', g2, herror)
      call h5io_wdble0(file_id, 'pwl-index', qind, herror)
      call h5io_wdble0(file_id, 'Theta_e', theta_e, herror)
      call h5io_wdble0(file_id, 'zeta_e', zeta_e, herror)
      call h5io_wdble0(file_id, 't_acc', tacc, herror)
      call h5io_wdble0(file_id, 't_esc', tesc, herror)
      call h5io_wdble1(file_id, 'time', t(1:), herror)
      call h5io_wdble1(file_id, 'Ntot', Ntot, herror)
      call h5io_wdble1(file_id, 'gamma', g, herror)
      call h5io_wdble2(file_id, 'distrib', n(1:, :), herror)
      call h5io_closef(file_id, herror)
      call h5close_f(herror)

   end subroutine steady_state









#if 0
   subroutine no_diff_CG99
      
   end subroutine no_diff_CG99
   
   
   subroutine IC_catastroph_PPM14
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
      integer :: i, j, k, numbins, numdf, numdt, ios, time_grid, herror
      real(dp) :: uB, uext, urad, R, L_j, gmin, gmax, numin, numax, qind, B, &
      tacc, g1, g2, tstep, zetae, Qth, Qnth, theta_e, tmax, d_lum, z, D, &
      gamma_bulk, theta_obs, R0, b_index, mu_obs, mu_com, nu_ext, tesc, kappa, &
      volume, sigma, beta_bulk, eps_e, L_B, mu_mag, Iind
      real(dp), allocatable, dimension(:) :: freqs, t, Ntot, Inu, gg, sen_lum, &
      dt, nu_obs, t_obs, dg
      real(dp), allocatable, dimension(:, :) :: nu0, nn, jnut, jmbs, jssc, jeic, &
      ambs, anut, Qinj, Ddif
      
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
      Qinj(numdt, numbins), Ddif(numdt, numbins))
      
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
      nu0 = 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)
      
      ! --->   Injection of particles
      volume = 4d0 * pi * R**3 / 3d0
      tesc = 1.5d0 * R / cLight
      tacc = 1d0 / (nu0(1, 1) * 10d0**4.5d0)!tesc
      if ( hyb_dis ) then
         kappa = 3d0 * theta_e + dexp(K1_func(-dlog(theta_e))) / dexp(K2_func(-dlog(theta_e)))
         Qth = (1d0 - zetae) * (L_j - L_B) / (kappa * volume * mass_e * cLight**2)
         Qnth = zetae * Qth * pwl_norm(1d0 - zetae, qind, g1, g2)
      else
         kappa = 1d0
         Qth = 0d0
         Qnth = 0d0 !eps_e * (L_j - L_B) * pwl_norm(volume * mass_e * cLight**2, qind, g1, g2)
      end if
      
      ! ----->>   Viewing angle
      mu_obs = dcos(theta_obs * pi / 180d0)
      mu_com = mu_com_f(gamma_bulk, mu_obs)
      D = Doppler(gamma_bulk, mu_obs)
      
      build_f: do j=1,numdf
         nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numdf - 1)) )
         freqs(j) = nu_com_f(nu_obs(j), z, gamma_bulk, mu_com)
      end do build_f
      
      t(0) = 0d0
      jnut = 0d0
      
      time_loop: do i = 1, numdt
         
         t_obs(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numdt - 1)) )
         t(i) = D * t_obs(i) / (1d0 + z)
         dt(i) = t(i) - t(i - 1)
         
         call mbs_emissivity(jmbs(:, i), freqs, gg, nn(i, :), B)
         call mbs_absorption(ambs(:, i), freqs, gg, nn(i, :), B)
         call RadTrans_blob(Inu, R, jmbs(:, i), ambs(:, i))
         call SSC_pwlEED(jssc(:, i), freqs, Inu, nn(i, :), gg)
         call EIC_pwlEED(jeic(:, i), freqs, uext, nu_ext, nn(i, :), gg)
         
         if ( b_index /= 0d0 ) then
            uB = 0.125 * (B * (1d0 + (cLight * gamma_bulk * t(i) / R0))**(-b_index))**2 / pi
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
      end do time_loop
      
      
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
      call h5io_wdble2(file_id, 'distrib', nn(1:, :), herror)
      call h5io_wdble2(file_id, 'nu0_tot', nu0(1:, :), herror)
      
      ! ------  Closing output file  ------
      call h5io_closef(file_id, herror)
      call h5close_f(herror)
      
   end subroutine IC_catastroph_PPM14

#endif


end program tests
