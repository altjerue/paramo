subroutine Paramo(params_file, output_file, with_cool, with_abs, with_ssc)
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
   implicit none
   character(len=*), intent(in) :: output_file, params_file
   logical, intent(in) :: with_cool, with_abs, with_ssc
   integer, parameter :: nmod = 10
   character(len=*), parameter :: screan_head = &
      '| Iteration |        Time |   Time step |    nu_0(g2) |       N_tot |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------', &
      on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
   integer :: i, j, k, numbins, numdf, numdt, ios, time_grid, &
      cool_kind, herror, tstop
   integer(HID_T) :: file_id, group_id
   real(dp) :: uB, R, Q0, gmin, gmax, numin, numax, qind, B, dtacc, g1, &
      g2, nu0_B, tstep, zetae, Qth, Qnth, theta_e, tmax, d_lum, z, D, &
      gamma_bulk, theta_obs, R0, Rinit, b_index, mu_obs, mu_com, u_ext, &
      nu_ext, tesc
   real(dp), allocatable, dimension(:) :: freqs, t, dg, Ntot, Imbs, gg, &
      aux_zero_arr, sen_lum, dfreqs, dtimes, dt, nu_obs, t_obs, to_com
   real(dp), allocatable, dimension(:, :) :: nu0, nn, jnut, jmbs, &
      jssc, jeic, ambs, anut, Qinj, Ddif

   open(unit=77, file=params_file, iostat=ios)
   if ( ios /= 0 ) call an_error("Paramo: Parameter file "//params_file//" could not be opened")
   read(77,*) R
   read(77,*) R0
   read(77,*) Rinit
   read(77,*) d_lum
   read(77,*) z
   read(77,*) gamma_bulk
   read(77,*) theta_obs
   read(77,*) B
   read(77,*) b_index
   read(77,*) theta_e
   read(77,*) dtacc
   read(77,*) tstep
   read(77,*) tmax
   read(77,*) Q0
   read(77,*) g1
   read(77,*) g2
   read(77,*) gmin
   read(77,*) gmax
   read(77,*) qind
   read(77,*) numin
   read(77,*) numax
   read(77,*) numbins
   read(77,*) numdt
   read(77,*) numdf
   read(77,*) cool_kind
   read(77,*) time_grid
   close(77)

   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   allocate(t(0:numdt), freqs(numdf), dg(numbins), Ntot(numdt),&
      aux_zero_arr(numbins - 1), dfreqs(numdf), dtimes(numdt), Imbs(numdf), &
      sen_lum(numdt), dt(numdt), nu_obs(numdf), t_obs(numdt), to_com(numdt))
   allocate(nn(0:numdt, numbins), nu0(0:numdt, numbins), gg(numbins), &
      ambs(numdf, numdt), jmbs(numdf, numdt), jnut(numdf, numdt), &
      jssc(numdf, numdt), anut(numdf, numdt), jeic(numdf, numdt), &
      Qinj(numdt, numbins), Ddif(numdt, numbins))

   dfreqs = 1d0
   dtimes = 1d0

   !
   !   # #    # # #####     ####   ####  #    # #####
   !   # ##   # #   #      #    # #    # ##   # #    #
   !   # # #  # #   #      #      #    # # #  # #    #
   !   # #  # # #   #      #      #    # #  # # #    #
   !   # #   ## #   #      #    # #    # #   ## #    #
   !   # #    # #   #       ####   ####  #    # #####
   !
   nu_ext = 1e14 * gamma_bulk
   u_ext = 1e-4 * gamma_bulk**2
   uB = B**2 / (8d0 * pi)
   nu0_B = 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)
   if ( with_cool .and. (cool_kind == 3) ) then
      nu0 = nu0_B * (Rinit / R0)**(-b_index)
   else
      nu0 = nu0_B
   end if
   
#ifdef HYB
   zetae = 0.99d0
#else
   zetae = 1.0d0
#endif
   Qnth = zetae * Q0
   Qth = (1d0 - zetae) * Q0
   Ddif = 1d-200
   tesc = 1.5 * R / cLight

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
   aux_zero_arr = 0d0

   write(*, *) '---> Calculating the emission'
   write(*, *) ''
   write(*, "(' Using tstep = ', F5.3)") tstep
   write(*, "(' Acceleration period = ',ES15.7)") dtacc
   write(*, "(' Initial synchrotron cooling time scale:', ES15.7)") 1d0 / (nu0_B * g2)
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
         t(i) = ((tmax - tstep)) * (dble(i) / dble(numdt - 1))
      case default
         write(*,*) "Wrong time-grid selection"
      end select
      dt(i) = t(i) - t(i - 1)
      dtimes(i) = 1d0 / dt(i)

      !  ###### ###### #####
      !  #      #      #    #
      !  #####  #####  #    #
      !  #      #      #    #
      !  #      #      #    #
      !  ###### ###### #####
      Qinj(i, :) = injection(t(i), dtacc, gg, g1, g2, qind, theta_e, Qth, Qnth)
      ! call FP_FinDif_difu(dt(i), gg, nn(i - 1, :), nn(i, :), nu0(i - 1, :), Ddif(i, :), Qinj(i, :), tesc)
      call FP_FinDif_cool(dt(i), gg, nn(i - 1, :), nn(i, :), nu0(i - 1, :), Qinj(i, :), tesc)
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

      if ( with_abs ) then
         ! ----->>   Absorption   <<-----
         call mbs_absorption(ambs(:, i), freqs, gg, nn(i, :), B)
         ! ----->>   Entering self-Compton   <<-----
         if ( with_ssc ) then
            call RadTrans_blob(Imbs, R, jmbs(:, i), ambs(:, i))
            call SSC_pwlEED(jssc(:, i), freqs, Imbs, nn(i, :), gg)
            call EIC_pwlEED(jeic(:, i), freqs, u_ext, nu_ext, nn(i, :), gg)
            !!! NOTE: The radius R is used following Chiaberge & Ghisellini (2009)
            ! call RadTrans(Imbs(:, i), R, jnu=jmbs(:, i), anu=ambs(:, i))
            ! call ssc_emissivity(freqs, gg, nn(i, :), Imbs, jssc(:, i))
         else
            jeic(:, i) = 0d0
            jssc(:, i) = 0d0
         end if
         ! ----->>   Totals   <<-----
         jnut(:, i) = jmbs(:, i) + jssc(:, i) + jeic(:, i)
         anut(:, i) = ambs(:, i)
      else
         ! ----->>   Absorption   <<-----
         ambs(:, i) = 0d0
         ! ----->>   Entering self-Compton   <<-----
         if ( with_ssc ) then
            call RadTrans_blob(Imbs, R, jmbs(:, i), ambs(:, i))
            call SSC_pwlEED(jssc(:, i), freqs, Imbs, nn(i, :), gg)
            call EIC_pwlEED(jeic(:, i), freqs, u_ext, nu_ext, nn(i, :), gg)
            !!! NOTE: The radius R is used following Chiaberge & Ghisellini (2009)
            ! call RadTrans(Imbs(:, i), R, jnu=jmbs(:, i))
            ! call ssc_emissivity(freqs, gg, nn(i, :), Imbs, jssc(:, i))
         else
            jeic(:, i) = 0d0
            jssc(:, i) = 0d0
         end if
         ! ----->>   Totals   <<-----
         jnut(:, i) = jmbs(:, i) + jssc(:, i) + jeic(:, i)
         anut(:, i) = ambs(:, i)
      end if


      !   ####   ####   ####  #      # #    #  ####
      !  #    # #    # #    # #      # ##   # #    #
      !  #      #    # #    # #      # # #  # #
      !  #      #    # #    # #      # #  # # #  ###
      !  #    # #    # #    # #      # #   ## #    #
      !   ####   ####   ####  ###### # #    #  ####
      if ( i < numdt ) then
         if ( with_cool ) then
            select case (cool_kind)
            case (1)
               !!!   ---{   Numerical radiative cooling   }---
               nu0(i, :) = IC_cool(nu0_B, gg, freqs, R * jssc(:, i))
            case (2)
               !!!   ---{   Zacharias & Schlickeiser, 2012, 761, 110   }---
               if (i == 1) write(*,*) "Using ZS12 cooling"
               nu0(i, :) = ssccZS12(gg, nn(i, :), B, R)
            case (3)
               !!!   ---{   Uhm & Zhang, 2014, Nat. Phys., 10, 351   }---
               if (i == 1) write(*, *) "Using power-law dacaying magnetic field"
               nu0(i, :) = nu0_B * ((Rinit + cLight * gamma_bulk * t(i)) / R0)**(-2d0 * b_index)
            case default
               write(*, *) 'Wrong cooling selection'
            end select
         else
            nu0(i, :) = nu0_B
         end if
      end if

      Ntot(i) = sum(nn(i, :) * dg, mask=(nn(i, :) > 1d-100))

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
   call h5io_wint0(group_id, 'cooling-kind', cool_kind, herror)
   call h5io_wint0(group_id, 'time-grid', time_grid, herror)
   call h5io_wdble0(group_id, 't_max', tmax, herror)
   call h5io_wdble0(group_id, 'tstep', tstep, herror)
   call h5io_wdble0(group_id, 'Dt_inj', dtacc, herror)
   call h5io_wdble0(group_id, 'Bfield', B, herror)
   call h5io_wdble0(group_id, 'R', R, herror)
   call h5io_wdble0(group_id, 'R0', R0, herror)
   call h5io_wdble0(group_id, 'Rinit', Rinit, herror)
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
   call h5io_wdble0(group_id, 'Q0', Q0, herror)
   call h5io_wdble0(group_id, 'Theta_e', theta_e, herror)
   call h5io_wdble0(group_id, 'zeta_e', zetae, herror)
   call h5io_wdble0(group_id, 'nu_min', numin, herror)
   call h5io_wdble0(group_id, 'nu_max', numax, herror)
   call h5io_closeg(group_id, herror)

   ! ------  Saving data  ------
   call h5io_wint0(file_id, 't_stop', tstop, herror)

   call h5io_wdble0(file_id, 'Q_th', Qth, herror)
   call h5io_wdble0(file_id, 'Q_nth', Qnth, herror)
   call h5io_wdble0(file_id, 'u_B', uB, herror)
   call h5io_wdble0(file_id, 'nu0_B', nu0_B, herror)

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
