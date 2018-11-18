subroutine Paramo(params_file, output_file, with_cool, with_abs, with_ssc, mbs_or_syn)
   use data_types
   use constants
   use misc
   use pwl_integ
   use hdf5
   use h5_inout
   use SRtoolkit
   use anaFormulae
#ifdef MBS
   use magnetobrem
#endif
   use radiation
   implicit none

   character(len=*), intent(in) :: output_file, params_file
   logical, intent(in) :: with_cool, with_abs, with_ssc, mbs_or_syn

   integer, parameter :: nmod=100
   character(len=*), parameter :: on_screen = &
      "(' Time iteration:', I6, ' | Time =', ES15.7, ' | Time step =', ES15.7, ' | nu_0(g_max) =', ES15.7, ' | Ntot =', ES15.7)"

   integer :: i, j, k, kg2, numbins, numdf, numdt, ii, ios
   integer :: time_grid, cool_kind, i_start, i_edge, herror
   integer(HID_T) :: file_id, group_id

   real(dp) :: uB, R, Q0, gmin, gmax, numin, numax, qind, B, dtacc, g1, &
      g2, nu0_B, tstep, zetae, Qth, Qnth, theta_e, tmax, d_lum, z, sind, D, &
      gamma_bulk, theta_obs, R0, Rinit, b_index, mu_obs, mu_com, tob_max, tob_min
   real(dp), allocatable, dimension(:) :: freqs, gg, t, dg, Ntot, &
      aux_zero_arr, sen_lum, dfreqs, dtimes, dt, nu_obs, t_obs, to_com
   real(dp), allocatable, dimension(:, :) :: nu0, nn, Qinj, jnut, jmbs, &
      jssc, ambs,anut,Imbs,Inut,Issc,Iobs

   logical :: at_the_edge

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

   allocate(t(0:numdt), gg(numbins), freqs(numdf), dg(numbins), Ntot(numdt),&
      aux_zero_arr(numbins - 1), dfreqs(numdf), dtimes(numdt), &
      sen_lum(numdt), dt(numdt), nu_obs(numdf), t_obs(numdt), to_com(numdt))
   allocate(nn(numdt,numbins), Qinj(numdt,numbins), nu0(numdt,numbins), &
      ambs(numdf,numdt), jmbs(numdf,numdt), jnut(numdf,numdt), &
      jssc(numdf,numdt), Imbs(numdf,numdt), Inut(numdf,numdt), &
      Issc(numdf,numdt), anut(numdf,numdt), Iobs(numdf,numdt))

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
   uB = B**2 / (8d0 * pi)
   nu0_B = 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)
   if ( with_cool .and. (cool_kind == 3) ) then
      nu0 = nu0_B * (Rinit / R0)**(-b_index)
   else
      nu0 = nu0_B
   end if
   
#ifdef HYB
   zetae = 0.75d0
#else
   zetae = 1.0d0
#endif
   Qnth = zetae * Q0
   Qth = (1d0 - zetae) * Q0
   print*, Qnth, Qth

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

   ! build_t: do i=1, numdt
   !    t_obs = tstep * ( (tmax / tstep)**(dble(i - 2) / dble(numdt - 2)) )
   !    ! t(i) = dble(i) * tmax / dble(numdt)
   ! end do build_t
   ! do i=1,numdt
   !    t_obs(i) = t(1) * ( t(numdt) / t(1))**(dble(i - 1) / dble(numdt - 1))
   ! end do

   t(0) = 0d0
   nn = 1d-200
   jnut = 0d0
   Inut = 0d0
   aux_zero_arr = 0d0
   kg2 = locate(gg, g2, .true.)

   write(*, *) '---> Calculating the emission'
   write(*, *) ''
   write(*, "(' Using tstep = ', F5.3)") tstep
   write(*, "(' Acceleration period = ',ES12.6)") dtacc
   write(*, "(' Initial synchrotron cooling time scale:', ES15.7)") 1d0 / (nu0_B * g2)
   write(*, *) 'Wrting data in: ', trim(output_file)
   write(*, *) ''


   !
   ! ###### #    #  ####  #      #    # ##### #  ####  #    #
   ! #      #    # #    # #      #    #   #   # #    # ##   #
   ! #####  #    # #    # #      #    #   #   # #    # # #  #
   ! #      #    # #    # #      #    #   #   # #    # #  # #
   ! #       #  #  #    # #      #    #   #   # #    # #   ##
   ! ######   ##    ####  ######  ####    #   #  ####  #    #
   !
   time_loop: do i=1, numdt

      select case(time_grid)
      case(1)
         t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numdt - 1)) )
      case(2)
         t(i) = t(i - 1) + tstep / (nu0(max0(1, i - 1), kg2) * g2)
      case(3)
         t(i) = ((tmax - tstep)) * (dble(i) / dble(numdt - 1))
      case default
         write(*,*) "Wrong time-grid selection"
      end select
      dt(i) = t(i) - t(i - 1)
      dtimes(i) = 1d0 / dt(i)
      Qinj(i, :) = injection(t(i), dtacc, gg, g1, g2, qind, theta_e, Qth, Qnth)


      !  ###### ###### #####
      !  #      #      #    #
      !  #####  #####  #    #
      !  #      #      #    #
      !  #      #      #    #
      !  ###### ###### #####
      ! call tridag_ser(aux_zero_arr, &
      !    1d0 + nu0(i - 1,:) * gg**2 * dt / dg, &
      !    - nu0(i - 1,2:) * gg(2:)**2 * dt / dg(:numbins - 1), &
      !    nn(i - 1,:) + dt * Qinj(i - 1,:), &
      !    nn(i,:))
      call cooling_lines(nn(i, :), gg, nu0(:i, :), t(1:i), dtacc, Qinj(:i, :))


      !   ----->   Then we compute the light path
      sen_lum(i) = sum(dt(:i)) * cLight * mu_com
      ! sen_lum(i) = dmin1( sum(dt(:i)) * cLight, 2d0 * R )* mu_com

      !  #####    ##   #####  #   ##   ##### #  ####  #    #
      !  #    #  #  #  #    # #  #  #    #   # #    # ##   #
      !  #    # #    # #    # # #    #   #   # #    # # #  #
      !  #####  ###### #    # # ######   #   # #    # #  # #
      !  #   #  #    # #    # # #    #   #   # #    # #   ##
      !  #    # #    # #####  # #    #   #   #  ####  #    #

      !!!*   First we compute the emissivity
      jmbs(:, i) = mbs_emissivity(freqs, gg, nn(i, :), B, mbs_or_syn)

      !!!*   We first compute (or not) the MBS absorption
      if ( with_abs ) then

         ambs(:, i) = mbs_absorption(freqs, gg, nn(i, :), B, mbs_or_syn)
         anut(:, i) = ambs(:, i)

         !!!*   Then we compute the MBS radiation field
         call solRadTrans_pathIndep(Imbs(:, i), dt(i) * cLight, jnu=jmbs(:, i), anu=ambs(:, i))

         !!!*   Then we compute (or not) the MBS self-Compton
         if ( with_ssc ) then
            !!!*   We compute the MSC locally
            call ssc_emissivity(freqs, gmin, g2, gg, nn(i, :), Imbs(:, i), jssc(:, i))
            !!!*   Then we compute the MSC radiation field
            call solRadTrans_pathIndep(Issc(:, i), dt(i) * cLight, jnu=jssc(:, i), anu=anut(:, i))
         else
            jssc(:, i) = 0d0
            Issc(:, i) = 0d0
         end if

         !!!*   Then we compute the total emissivity
         jnut(:, i) = jmbs(:, i) + jssc(:, i)

         !!!*   Then we solve the radiative transfer equation for the total intensity
         call solRadTrans_pathIndep(Inut(:, i), dt(i) * cLight, jnu=jnut(:, i), anu=anut(:, i))

      else

         ambs(:, i) = 0d0
         anut(:, i) = 0d0

         !!!*   Then we compute the MBS radiation field
         call solRadTrans_pathIndep(Imbs(:, i), dt(i) * cLight, jnu=jmbs(:, i))

         !!!*   Then we compute (or not) the MBS self-Compton
         if ( with_ssc ) then
            !!!*   We compute the MSC locally
            call ssc_emissivity(freqs, gmin, g2, gg, nn(i,: ), Imbs(:, i), jssc(:, i))
            !!!*   Then we compute the MSC intensity locally
            call solRadTrans_pathIndep(Issc(:, i), dt(i) * cLight, jnu=jssc(:, i))
         else
            jssc(:, i) = 0d0
            Issc(:, i) = 0d0
         end if

         !!!*   Then we compute the total emissivity
         jnut(:, i) = jmbs(:, i) + jssc(:, i)

         !!!*   Then we solve the radiative transfer equation for the total intensity
         call solRadTrans_pathIndep(Inut(:, i), dt(i) * cLight, jnu=jnut(:, i))

      end if


      !!!*   With or without cooling
      if ( i < numdt ) then
         if ( with_cool ) then
            select case (cool_kind)
            case (1)
               !!!   ---{   Numerical radiative cooling
               nu0(i, :) = ssc_cool_coef(nu0_B, gg, freqs, Inut(:,i))
               !!!   }---
            case (2)
               !!!   ---{   Zacharias & Schlickeiser, 2012, 761, 110
               if (i == 1) write(*,*) "Using ZS12 cooling"
               nu0(i, :) = ssccZS12(gg, nn(i,:), B, R)
               !!!   }---
            case (3)
               !!!   ---{   Uhm & Zhang, 2014, Nat. Phys., 10, 351
               if (i == 1) write(*,*) "Using power-law dacaying magnetic field"
               nu0(i, :) = nu0_B * ((Rinit + cLight * gamma_bulk * t(i)) / R0)**(-2d0 * b_index)
               !!!   }---
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


   !
   !  ####    ##   #    # # #    #  ####
   ! #       #  #  #    # # ##   # #    #
   !  ####  #    # #    # # # #  # #
   !      # ###### #    # # #  # # #  ###
   ! #    # #    #  #  #  # #   ## #    #
   !  ####  #    #   ##   # #    #  ####
   !
   write(*, *) ""
   ! call system("[ -f "//trim(output_file)//" ] && rm "//trim(output_file)//" && echo 'output: '"//trim(output_file)//" || echo 'output: '"//trim(output_file))

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
   call h5io_wdble0(file_id, 'Q_th', Qth, herror)
   call h5io_wdble0(file_id, 'Q_nth', Qnth, herror)
   call h5io_wdble0(file_id, 'u_B', uB, herror)
   call h5io_wdble0(file_id, 'nu0_B', nu0_B, herror)
   call h5io_wdble1(file_id, 'time', t(1:), herror)
   call h5io_wdble1(file_id, 't_obs', t_obs, herror)
   call h5io_wdble1(file_id, 'Ntot', Ntot, herror)
   call h5io_wdble1(file_id, 'sen_lum', sen_lum, herror)
   call h5io_wdble1(file_id, 'frequency', freqs, herror)
   call h5io_wdble1(file_id, 'nu_obs', nu_obs, herror)
   call h5io_wdble1(file_id, 'gamma', gg, herror)
   call h5io_wdble2(file_id, 'jnut', jnut, herror)
   ! call h5io_wdble2(file_id, 'jmbs', jmbs, herror)
   ! call h5io_wdble2(file_id, 'jssc', jssc, herror)
   ! call h5io_wdble2(file_id, 'ambs', ambs, herror)
   call h5io_wdble2(file_id, 'anut', anut, herror)
   ! call h5io_wdble2(file_id, 'Imbs', Imbs, herror)
   ! call h5io_wdble2(file_id, 'Issc', Issc, herror)
   call h5io_wdble2(file_id, 'Inut', Inut, herror)
   call h5io_wdble2(file_id, 'Iobs', Iobs, herror)
   call h5io_wdble2(file_id, 'Qinj', Qinj, herror)
   call h5io_wdble2(file_id, 'distrib', nn, herror)
   call h5io_wdble2(file_id, 'nu0_tot', nu0, herror)

   ! ------  Closing output file  ------
   call h5io_closef(file_id, herror)
   call h5close_f(herror)

   write(*,*) '=======  FINISHED  ======='
   write(*,*) ''

end subroutine Paramo
