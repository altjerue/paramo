subroutine paramo
   use data_types
   use constants
   use misc
   use pwl_integ
   use hdf5
   use SRtoolkit
   use anaFormulae
   use magnetobrem
   use radiation
   implicit none

   ! TODO get rid of all the things that are not mine

   character(len=*), parameter :: on_screen = &
      "(' Time iteration:', I6, ' | Time =', ES15.7, ' | Time step =', ES15.7, ' | nu_0(g_max) =', ES15.7, ' | Ntot =', ES15.7)"

   integer :: i, j, k, kg2, numbins, numdf, numdt, nmod=100, numArgs, ios, ii
   integer :: time_grid, cool_kind, i_start, i_edge
   integer(HID_T) :: file_id, group_id

   real(dp) :: uB, R, Q0, gmin, gmax, numin, numax, qind, B, dtacc, g1, &
      g2, nu0_B, tstep, zetae, Qth, Qnth, theta_e, tmax, d_lum, z, sind, D, &
      gamma_bulk, theta_obs, R0, Rinit, b_index, mu_obs, mu_com, tob_max, tob_min
   real(dp), allocatable, dimension(:) :: freqs, gg, t, dg, Ntot, &
      aux_zero_arr, sen_lum, dfreqs, dtimes, dt, nu_obs, t_obs, to_com
   real(dp), allocatable, dimension(:, :) :: nu0, nn, Qinj, jnut, jmbs, &
      jssc, ambs,anut,Imbs,Inut,Issc,Iobs

   character(len = 256) :: output_file, program_name, wCool, wMBSabs, &
      params_file, file_tail, wSSCem

   logical :: with_cool, with_abs, with_ssc, mbs_or_syn, at_the_edge

   !
   !  ####  #####   ##   #####  ##### #    # #####
   ! #        #    #  #  #    #   #   #    # #    #
   !  ####    #   #    # #    #   #   #    # #    #
   !      #   #   ###### #####    #   #    # #####
   ! #    #   #   #    # #   #    #   #    # #
   !  ####    #   #    # #    #   #    ####  #
   !
   numArgs = command_argument_count()
   call get_command_argument(0, program_name)

   if (numArgs /= 4) call an_error(&
      'Usage: '//new_line('A')//&
      '  '//trim(program_name)//' params_file with-cooling'//new_line('A')//&
      'Options:'//new_line('A')//&
      '  params_file: Parameters file'//new_line('A')//&
      '  with cooling: T/F (True/False)'//new_line('A')//&
      '  with SSC: T/F (True/False)'//new_line('A')//&
      '  with MBS self-absorption: T/F (True/False)')

   !
   !   ::::::::::   Opening parameters file   ::::::::::
   !
   call get_command_argument(1, params_file)
   open(unit=77, file=params_file, iostat=ios)
   if ( ios /= 0 ) call an_error(&
      'Problem opening '//trim(params_file)//'. Usage: '//new_line('A')//&
      '  '//trim(program_name)//' params_file with-cooling'//new_line('A')//&
      'Options:'//new_line('A')//&
      '  params_file: Parameters file'//new_line('A')//&
      '  with cooling: T/F (True/False)'//new_line('A')//&
      '  with SSC: T/F (True/False)'//new_line('A')//&
      '  with MBS self-absorption: T/F (True/False)')
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
   read(77,*) file_tail
   close(77)

   !
   !    -----> With or without cooling
   !
   call get_command_argument(2, wCool)
   if ( wCool == 'T' ) then
      with_cool = .true.
      output_file = 'V'!//trim(dirs(dir_index))
   elseif ( wCool == 'F' ) then
      with_cool = .false.
      output_file = 'C'!//trim(dirs(dir_index))
   else
      call an_error(&
         'Wrong cooling. Value given: '//trim(wCool)//'. Usage: '//new_line('A')//&
         '  '//trim(program_name)//' params_file with-cooling'//new_line('A')//&
         'Options:'//new_line('A')//&
         '  params_file: Parameters file'//new_line('A')//&
         '  with cooling: T/F (True/False)'//new_line('A')//&
         '  with SSC: T/F (True/False)'//new_line('A')//&
         '  with MBS self-absorption: T/F (True/False)')
   end if

   !
   !    -----> With or without absorption
   !
   call get_command_argument(3, wMBSabs)
   if ( wMBSabs == 'T' ) then
      with_abs = .true.
      output_file = trim(output_file)//'O'!//trim(dirs(dir_index))
   elseif ( wMBSabs == 'F' ) then
      with_abs = .false.
      output_file = trim(output_file)//'T'!//trim(dirs(dir_index))
   else
      call an_error(&
         'Wrong absorption. Value given: '//trim(wMBSabs)//'. Usage: '//new_line('A')//&
         '  '//trim(program_name)//' params_file with-cooling'//new_line('A')//&
         'Options:'//new_line('A')//&
         '  params_file: Parameters file'//new_line('A')//&
         '  with cooling: T/F (True/False)'//new_line('A')//&
         '  with SSC: T/F (True/False)'//new_line('A')//&
         '  with MBS self-absorption: T/F (True/False)')
   end if

   !
   !    -----> With or without SSC emissivity
   !
   call get_command_argument(4, wSSCem)
   if ( wSSCem == 'T' ) then
      with_ssc = .true.
      output_file = trim(output_file)//'wSSC'!//trim(dirs(dir_index))
   elseif ( wSSCem == 'F' ) then
      with_ssc = .false.
      output_file = trim(output_file)//'oSSC'!//trim(dirs(dir_index))
   else
      call an_error(&
         'Wrong SSC emissivity. Value given: '//trim(wSSCem)//'. Usage: '//new_line('A')//&
         '  '//trim(program_name)//' params_file with-cooling'//new_line('A')//&
         'Options:'//new_line('A')//&
         '  params_file: Parameters file'//new_line('A')//&
         '  with cooling: T/F (True/False)'//new_line('A')//&
         '  with SSC: T/F (True/False)'//new_line('A')//&
         '  with MBS self-absorption: T/F (True/False)')
   end if


   !
   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   !
   call K2_init
   call RF_init
#ifndef MBS
   call sint_load_table("uinterp.h5")
   mbs_or_syn = .false.
   output_file = 'S'//trim(output_file)
#else
   call sint_load_table("uinterp.h5")
   call load_mb_table("disTable.h5")
   mbs_or_syn = .true.
   output_file = 'M'//trim(output_file)
   print*, globgmax, chunche_c100g20, chunche_c100g100
#endif

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
   nu0_B = 4d0 * sigmaT * uB / (3d0 * me * cspeed)
   if ( with_cool .and. (cool_kind == 3) ) then
      nu0 = nu0_B * (Rinit / R0)**(-b_index)
   else
      nu0 = nu0_B
   end if
   
#ifdef HYB
   zetae = 0.75d0
   output_file = 'H'//trim(output_file)
#else
   zetae = 1.0d0
   output_file = 'P'//trim(output_file)
#endif
   Qnth = zetae * Q0
   Qth = (1d0 - zetae) * Q0
   print*, Qnth, Qth

   ! ----->>   Viewing angle
   mu_obs = dcos(theta_obs * pi / 180d0)
   mu_com = mu_com_f(gamma_bulk, mu_obs)
   D = Doppler(gamma_bulk, mu_obs)

   ! output_file = trim(output_file)//'_FinDif-'//trim(file_tail)//'.h5'
   output_file = trim(output_file)//'-'//trim(file_tail)//'.h5'

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


      !
      !    :::::  Solving matrix  :::::
      !
      ! call tridag_ser(aux_zero_arr, &
      !    1d0 + nu0(i - 1,:) * gg**2 * dt / dg, &
      !    - nu0(i - 1,2:) * gg(2:)**2 * dt / dg(:numbins - 1), &
      !    nn(i - 1,:) + dt * Qinj(i - 1,:), &
      !    nn(i,:))
      call cooling_lines(nn(i, :), gg, nu0(:i, :), t(1:i), dtacc, Qinj(:i, :))


      !!!   ---{   Then we compute the light path
      sen_lum(i) = sum(dt(:i)) * cspeed * mu_com
      ! sen_lum(i) = dmin1( sum(dt(:i)) * cspeed, 2d0 * R )* mu_com
      !!!   }---

      !
      !   ####   ####  #    # #####  #    # #####    ###### #    # #  ####   ####
      !  #    # #    # ##  ## #    # #    #   #      #      ##  ## # #      #
      !  #      #    # # ## # #    # #    #   #      #####  # ## # #  ####   ####
      !  #      #    # #    # #####  #    #   #      #      #    # #      #      #
      !  #    # #    # #    # #      #    #   #      #      #    # # #    # #    #
      !   ####   ####  #    # #       ####    #      ###### #    # #  ####   ####
      !

      !!!*   First we compute the emissivity
      jmbs(:, i) = mbs_emissivity(freqs, gg, nn(i, :), B, mbs_or_syn)

      !!!*   We first compute (or not) the MBS absorption
      if ( with_abs ) then

         ambs(:, i) = mbs_absorption(freqs, gg, nn(i, :), B, mbs_or_syn)
         anut(:, i) = ambs(:, i)

         !!!*   Then we compute the MBS radiation field
         call solRadTrans_pathIndep(Imbs(:, i), dt(i) * cspeed, jnu=jmbs(:, i), anu=ambs(:, i))

         !!!*   Then we compute (or not) the MBS self-Compton
         if ( with_ssc ) then
            !!!*   We compute the MSC locally
            call ssc_emissivity(freqs, gmin, g2, gg, nn(i, :), Imbs(:, i), jssc(:, i))
            !!!*   Then we compute the MSC radiation field
            call solRadTrans_pathIndep(Issc(:, i), dt(i) * cspeed, jnu=jssc(:, i), anu=anut(:, i))
         else
            jssc(:, i) = 0d0
            Issc(:, i) = 0d0
         end if

         !!!*   Then we compute the total emissivity
         jnut(:, i) = jmbs(:, i) + jssc(:, i)

         !!!*   Then we solve the radiative transfer equation for the total intensity
         call solRadTrans_pathIndep(Inut(:, i), dt(i) * cspeed, jnu=jnut(:, i), anu=anut(:, i))

      else

         ambs(:, i) = 0d0
         anut(:, i) = 0d0

         !!!*   Then we compute the MBS radiation field
         call solRadTrans_pathIndep(Imbs(:, i), dt(i) * cspeed, jnu=jmbs(:, i))

         !!!*   Then we compute (or not) the MBS self-Compton
         if ( with_ssc ) then
            !!!*   We compute the MSC locally
            call ssc_emissivity(freqs, gmin, g2, gg, nn(i,: ), Imbs(:, i), jssc(:, i))
            !!!*   Then we compute the MSC intensity locally
            call solRadTrans_pathIndep(Issc(:, i), dt(i) * cspeed, jnu=jssc(:, i))
         else
            jssc(:, i) = 0d0
            Issc(:, i) = 0d0
         end if

         !!!*   Then we compute the total emissivity
         jnut(:, i) = jmbs(:, i) + jssc(:, i)

         !!!*   Then we solve the radiative transfer equation for the total intensity
         call solRadTrans_pathIndep(Inut(:, i), dt(i) * cspeed, jnu=jnut(:, i))

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
               nu0(i, :) = nu0_B * ((Rinit + cspeed * gamma_bulk * t(i)) / R0)**(-2d0 * b_index)
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

   Iobs = cspeed * Iobs

   !
   !  ####    ##   #    # # #    #  ####
   ! #       #  #  #    # # ##   # #    #
   !  ####  #    # #    # # # #  # #
   !      # ###### #    # # #  # # #  ###
   ! #    # #    #  #  #  # #   ## #    #
   !  ####  #    #   ##   # #    #  ####
   !
   write(*, *) ""
   call system("[ -f "//trim(output_file)//" ] && rm "//trim(output_file)//" && echo 'output: '"//trim(output_file)//" || echo 'output: '"//trim(output_file))

   call hm_init
   call hm_create(output_file, file_id)

   ! ------  Saving initial parameters  ------
   call hm_gcreate(file_id, "Params", group_id)

   call hm_write0_int(group_id, numdt,     'numdt')
   call hm_write0_int(group_id, numdf,     'numdf')
   call hm_write0_int(group_id, numbins,   'numbins')
   call hm_write0_int(group_id, cool_kind, 'cooling kind')
   call hm_write0_int(group_id, time_grid, 'time grid')

   call hm_write0_double(group_id, tmax,       't_max')
   call hm_write0_double(group_id, tstep,      'tstep')
   call hm_write0_double(group_id, dtacc,      'Dt_inj')
   call hm_write0_double(group_id, B,          'Bfield')
   call hm_write0_double(group_id, R,          'R')
   call hm_write0_double(group_id, R0,         'R0')
   call hm_write0_double(group_id, Rinit,      'Rinit')
   call hm_write0_double(group_id, d_lum,      'd_lum')
   call hm_write0_double(group_id, z,          'redshift')
   call hm_write0_double(group_id, gamma_bulk, 'gamma_bulk')
   call hm_write0_double(group_id, theta_obs,  'view-angle')
   call hm_write0_double(group_id, b_index,    'Bfield-index')
   call hm_write0_double(group_id, gmin,       'gamma_min')
   call hm_write0_double(group_id, gmax,       'gamma_max')
   call hm_write0_double(group_id, g1,         'gamma_1')
   call hm_write0_double(group_id, g2,         'gamma_2')
   call hm_write0_double(group_id, qind,       'pwl-index')
   call hm_write0_double(group_id, Q0,         'Q0')
   call hm_write0_double(group_id, theta_e,    'Theta_e')
   call hm_write0_double(group_id, zetae,      'zeta_e')
   call hm_write0_double(group_id, numin,      'nu_min')
   call hm_write0_double(group_id, numax,      'nu_max')

   call hm_gclose(group_id)

   ! ------  Saving data  ------
   call hm_write0_double(file_id, Qth,   'Q_th')
   call hm_write0_double(file_id, Qnth,  'Q_nth')
   call hm_write0_double(file_id, uB,    'u_B')
   call hm_write0_double(file_id, nu0_B, 'nu0_B')

   call hm_write1_double(file_id, numdt,   t(1:),   'time')
   call hm_write1_double(file_id, numdt,   t_obs,   't_obs')
   call hm_write1_double(file_id, numdt,   Ntot,    'Ntot')
   call hm_write1_double(file_id, numdt,   sen_lum, 'sen_lum')
   call hm_write1_double(file_id, numdf,   freqs,   'frequency')
   call hm_write1_double(file_id, numdf,   nu_obs,  'nu_obs')
   call hm_write1_double(file_id, numbins, gg,      'gamma')

   call hm_write2_double(file_id, numdf, numdt,   jnut, 'jnut')
   ! call hm_write2_double(file_id, numdf, numdt,   jmbs, 'jmbs')
   ! call hm_write2_double(file_id, numdf, numdt,   jssc, 'jssc')
   ! call hm_write2_double(file_id, numdf, numdt,   ambs, 'ambs')
   call hm_write2_double(file_id, numdf, numdt,   anut, 'anut')
   ! call hm_write2_double(file_id, numdf, numdt,   Imbs, 'Imbs')
   ! call hm_write2_double(file_id, numdf, numdt,   Issc, 'Issc')
   call hm_write2_double(file_id, numdf, numdt,   Inut, 'Inut')
   call hm_write2_double(file_id, numdf, numdt,   Iobs, 'Iobs')
   call hm_write2_double(file_id, numdt, numbins, Qinj, 'Qinj')
   call hm_write2_double(file_id, numdt, numbins, nn,   'distrib')
   call hm_write2_double(file_id, numdt, numbins, nu0,  'nu0_tot')

   call hm_close(file_id)
   call hm_finalize

   write(*,*) '=======  FINISHED  ======='
   write(*,*) ''

end subroutine paramo
