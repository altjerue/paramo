subroutine interJets(params_file, output_file)
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
   ! logical, intent(in) :: cool_withKN, with_abs
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
   integer :: i, j, k, NG, NT, NF
   real(dp) :: uB, uext, R, gmin, gmax, numin, numax, pind, B, D, g1, g2, &
         tstep, Qnth, tmax, d_lum, z, tinj, gamma_bulk, theta_obs, &
         mu_obs, nu_ext, tesc, tlc, L_jet, volume, sigma, beta_bulk, L_B, &
         eps_B, f_rec, urad_const, f_esc, eps_acc, L_e
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
   gmin = par_gmin
   gmax = par_gmax
   pind = par_pind
   numin = par_numin
   numax = par_numax
   NG = par_NG
   NT = par_NT
   NF = par_NF
   f_esc = par_fesc
   with_cool = .false.

   allocate(t(0:NT),freqs(NF),Ntot(NT),Inu(NF),dt(NT),nu_obs(NF),t_obs(NT),&
         dg(NG),urad(NG))
   allocate(nn(NG,0:NT),dotg(NG,0:NT),gg(NG),ambs(NF,NT),jmbs(NF,NT),&
         jnut(NF,NT),jssc(NF,NT),anut(NF,NT),jeic(NF,NT),Qinj(NG,NT),&
         Ddiff(NG,0:NT))

   !   # #    # # #####     ####   ####  #    # #####
   !   # ##   # #   #      #    # #    # ##   # #    #
   !   # # #  # #   #      #      #    # # #  # #    #
   !   # #  # # #   #      #      #    # #  # # #    #
   !   # #   ## #   #      #    # #    # #   ## #    #
   !   # #    # #   #       ####   ####  #    # #####

   !-----> Blob properties
   !NOTE: We need the bulk properties of the blob in both phases. This in order to transform the properties of the fluid into its reference frame
   beta_bulk = bofg(gamma_bulk)           ! Beta bulk
   R = par_R                              ! Radius of the blob

   !----->    Magnetic field
   L_B = sigma * L_jet / (1d0 + sigma)
   uB = L_B / (2d0 * pi * R**2 * cLight * beta_bulk * gamma_bulk * (gamma_bulk - 1d0))
   B = dsqrt(uB * 8d0 * pi)

   !----->   Injection of particles
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
   !    g2 = dsqrt(6d0 * pi * eCharge * 1d-3 / (sigmaT * B))
   end if

   volume = 4d0 * pi * R**3 / 3d0
   tlc = R / cLight
   tesc = f_esc * tlc
   tinj = tlc

   L_e = f_rec * L_B / (1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0))
   ! Qnth = f_rec * uB * pwl_norm(tlc * energy_e, pind - 1d0, g1, g2)
   Qnth = f_rec * L_B * pwl_norm(1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0) * volume * energy_e, pind - 1d0, g1, g2)
   ! Qnth = f_rec * L_B * pwl_norm(volume * energy_e, pind - 1d0, g1, g2)
   ! Qnth = f_rec * L_B / ( (g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6)) * 1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0) * volume * energy_e )

   ! ----->   Output with initial setup
   write(*, "('--> Simulation setup')")
   write(*, "('theta_obs =', ES15.7)") par_theta_obs
   write(*, "('Doppler   =', ES15.7)") D
   write(*, "('gamma_1   =', ES15.7)") g1
   write(*, "('gamma_2   =', ES15.7)") g2
   write(*, "('L_jet     =', ES15.7)") L_jet
   write(*, "('L_B       =', ES15.7)") L_B
   write(*, "('L_e       =', ES15.7)") L_e
   write(*, "('Q_nth     =', ES15.7)") Qnth
   write(*, "('u_B       =', ES15.7)") uB
   write(*, "('B         =', ES15.7)") B
   write(*, "('u_ext     =', ES15.7)") par_uext
   write(*, "('nu_ext    =', ES15.7)") par_nu_ext
   write(*, "('sigma     =', ES15.7)") sigma
   write(*, "('Gamma     =', ES15.7)") gamma_bulk
   write(*, "('t_dyn     =', ES15.7)") tlc
   write(*, "('t_esc     =', ES15.7)") tesc
   write(*, "('t_inj     =', ES15.7)") tinj
   write(*, "('R_b       =', ES15.7)") R

   build_f: do j=1, NF
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(NF - 1)) )
      freqs(j) = nu_com_f(nu_obs(j), z, D)
   end do build_f

   build_g: do k = 1, NG
      gg(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(NG - 1))
      ! gg(k) = (gmin - 1d0) * ((gmax - 1d0) / (gmin - 1d0))**(dble(k - 1) / dble(NG - 1)) + 1d0
      if ( k > 1 ) dg(k) = gg(k) - gg(k - 1)
   end do build_g
   dg(1) = dg(2)

   t(0) = 0d0
   dotg(:, 0) = urad_const * (uB + uext) * pofg(gg)**2
   Ddiff(:, 0) = 1d-200

   Qinj(:, 0) = injection_pwl(0d0, tinj, gg, g1, g2, pind, Qnth)
   nn(:, 0) = Qinj(:, 0)

   write(*, "('--> Calculating the emission')")
   write(*, *) ''
   write(*, "('Using tstep = ', F5.3)") tstep
   write(*, "('Wrting data in: ', A)") trim(output_file)
   write(*, *) ''
   write(*, *) screan_head

   call plasmoid_phaseI  !TODO: arguments
   call plasmoid_phaseII !TODO: arguments


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
   call h5io_wint0(group_id, 'NT', NT, herror)
   call h5io_wint0(group_id, 'NF', NF, herror)
   call h5io_wint0(group_id, 'NG', NG, herror)
   call h5io_wdble0(group_id, 't_max', tmax, herror)
   call h5io_wdble0(group_id, 'tstep', tstep, herror)
   call h5io_wdble0(group_id, 'R_b', R, herror)
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
   call h5io_closeg(group_id, herror)
   !----->   Saving data Phase I
   call h5io_createg(file_id, "Phase_I", group_id, herror)
   call h5io_closeg(group_id, herror)
   !----->   Saving data Phase II
   call h5io_createg(file_id, "Phase_II", group_id, herror)
   call h5io_wdble0(file_id, 't_inj', tinj, herror)
   call h5io_wdble0(file_id, 't_esc', tesc, herror)
   call h5io_wdble0(file_id, 'Bfield', B, herror)
   call h5io_wdble0(file_id, 'L_B', L_B, herror)          ! Eq. (8)
   call h5io_wdble0(file_id, 'uB', uB, herror)            ! Eq. (9)
   call h5io_wdble0(file_id, 'L_e', L_e, herror)          ! Eq. (10)
   call h5io_wdble0(file_id, 'Q_nth', Qnth, herror)       ! Eq. (13)
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
   call h5io_closeg(group_id, herror)
   !----->   Closing output file
   call h5io_closef(group_id, herror)
   call h5close_f(herror)
#endif
   write(*,*) '=======  FINISHED  ======='
   write(*,*) ''

end subroutine interJets




! ######                                 ###
! #     # #    #   ##    ####  ######     #
! #     # #    #  #  #  #      #          #
! ######  ###### #    #  ####  #####      #
! #       #    # ######      # #          #
! #       #    # #    # #    # #          #
! #       #    # #    #  ####  ######    ###
!
!TODO: Phse I: solve kinetic equation for dynamic cooling and with injection of particles
subroutine plasmoid_phaseI !(params_file, output_file)
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

   ! logical, intent(in) :: cool_withKN, with_abs
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
   integer :: i, j, k, NG, NT, NF
   real(dp) :: uB, uext, R, gmin, gmax, numin, numax, pind, B, D, g1, g2, &
         tstep, Qnth, tmax, d_lum, z, tinj, gamma_bulk, theta_obs, &
         mu_obs, nu_ext, tesc, tlc, L_jet, volume, sigma, beta_bulk, L_B, &
         eps_B, f_rec, urad_const, f_esc, eps_acc, L_e
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
   gmin = par_gmin
   gmax = par_gmax
   pind = par_pind
   numin = par_numin
   numax = par_numax
   NG = par_NG
   NT = par_NT
   NF = par_NF
   f_esc = par_fesc
   with_cool = .false.

   allocate(t(0:NT),freqs(NF),Ntot(NT),Inu(NF),dt(NT),nu_obs(NF),t_obs(NT),&
         dg(NG),urad(NG))
   allocate(nn(NG,0:NT),dotg(NG,0:NT),gg(NG),ambs(NF,NT),jmbs(NF,NT),&
         jnut(NF,NT),jssc(NF,NT),anut(NF,NT),jeic(NF,NT),Qinj(NG,NT),&
         Ddiff(NG,0:NT))

   !   # #    # # #####     ####   ####  #    # #####
   !   # ##   # #   #      #    # #    # ##   # #    #
   !   # # #  # #   #      #      #    # # #  # #    #
   !   # #  # # #   #      #      #    # #  # # #    #
   !   # #   ## #   #      #    # #    # #   ## #    #
   !   # #    # #   #       ####   ####  #    # #####

   !-----> Blob properties
   !NOTE: We need the bulk properties of the blob in both phases. This in order to transform the properties of the fluid into its reference frame
   beta_bulk = bofg(gamma_bulk)           ! Beta bulk
   R = par_R                              ! Radius of the blob

   !----->    Magnetic field
   L_B = sigma * L_jet / (1d0 + sigma)
   uB = L_B / (2d0 * pi * R**2 * cLight * beta_bulk * gamma_bulk * (gamma_bulk - 1d0))
   B = dsqrt(uB * 8d0 * pi)

   !----->   Injection of particles
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
      ! g2 = dsqrt(6d0 * pi * eCharge * 1d-3 / (sigmaT * B))
   end if

   volume = 4d0 * pi * R**3 / 3d0
   tlc = R / cLight
   tesc = f_esc * tlc
   tinj = tlc

   L_e = f_rec * L_B / (1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0))
!   Qnth = f_rec * uB * pwl_norm(tlc * energy_e, pind - 1d0, g1, g2)
   Qnth = f_rec * L_B * pwl_norm(1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0) * volume * energy_e, pind - 1d0, g1, g2)
   ! Qnth = f_rec * L_B * pwl_norm(volume * energy_e, pind - 1d0, g1, g2)
   ! Qnth = f_rec * L_B / ( (g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6)) * 1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0) * volume * energy_e )

   ! ----->   Output with initial setup
   write(*, "('--> Simulation setup')")
   write(*, "('theta_obs =', ES15.7)") par_theta_obs
   write(*, "('Doppler   =', ES15.7)") D
   write(*, "('gamma_1   =', ES15.7)") g1
   write(*, "('gamma_2   =', ES15.7)") g2
   write(*, "('L_jet     =', ES15.7)") L_jet
   write(*, "('L_B       =', ES15.7)") L_B
   write(*, "('L_e       =', ES15.7)") L_e
   write(*, "('Q_nth     =', ES15.7)") Qnth
   write(*, "('u_B       =', ES15.7)") uB
   write(*, "('B         =', ES15.7)") B
   write(*, "('u_ext     =', ES15.7)") par_uext
   write(*, "('nu_ext    =', ES15.7)") par_nu_ext
   write(*, "('sigma     =', ES15.7)") sigma
   write(*, "('Gamma     =', ES15.7)") gamma_bulk
   write(*, "('t_dyn     =', ES15.7)") tlc
   write(*, "('t_esc     =', ES15.7)") tesc
   write(*, "('t_inj     =', ES15.7)") tinj
   write(*, "('R_b       =', ES15.7)") R

   build_f: do j=1, NF
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(NF - 1)) )
      freqs(j) = nu_com_f(nu_obs(j), z, D)
   end do build_f

   build_g: do k = 1, NG
      gg(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(NG - 1))
      ! gg(k) = (gmin - 1d0) * ((gmax - 1d0) / (gmin - 1d0))**(dble(k - 1) / dble(NG - 1)) + 1d0
      if ( k > 1 ) dg(k) = gg(k) - gg(k - 1)
   end do build_g
   dg(1) = dg(2)

   t(0) = 0d0
   dotg(:, 0) = urad_const * (uB + uext) * pofg(gg)**2
   Ddiff(:, 0) = 1d-200

   Qinj(:, 0) = injection_pwl(0d0, tinj, gg, g1, g2, pind, Qnth)
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
   time_loop: do i = 1, NT

      t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(NT - 1)) )
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
      !$omp parallel do collapse(1) schedule(auto) default(shared) &
      !$omp& private(j)
      do j = 1, NF
         call syn_emissivity(jmbs(j, i), freqs(j), gg, nn(:, i), B)
         call syn_absorption(ambs(j, i), freqs(j), gg, nn(:, i), B)
      end do
      !$omp end parallel do

      call RadTrans_blob(Inu, R, jmbs(:, i), ambs(:, i))

      !$omp parallel do collapse(1) schedule(auto) default(shared) &
      !$omp& private(j)
      do j = 1, NF
         call IC_iso_powlaw(jssc(j, i), freqs(j), freqs, Inu, nn(:, i), gg)
         jnut(j, i) = jmbs(j, i) + jssc(j, i)
         anut(j, i) = ambs(j, i)
      end do
      !$omp end parallel do


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
         call rad_cool_pwl(dotg(:, i), gg, freqs, 4d0 * pi * Inu / cLight, .false.)
      end if
      dotg(:, i) = dotg(:, i) + urad_const * (uB + uext) * pofg(gg)**2

   end do time_loop

end subroutine plasmoid_phaseI




! ######                                 ### ###
! #     # #    #   ##    ####  ######     #   #
! #     # #    #  #  #  #      #          #   #
! ######  ###### #    #  ####  #####      #   #
! #       #    # ######      # #          #   #
! #       #    # #    # #    # #          #   #
! #       #    # #    #  ####  ######    ### ###
!
!TODO: Phase II: solve kinetic equation without injection, and synchrotron, SSC and adiabatic cooling
subroutine plasmoid_phaseII!(params_file, output_file)
   use data_types
   use constants
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

   ! logical, intent(in) :: cool_withKN, with_abs
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
   integer :: i, j, k, NG, NT, NF
   real(dp) :: uB, uext, R, gmin, gmax, numin, numax, pind, B, D, g1, g2, &
         tstep, Qnth, tmax, d_lum, z, tinj, gamma_bulk, theta_obs, &
         mu_obs, nu_ext, tesc, tlc, L_jet, volume, sigma, beta_bulk, L_B, &
         eps_B, f_rec, urad_const, f_esc, eps_acc, L_e
   real(dp), allocatable, dimension(:) :: freqs, t, Ntot, Inu, gg, dt, nu_obs, &
         t_obs, dg, urad
   real(dp), allocatable, dimension(:,:) :: dotg, nn, jnut, jmbs, jssc, jeic, &
         ambs, anut, Qinj, Ddiff
   logical :: with_cool

   !TODO: Phse I: solve kinetic equation for dynamic cooling and with injection of particles
   !TODO: Phase II: solve kinetic equation without injection, and synchrotron, SSC and adiabatic cooling

   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   allocate(t(0:NT),freqs(NF),Ntot(NT),Inu(NF),dt(NT),nu_obs(NF),t_obs(NT),&
         dg(NG),urad(NG))
   allocate(nn(NG,0:NT),dotg(NG,0:NT),gg(NG),ambs(NF,NT),jmbs(NF,NT),&
         jnut(NF,NT),jssc(NF,NT),anut(NF,NT),jeic(NF,NT),Qinj(NG,NT),&
         Ddiff(NG,0:NT))

   !   # #    # # #####     ####   ####  #    # #####
   !   # ##   # #   #      #    # #    # ##   # #    #
   !   # # #  # #   #      #      #    # # #  # #    #
   !   # #  # # #   #      #      #    # #  # # #    #
   !   # #   ## #   #      #    # #    # #   ## #    #
   !   # #    # #   #       ####   ####  #    # #####

   !-----> Blob properties
   !NOTE: We need the bulk properties of the blob in both phases. This in order to transform the properties of the fluid into its reference frame
   beta_bulk = bofg(gamma_bulk)           ! Beta bulk
   R = par_R                              ! Radius of the blob

   !----->    Magnetic field
   L_B = sigma * L_jet / (1d0 + sigma)
   uB = L_B / (2d0 * pi * R**2 * cLight * beta_bulk * gamma_bulk * (gamma_bulk - 1d0))
   B = dsqrt(uB * 8d0 * pi)

   !----->   Injection of particles
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
   !    g2 = dsqrt(6d0 * pi * eCharge * 1d-3 / (sigmaT * B))
   end if

   volume = 4d0 * pi * R**3 / 3d0
   tlc = R / cLight
   tesc = f_esc * tlc
   tinj = tlc

   L_e = f_rec * L_B / (1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0))
!   Qnth = f_rec * uB * pwl_norm(tlc * energy_e, pind - 1d0, g1, g2)
   Qnth = f_rec * L_B * pwl_norm(1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0) * volume * energy_e, pind - 1d0, g1, g2)
   ! Qnth = f_rec * L_B * pwl_norm(volume * energy_e, pind - 1d0, g1, g2)
   ! Qnth = f_rec * L_B / ( (g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6)) * 1.5d0 * beta_bulk * gamma_bulk * (gamma_bulk - 1d0) * volume * energy_e )

   ! ----->   Output with initial setup
   write(*, "('--> Simulation setup')")
   write(*, "('theta_obs =', ES15.7)") par_theta_obs
   write(*, "('Doppler   =', ES15.7)") D
   write(*, "('gamma_1   =', ES15.7)") g1
   write(*, "('gamma_2   =', ES15.7)") g2
   write(*, "('L_jet     =', ES15.7)") L_jet
   write(*, "('L_B       =', ES15.7)") L_B
   write(*, "('L_e       =', ES15.7)") L_e
   write(*, "('Q_nth     =', ES15.7)") Qnth
   write(*, "('u_B       =', ES15.7)") uB
   write(*, "('B         =', ES15.7)") B
   write(*, "('u_ext     =', ES15.7)") par_uext
   write(*, "('nu_ext    =', ES15.7)") par_nu_ext
   write(*, "('sigma     =', ES15.7)") sigma
   write(*, "('Gamma     =', ES15.7)") gamma_bulk
   write(*, "('t_dyn     =', ES15.7)") tlc
   write(*, "('t_esc     =', ES15.7)") tesc
   write(*, "('t_inj     =', ES15.7)") tinj
   write(*, "('R_b       =', ES15.7)") R

   build_f: do j=1, NF
      nu_obs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(NF - 1)) )
      freqs(j) = nu_com_f(nu_obs(j), z, D)
   end do build_f

   build_g: do k = 1, NG
      gg(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(NG - 1))
      ! gg(k) = (gmin - 1d0) * ((gmax - 1d0) / (gmin - 1d0))**(dble(k - 1) / dble(NG - 1)) + 1d0
      if ( k > 1 ) dg(k) = gg(k) - gg(k - 1)
   end do build_g
   dg(1) = dg(2)

   t(0) = 0d0
   dotg(:, 0) = urad_const * (uB + uext) * pofg(gg)**2
   Ddiff(:, 0) = 1d-200

   Qinj(:, 0) = injection_pwl(0d0, tinj, gg, g1, g2, pind, Qnth)
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
   time_loop: do i = 1, NT

      t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(NT - 1)) )
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
      !$omp parallel do collapse(1) schedule(auto) default(shared) &
      !$omp& private(j)
      do j = 1, NF
         call syn_emissivity(jmbs(j, i), freqs(j), gg, nn(:, i), B)
         call syn_absorption(ambs(j, i), freqs(j), gg, nn(:, i), B)
      end do
      !$omp end parallel do

      call RadTrans_blob(Inu, R, jmbs(:, i), ambs(:, i))

      !$omp parallel do collapse(1) schedule(auto) default(shared) &
      !$omp& private(j)
      do j = 1, NF
         call IC_iso_powlaw(jssc(j, i), freqs(j), freqs, Inu, nn(:, i), gg)
         jnut(j, i) = jmbs(j, i) + jssc(j, i)
         anut(j, i) = ambs(j, i)
      end do
      !$omp end parallel do


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
         call rad_cool_pwl(dotg(:, i), gg, freqs, 4d0 * pi * Inu / cLight, .false.)
      end if
      dotg(:, i) = dotg(:, i) + urad_const * (uB + uext) * pofg(gg)**2

   end do time_loop


end subroutine plasmoid_phaseII
