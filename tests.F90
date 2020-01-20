program tests
   use data_types
   use constants
   use params
   use misc
   use pwl_integ
   use hdf5
   use h5_inout
   use SRtoolkit
   use Aglow_models
   use anaFormulae
   use dist_evol
   use radiation
   use K1
   use K2
   implicit none
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! TODO:
   !   [ ] dir as an input parameter
   !   [ ] blast wave test: Numerical vs Sari, Piran & Narayan (1998)
   !   [ ] choise of test as input
   !       [ ] Modify Comala
   !   [ ]
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call steady_state
   ! call rad_procs

   ! write(*, *) '=======  FINISHED  ======='
   write(*, *) ''

contains

   subroutine steady_state
      ! ************************************************************************
      !
      !   This subroutine produces the steady state solutions of the kinetic
      !   equation.
      !
      ! ************************************************************************
      implicit none
      integer(HID_T) :: file_id
      integer :: i, k, numg, numt, herror, numf, j
      real(dp) :: g1, g2, gmin, gmax, tmax, tstep, qind, tacc, tesc, R, B, numax, numin, uB
      real(dp), allocatable, dimension(:) :: t, g, Q0, D0, C0, aux0, dt, dg, &
         Ntot1, Ntot2, Ntot3, Ntot4, Ntot5, Ntot6, zero1, zero2, freqs, &
         Inu1, Inu4, Inu5, Inu6
      real(dp), allocatable, dimension(:, :) :: n1, n2, n3, n4, n5, n6
      real(dp), allocatable, dimension(:, :) :: jmbs1, ambs1, jssc1, jmbs4, ambs4, jssc4, jmbs5, ambs5, jssc5, jmbs6, ambs6, jssc6

      numg = 128
      numt = 300
      numf = 192
      g1 = 1e2
      g2 = 1e6
      gmin = 1.01d0
      gmax = 1.5d0 * g2
      numin = 1d10
      numax = 1d27
      tstep = 1e0
      tmax = 1e7
      qind = 0d0
      R = 1e16
      B = 1d0
      uB = B**2 / (8d0 * pi)

      allocate(g(numg), Q0(numg), D0(numg), t(0:numt), dt(numt), dg(numg), &
         Ntot1(numt), Ntot2(numt), Ntot3(numt), Ntot4(numt), Ntot5(numt), &
         Ntot6(numt), C0(numg), aux0(numg), zero1(numg), zero2(numg), &
         freqs(numf), Inu1(numf), Inu4(numf), Inu5(numf), Inu6(numf))
      allocate(n1(0:numt, numg), n2(0:numt, numg), n3(0:numt, numg), &
         n4(0:numt, numg), n5(0:numt, numg), n6(0:numt, numg), &
         jmbs1(numf, numt), ambs1(numf, numt), jssc1(numf, numt), &
         jmbs4(numf, numt), ambs4(numf, numt), jssc4(numf, numt), &
         jmbs5(numf, numt), ambs5(numf, numt), jssc5(numf, numt), &
         jmbs6(numf, numt), ambs6(numf, numt), jssc6(numf, numt))

      build_f: do j=1,numf
         freqs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numf - 1)) )
      end do build_f

      build_g: do k = 1, numg
         g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numg - 1))
         if ( k > 1 ) dg(k) = g(k) - g(k - 1)
      end do build_g
      dg(1) = dg(2)

      t(0) = 0d0
      zero1 = 1d-200
      zero2 = 0d0
      C0 = 3.48d-11 ! 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)
      tacc = 1d0 / (C0(1) * 10d0**4.5d0) !tesc
      tesc = tacc ! 1d200 ! 1.5d0 * R / cLight !
      D0 = 0.5d0 * pofg(g)**2 / tacc
      n1(0, :) = injection_pwl(1d0, tacc, g, g1, g2, qind, 1d0)
      n2(0, :) = n1(0, :)
      n3(0, :) = n1(0, :)
      n4(0, :) = n1(0, :)
      n5(0, :) = n1(0, :)
      n6(0, :) = n1(0, :)

      write(*, "(4ES15.7)") C0(1), D0(1), tacc, tmax

      time_loop: do i = 1, numt

         t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numt - 1)) )
         dt(i) = t(i) - t(i - 1)

         Q0 = injection_pwl(t(i), tacc, g, g1, g2, qind, 1d0)
         ! call FP_FinDif_cool(dt(i), g, n1(i - 1, :), n1(i, :), C0, zero2, 1d200)
         ! call FP_FinDif_cool(dt(i), g, n2(i - 1, :), n2(i, :), C0, Q0, 1d200)
         ! call FP_FinDif_cool(dt(i), g, n3(i - 1, :), n3(i, :), C0, Q0, tesc)
         call FP_FinDif_difu(dt(i), g, n1(i - 1, :), n1(i, :), C0 * pofg(g)**2, zero1, zero2, 1d200, R / cLight)
         call FP_FinDif_difu(dt(i), g, n2(i - 1, :), n2(i, :), C0 * pofg(g)**2, zero1, Q0,    1d200, R / cLight)
         call FP_FinDif_difu(dt(i), g, n3(i - 1, :), n3(i, :), C0 * pofg(g)**2, zero1, Q0,    tesc,  R / cLight)
         call FP_FinDif_difu(dt(i), g, n4(i - 1, :), n4(i, :), C0 * pofg(g)**2, D0,    zero2, 1d200, R / cLight)
         call FP_FinDif_difu(dt(i), g, n5(i - 1, :), n5(i, :), C0 * pofg(g)**2, D0,    Q0,    1d200, R / cLight)
         call FP_FinDif_difu(dt(i), g, n6(i - 1, :), n6(i, :), C0 * pofg(g)**2, D0,    Q0,    tesc,  R / cLight)

         Ntot1(i) = sum(n1(i, :) * dg)
         Ntot2(i) = sum(n2(i, :) * dg)
         Ntot3(i) = sum(n3(i, :) * dg)
         Ntot4(i) = sum(n4(i, :) * dg)
         Ntot5(i) = sum(n5(i, :) * dg)
         Ntot6(i) = sum(n6(i, :) * dg)

         !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
         !$OMP& PRIVATE(j)
         do j = 1, numf
            call mbs_emissivity(jmbs1(j, i), freqs(j), g, n1(i, :), B)
            call mbs_emissivity(jmbs4(j, i), freqs(j), g, n4(i, :), B)
            call mbs_emissivity(jmbs5(j, i), freqs(j), g, n5(i, :), B)
            call mbs_emissivity(jmbs6(j, i), freqs(j), g, n6(i, :), B)
            call mbs_absorption(ambs1(j, i), freqs(j), g, n1(i, :), B)
            call mbs_absorption(ambs4(j, i), freqs(j), g, n4(i, :), B)
            call mbs_absorption(ambs5(j, i), freqs(j), g, n5(i, :), B)
            call mbs_absorption(ambs6(j, i), freqs(j), g, n6(i, :), B)
         end do
         !$OMP END PARALLEL DO

         call RadTrans_blob(Inu1, R, jmbs4(:, i), ambs1(:, i))
         call RadTrans_blob(Inu4, R, jmbs4(:, i), ambs4(:, i))
         call RadTrans_blob(Inu5, R, jmbs4(:, i), ambs5(:, i))
         call RadTrans_blob(Inu6, R, jmbs6(:, i), ambs6(:, i))

         !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
         !$OMP& PRIVATE(j)
         do j = 1, numf
            call IC_iso_powlaw(jssc1(j, i), freqs(j), freqs, Inu1, n1(i, :), g)
            call IC_iso_powlaw(jssc4(j, i), freqs(j), freqs, Inu4, n4(i, :), g)
            call IC_iso_powlaw(jssc5(j, i), freqs(j), freqs, Inu5, n5(i, :), g)
            call IC_iso_powlaw(jssc6(j, i), freqs(j), freqs, Inu6, n6(i, :), g)
         end do
         !$OMP END PARALLEL DO

      end do time_loop

      call h5open_f(herror)
      call h5io_createf("SSsol_dbg.h5", file_id, herror)
      call h5io_wint0(file_id, 'numdt', numt, herror)
      call h5io_wint0(file_id, 'numbins', numg, herror)
      call h5io_wdble0(file_id, 't_max', tmax, herror)
      call h5io_wdble0(file_id, 'tstep', tstep, herror)
      call h5io_wdble0(file_id, 'radius', R, herror)
      call h5io_wdble0(file_id, 'gamma_min', gmin, herror)
      call h5io_wdble0(file_id, 'gamma_max', gmax, herror)
      call h5io_wdble0(file_id, 'gamma_1', g1, herror)
      call h5io_wdble0(file_id, 'gamma_2', g2, herror)
      call h5io_wdble0(file_id, 'pwl-index', qind, herror)
      call h5io_wdble0(file_id, 't_acc', tacc, herror)
      call h5io_wdble0(file_id, 't_esc', tesc, herror)
      call h5io_wdble1(file_id, 'time', t(1:), herror)
      call h5io_wdble1(file_id, 'gamma', g, herror)
      call h5io_wdble1(file_id, 'freqs', freqs, herror)
      call h5io_wdble2(file_id, 'dist1', n1(1:, :), herror)
      call h5io_wdble2(file_id, 'dist2', n2(1:, :), herror)
      call h5io_wdble2(file_id, 'dist3', n3(1:, :), herror)
      call h5io_wdble2(file_id, 'dist4', n4(1:, :), herror)
      call h5io_wdble2(file_id, 'dist5', n5(1:, :), herror)
      call h5io_wdble2(file_id, 'dist6', n6(1:, :), herror)
      call h5io_wdble2(file_id, 'jsyn1', jmbs1, herror)
      call h5io_wdble2(file_id, 'jsyn4', jmbs4, herror)
      call h5io_wdble2(file_id, 'jsyn5', jmbs5, herror)
      call h5io_wdble2(file_id, 'jsyn6', jmbs6, herror)
      call h5io_wdble2(file_id, 'asyn1', ambs1, herror)
      call h5io_wdble2(file_id, 'asyn4', ambs4, herror)
      call h5io_wdble2(file_id, 'asyn5', ambs5, herror)
      call h5io_wdble2(file_id, 'asyn6', ambs6, herror)
      call h5io_wdble2(file_id, 'jssc1', jssc1, herror)
      call h5io_wdble2(file_id, 'jssc4', jssc4, herror)
      call h5io_wdble2(file_id, 'jssc5', jssc5, herror)
      call h5io_wdble2(file_id, 'jssc6', jssc6, herror)
      call h5io_wdble1(file_id, 'Ntot1', Ntot1, herror)
      call h5io_wdble1(file_id, 'Ntot2', Ntot2, herror)
      call h5io_wdble1(file_id, 'Ntot3', Ntot3, herror)
      call h5io_wdble1(file_id, 'Ntot4', Ntot4, herror)
      call h5io_wdble1(file_id, 'Ntot5', Ntot5, herror)
      call h5io_wdble1(file_id, 'Ntot6', Ntot6, herror)
      call h5io_closef(file_id, herror)
      call h5close_f(herror)

      write(*, "('--> Fokker-Planck solver test')")

   end subroutine steady_state





#if 0
   !  #####    ##   #####     #####  #####   ####  #####
   !  #    #  #  #  #    #    #    # #    # #    # #    #
   !  #    # #    # #    #    #    # #    # #    # #####
   !  #####  ###### #    #    #####  #####  #    # #    #
   !  #   #  #    # #    #    #      #   #  #    # #    #
   !  #    # #    # #####     #      #    #  ####  #####
   subroutine rad_procs
      implicit none
      integer(HID_T) :: file_id
      integer :: numf, numg, j, herror
      real(dp) :: numin, numax, u_ext, u_e, u_0, g1, g2, g3, p, q, alpha, &
         nu1, nu2, nu_ext, B, k0, ke, R
      real(dp), allocatable, dimension(:) :: jmbs, ambs, jeic, jssc, Iin, n, &
         nuj, nua, g

      numf = 256
      numg = 128

      allocate(jmbs(numf), ambs(numf), jeic(numf), jssc(numf), Iin(numf), &
         n(numg), g(numg), nuj(numf), nua(numf))

      numin = 1d5
      numax = 1d30
      B = 50d0
      R = 1e15
      nu_ext = 1d-8 * mass_e * cLight**2 / hPlanck
      u_ext = 1d0
      u_0 = 1d0
      u_e = 1d4
      g1 = 5d0
      g2 = 1d9
      nu1 = 1d-8 * mass_e * cLight**2 / hPlanck
      nu2 = 1d0 * mass_e * cLight**2 / hPlanck
      p = 3.2d0
      alpha = 0.5d0

      ke = u_e * pwl_norm(mass_e * cLight**2, p - 1d0, g1, g2)
      do j = 1, numg
         g(j) = g1 * (g2 / g1)**(dble(j - 1) / dble(numg - 1))
         n(j) = ke * powlaw_dis(g(j), g1, g2, p)
      end do

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      k0 = u_0 * cLight * pwl_norm(1d0, alpha, nu1, nu2)
      do j = 1, numf
         nuj(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numf - 1)) )
         ! Iin(j) = k0 * powlaw_dis(nuj(j), nu1, nu2, alpha)
         call mbs_emissivity(jmbs(j), nuj(j), g, n, B)
         call mbs_absorption(ambs(j), nuj(j), g, n, B)
         call RadTrans(Iin(j), R, jmbs(j), ambs(j))
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numf
         call IC_emis_full(nuj(j), nuj, g, n, Iin, jssc(j))
         call IC_emis_full(nuj(j), nu_ext, g, n, u_ext * cLight / (4d0 * pi), jeic(j))
         ! call EIC_pwlEED(jeic, nuj, u_ext, nu_ext, n, g)
         write(*,*) nuj(j), nuj(j) * jmbs(j), nuj(j) * jssc(j), nuj(j) * jeic(j)
      end do
      !$OMP END PARALLEL DO

      call h5open_f(herror)
      call h5io_createf(dir//"Rad.h5", file_id, herror)
      call h5io_wint0(file_id, 'numdf', numf, herror)
      call h5io_wint0(file_id, 'numg', numg, herror)
      call h5io_wdble0(file_id, 'nu_min', numin, herror)
      call h5io_wdble0(file_id, 'nu_max', numax, herror)
      call h5io_wdble1(file_id, 'nuj', nuj, herror)
      call h5io_wdble1(file_id, 'gamma', g, herror)
      call h5io_wdble1(file_id, 'distrib', n, herror)
      call h5io_wdble1(file_id, 'jmbs', jmbs, herror)
      call h5io_wdble1(file_id, 'jssc', jssc, herror)
      call h5io_wdble1(file_id, 'jeic', jeic, herror)
! #if 0
      !  ----->   Absorption
      B = 1d3
      numin = 1d8
      numax = 1d21
      do j = 1, numf
         nua(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numf - 1)) )
      end do
      call h5io_wdble1(file_id, 'nua', nua, herror)

      g1 = 1d2
      g2 = 1d5
      p = 3.0d0
      ke = u_e * pwl_norm(mass_e * cLight**2, p - 1d0, g1, g2)
      do j = 1, numg
         g(j) = g1 * (g2 / g1)**(dble(j - 1) / dble(numg - 1))
         n(j) = ke * powlaw_dis(g(j), g1, g2, p)
      end do

      call mbs_absorption(ambs, nua, g, n, B)
      call h5io_wdble1(file_id, 'ambs1', ambs, herror)

      g1 = 1.01d0
      g2 = 1d2
      p = 1.3d0
      ke = u_e * pwl_norm(mass_e * cLight**2, p - 1d0, g1, g2)
      do j = 1, numg
         g(j) = g1 * (g2 / g1)**(dble(j - 1) / dble(numg - 1))
         n(j) = ke * powlaw_dis(g(j), g1, g2, p)
      end do

      call mbs_absorption(ambs, nua, g, n, B)
      call h5io_wdble1(file_id, 'ambs2', ambs, herror)

      g1 = 1d1
      g2 = 1d3
      g3 = 1d5
      q = 1.3d0
      p = 3.0d0
      do j = 1, numg
         g(j) = g1 * (g3 / g1)**(dble(j - 1) / dble(numg - 1))
      end do
      ke = u_e * (g2**(q - p) * pwl_norm(mass_e * cLight**2, q - 1d0, g1, g2) + pwl_norm(mass_e * cLight**2, p - 1d0, g2, g3))
      do j = 1, numg
         if ( g(j) < g2 ) then
            n(j) = ke * g2**(q - p) * powlaw_dis(g(j), g1, g2, q)
         else
            n(j) = ke * powlaw_dis(g(j), g2, g3, p)
         end if
      end do

      call mbs_absorption(ambs, nua, g, n, B)
      call h5io_wdble1(file_id, 'ambs3', ambs, herror)
! #endif
      call h5io_closef(file_id, herror)
      call h5close_f(herror)

      ! write(*, "('--> Radiation processes test')")

   end subroutine rad_procs


! #if 0
   subroutine no_diff_CG99

      call K1_init
      call K2_init

      !   # #    # # #####     ####   ####  #    # #####
      !   # ##   # #   #      #    # #    # ##   # #    #
      !   # # #  # #   #      #      #    # # #  # #    #
      !   # #  # # #   #      #      #    # #  # # #    #
      !   # #   ## #   #      #    # #    # #   ## #    #
      !   # #    # #   #       ####   ####  #    # #####
      beta_bulk = bofg(gamma_bulk)

      ! --->    Magnetic field
      L_B = sigma * L_j / (1d0 + sigma)
      uB = L_B / (pi * cLight * beta_bulk * (gamma_bulk * R)**2) ! B**2 / (8d0 * pi)
      B = dsqrt(uB * 8d0 * pi)
      mu_mag = gamma_bulk * (1d0 + sigma)
      nu0 = 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)

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

      do i = 1, numt
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
      end do


   end subroutine no_diff_CG99


   subroutine syn_test
      implicit none
      integer(HID_T) :: file_id
      integer :: numf, numg, j, herror
      real(dp) :: numin, numax, u_ext, u_e, u_0, g1, g2, g3, p, q, alpha, &
         nu1, nu2, nu_ext, B, k0, ke
      real(dp), allocatable, dimension(:) :: jmbs, ambs, jeic, jssc, Iin, n, &
         nuj, nua, g

      numf = 256
      numg = 128

      allocate(jmbs(numf), ambs(numf), jeic(numf), jssc(numf), Iin(numf), &
         n(numg), g(numg), nuj(numf), nua(numf))

      numin = 1d8
      numax = 1d28
      B = 1d0
      nu_ext = 1d-8 * mass_e * cLight**2 / hPlanck
      u_ext = 1d0
      u_0 = 1d0
      u_e = 1d0
      g1 = 1d2
      g2 = 1d7
      nu1 = 1d-9 * mass_e * cLight**2 / hPlanck
      nu2 = 1d-6 * mass_e * cLight**2 / hPlanck
      p = 2.2d0
      alpha = 0.5d0

      ke = u_e * pwl_norm(mass_e * cLight**2, p - 1d0, g1, g2)
      do j = 1, numg
         g(j) = g1 * (g2 / g1)**(dble(j - 1) / dble(numg - 1))
         n(j) = ke * powlaw_dis(g(j), g1, g2, p)
      end do

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numf
         nuj(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numf - 1)) )
         call mbs_emissivity(jmbs(j), nuj(j), g, n, B)
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numf
         call EIC_pwlEED(jeic, nuj, u_ext, nu_ext, n, g)
      end do
      !$OMP END PARALLEL DO

      call h5open_f(herror)
      call h5io_createf(dir//"Syn.h5", file_id, herror)
      call h5io_wint0(file_id, 'numdf', numf, herror)
      call h5io_wint0(file_id, 'numg', numg, herror)
      call h5io_wdble0(file_id, 'nu_min', numin, herror)
      call h5io_wdble0(file_id, 'nu_max', numax, herror)
      call h5io_wdble1(file_id, 'nuj', nuj, herror)
      call h5io_wdble1(file_id, 'gamma', g, herror)
      call h5io_wdble1(file_id, 'distrib', n, herror)
      call h5io_wdble1(file_id, 'jmbs', jmbs, herror)
      call h5io_wdble1(file_id, 'jssc', jssc, herror)
      call h5io_wdble1(file_id, 'jeic', jeic, herror)

      !  ----->   Absorption
      B = 1d3
      numin = 1d8
      numax = 1d21
      do j = 1, numf
         nua(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numf - 1)) )
      end do
      call h5io_wdble1(file_id, 'nua', nua, herror)

      g1 = 1d2
      g2 = 1d5
      p = 3.0d0
      ke = u_e * pwl_norm(mass_e * cLight**2, p - 1d0, g1, g2)
      do j = 1, numg
         g(j) = g1 * (g2 / g1)**(dble(j - 1) / dble(numg - 1))
         n(j) = ke * powlaw_dis(g(j), g1, g2, p)
      end do

      call mbs_absorption(ambs, nua, g, n, B)
      call h5io_wdble1(file_id, 'ambs1', ambs, herror)

      g1 = 1.01d0
      g2 = 1d2
      p = 1.3d0
      ke = u_e * pwl_norm(mass_e * cLight**2, p - 1d0, g1, g2)
      do j = 1, numg
         g(j) = g1 * (g2 / g1)**(dble(j - 1) / dble(numg - 1))
         n(j) = ke * powlaw_dis(g(j), g1, g2, p)
      end do

      call mbs_absorption(ambs, nua, g, n, B)
      call h5io_wdble1(file_id, 'ambs2', ambs, herror)

      g1 = 1d1
      g2 = 1d3
      g3 = 1d5
      q = 1.3d0
      p = 3.0d0
      do j = 1, numg
         g(j) = g1 * (g3 / g1)**(dble(j - 1) / dble(numg - 1))
      end do
      ke = u_e * (g2**(q - p) * pwl_norm(mass_e * cLight**2, q - 1d0, g1, g2) + pwl_norm(mass_e * cLight**2, p - 1d0, g2, g3))
      do j = 1, numg
         if ( g(j) < g2 ) then
            n(j) = ke * g2**(q - p) * powlaw_dis(g(j), g1, g2, q)
         else
            n(j) = ke * powlaw_dis(g(j), g2, g3, p)
         end if
      end do

      call mbs_absorption(ambs, nua, g, n, B)
      call h5io_wdble1(file_id, 'ambs3', ambs, herror)

      call h5io_closef(file_id, herror)
      call h5close_f(herror)

      write(*, "('--> Radiation processes test')")

   end subroutine syn_test
#endif


end program tests
