#define STEADY_STATE  (1)
#define RAD_PROCS     (2)
#define BLACK_BODY    (3)
#define SYN_AFTERGLOW (4)
#define ODE_SOLVER    (5)
#define TEST_CHOICE (SYN_AFTERGLOW)

program tests
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
   use specialf
   use blastwave
   implicit none

#ifdef HDF5
   integer(HID_T) :: file_id, group_id
   integer :: herror
#endif

#if TEST_CHOICE == 1
   call steady_state
#elif TEST_CHOICE == 2
   call rad_procs
#elif TEST_CHOICE == 3
   call BlackBody_tests(.true.)
#elif TEST_CHOICE == 4
   call syn_afterglow
#elif TEST_CHOICE == 5
   call ode_solver_test(50, 2, 0d0, 5d0, (/ 0d0, 1d0 /), 2)
#endif

contains

   !> This subroutine produces the steady state solutions of the kinetic equation.
   subroutine steady_state
      implicit none
      integer :: i, k, numg, numt, numf, j
      real(dp) :: g1, g2, gmin, gmax, tmax, tstep, qind, tacc, tesc, R, B, numax, numin, uB
      real(dp), allocatable, dimension(:) :: t, g, Q0, D0, C0, aux0, dt, dg, &
         Ntot1, Ntot2, Ntot3, Ntot4, Ntot5, Ntot6, zero1, zero2, freqs, &
         Inu1, Inu4, Inu5, Inu6
      real(dp), allocatable, dimension(:,:) :: n1, n2, n3, n4, n5, n6
      real(dp), allocatable, dimension(:,:) :: jmbs1, ambs1, jssc1, jmbs4, ambs4, jssc4, jmbs5, ambs5, jssc5, jmbs6, ambs6, jssc6

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

         !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j)
         do j = 1, numf
            call syn_emissivity(jmbs1(j, i), freqs(j), g, n1(i, :), B)
            call syn_emissivity(jmbs4(j, i), freqs(j), g, n4(i, :), B)
            call syn_emissivity(jmbs5(j, i), freqs(j), g, n5(i, :), B)
            call syn_emissivity(jmbs6(j, i), freqs(j), g, n6(i, :), B)
            call syn_absorption(ambs1(j, i), freqs(j), g, n1(i, :), B)
            call syn_absorption(ambs4(j, i), freqs(j), g, n4(i, :), B)
            call syn_absorption(ambs5(j, i), freqs(j), g, n5(i, :), B)
            call syn_absorption(ambs6(j, i), freqs(j), g, n6(i, :), B)
         end do
         !$OMP END PARALLEL DO

         call RadTrans_blob(Inu1, R, jmbs4(:, i), ambs1(:, i))
         call RadTrans_blob(Inu4, R, jmbs4(:, i), ambs4(:, i))
         call RadTrans_blob(Inu5, R, jmbs4(:, i), ambs5(:, i))
         call RadTrans_blob(Inu6, R, jmbs6(:, i), ambs6(:, i))

         !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j)
         do j = 1, numf
            call IC_iso_powlaw(jssc1(j, i), freqs(j), freqs, Inu1, n1(i, :), g)
            call IC_iso_powlaw(jssc4(j, i), freqs(j), freqs, Inu4, n4(i, :), g)
            call IC_iso_powlaw(jssc5(j, i), freqs(j), freqs, Inu5, n5(i, :), g)
            call IC_iso_powlaw(jssc6(j, i), freqs(j), freqs, Inu6, n6(i, :), g)
         end do
         !$OMP END PARALLEL DO

      end do time_loop

      write(*, "('--> Fokker-Planck solver test')")

   end subroutine steady_state

   !> Tests with a Blackbody
   subroutine BlackBody_tests(with_kncool)
      implicit none
      logical, intent(in) :: with_kncool
      integer :: Ng, Nf, j, k
      real(dp) :: T, Theta, xi_c, gmin, gmax, fmin, fmax, fbb_max, ubb, g1, g2, p
      real(dp), allocatable, dimension(:) :: Ibb, dotg1, dotg2, nu, g, n, j1, j2
      Ng = 384
      Nf = 512
      allocate(nu(Nf), g(Ng), dotg1(Ng), dotg2(Ng), Ibb(Nf), n(Ng), j1(Nf), j2(Nf))
      gmin = 1.1d0
      gmax = 1d6
      g1 = 1d3
      g2 = 1d5
      p = 2d0
      fmin = 1d8
      fmax = 1d23
      T = 2.72d0
      Theta = 30d0!kBoltz * T_e / energy_e
      xi_c = 4d0 * hPlanck / energy_e
      fbb_max = 2.8214393721220788934d0 * kBoltz * T / hPlanck
      ubb = BBenergy_dens(T)
      write(*,"(3ES14.7)") T, fbb_max, ubb

      do k = 1, Ng
         g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(Ng - 1))
         n(k) = RMaxwell(g(k), Theta)!powlaw_dis(g(k), g1, g2, p)
      end do

      do j = 1, Nf
         nu(j) = fmin * (fmax / fmin)**(dble(j - 1) / dble(Nf - 1))
         Ibb(j) = BBintensity(nu(j), T)
      end do

      open(unit=77, file="ic_emiss.dat", action="write")
      do j=1, Nf
         call IC_iso_powlaw(j1(j), nu(j), nu, Ibb, n, g)
         call IC_iso_monochrom(j2(j), nu(j), ubb, fbb_max, n, g)
         write(77, "(3ES14.7)") nu(j), j1(j), j2(j)
      end do
      close(77)

      open(unit=1, file="cooling.dat", action="write")
      call rad_cool_pwl(dotg1, g, nu, 4 * pi * Ibb / cLight, with_kncool)
      call rad_cool_mono(dotg2, g, fbb_max, ubb, with_kncool)
      do k=1, Ng
         write(1, "(3ES14.7)") g(k), dotg1(k), dotg2(k)
      end do
      close(1)

   end subroutine BlackBody_tests


   subroutine ode_solver_test(num_steps, num_derivs, x0, xf, y0, which_rk)
      implicit none
      integer, intent(in) :: num_steps, num_derivs, which_rk
      real(dp), intent(in) :: x0, xf
      real(dp), dimension(:), intent(in) :: y0
      integer :: i
      real(dp) :: dx
      real(dp), dimension(0:num_steps) :: x
      real(dp), dimension(0:num_steps, num_derivs) :: y, dydx, yerr
      dx = (xf - x0) / dble(num_steps)
      x(0) = x0
      y(0, :) = y0
      call derivs_tests(x(0), y(0, :), dydx(0, :))
      yerr(0, :) = 0d0
      open(77, file="ode_solver.dat", action="write")
      write(77, "(5ES16.7)") x(0), y(0, 1), y(0, 2), yerr(0, 1), yerr(0, 2)
      do i=1, num_steps
         x(i) = x(i - 1) + dx
         call derivs_tests(x(i - 1), y(i - 1, :), dydx(i - 1, :))
         select case (which_rk)
         case (1)
            call rk4(y(i - 1, :), dydx(i - 1, :), x(i - 1), dx, y(i, :), derivs_tests)
            write(77, "(5ES16.7)") x(i), y(i, 1), y(i, 2)
         case (2)
            call rkck(y(i - 1, :), dydx(i - 1, :), x(i - 1), dx, y(i, :), yerr(i, :), derivs_tests)
            write(77, "(5ES16.7)") x(i), y(i, 1), y(i, 2), yerr(i, 1), yerr(i, 2)
         case default
            call an_error("ode_solver_test: choose RK4 (1) or RKCK (2)")
         end select
      end do
      close(77)
   end subroutine ode_solver_test

   subroutine derivs_tests(x, y, dydx)
      implicit none
      real(dp), intent(in) :: x
      real(dp), dimension(:), intent(in) :: y
      real(dp), dimension(:), intent(out) :: dydx
      integer :: ndum
      ndum = assert_eq((/ size(y), size(dydx) /), 'derivs_tests')
      dydx(1) = 3d0 * y(1) + 4d0 * y(2)
      dydx(2) = 3d0 * y(2) - 4d0 * y(1)
   end subroutine derivs_tests


   !> Synchrotron aferglow. Comparison between numerical synchrotron and the
   !! adiabatic blast-wave described in Sari, Piran & Narayan 1998 (SPN98).
   subroutine syn_afterglow
      implicit none
      integer, parameter :: nmod = 50
      character(len=*), parameter :: horiz_line = &
            & '---------------------------------------------------------------------',&
            &screan_head = horiz_line//new_line('A')//&
            &' | Iteration |  Lab. time  |   BW radius |  gamma_bulk |      Bfield |'&
            &//new_line('A')//' '//horiz_line, &
            &on_screen = "(' | ', I9, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' | ', ES11.4, ' |')"
      integer :: numg, numf, numt, i, j
      real(dp) :: n_ext, dr, dt, E0, eps_B, eps_e, B, g1, g2, g2_const, gmin, &
            gmax, eps_g2, G0, L_e, numax, numin, pind, Q0, tmax, tstep, tlc, &
            rb, b_const, r0, rmax
      real(dp), allocatable, dimension(:) :: t, r, Gbulk, Bbulk, gamma_e, nu, &
            Qinj, zeros, t_lab
      real(dp), allocatable, dimension(:,:) :: dotg, n_e, jnut, anut, jnut_b
      character(len=256) :: output_file

      output_file = "afterglow_test.h5"

      numg = 200
      numf = 300
      numt = 200

      allocate(gamma_e(numg), t(0:numt), Gbulk(0:numt), Bbulk(0:numt), &
            r(0:numt), nu(numf), t_lab(0:numt))
      allocate(zeros(numg))
      allocate(n_e(numg, 0:numt), dotg(numg, 0:numt), jnut(numf, numt), &
            anut(numf, numt), jnut_b(numf, numt))
      zeros = zeros1D(numg, .true.)

      ! SPN98 parameters
      pind = 2.5d0
      eps_g2 = 1d0
      G0 = 200d0
      E0 = 1d52
      tstep = 1d-3
      tmax = 1d5
      n_ext = 1d0
      eps_e = 1d-2
      eps_B = 1d-3
      r0 = 1e14
      rmax = 1e17
      r(0) = r0
      t(0) = 0d0
      t_lab(0) = 0d0
      dt = 0.1d0

      Gbulk(0) = adiab_bw(r(0), G0, E0, n_ext)
      Bbulk(0) = bofg(Gbulk(0))

      !--->  Magnetic field
      B = dsqrt(32d0 * pi * energy_p * eps_B * n_ext * (Gbulk(0) - 1d0) * Gbulk(0))

      !---> Minimum and maximum Lorentz factors of the particles distribution
      g1 = eps_e * mass_p * (pind - 2d0) * (Gbulk(0) - 1d0) / ((pind - 1d0) * mass_e)
      g2_const = dsqrt(6d0 * pi * eCharge * eps_g2 / sigmaT)
      b_const = cLight * sigmaT / (6d0 * pi * energy_e)
      g2 = g2_const / dsqrt(B)

      gmin = 1.0001d0
      gmax = g2 * 10d0
      numin = 1e6
      numax = 1e21

      build_f: do i = 1, numf
         nu(i) = numin * (numax / numin)**(dble(i - 1) / dble(numf - 1))
      end do build_f

      build_g: do i = 1, numg
         gamma_e(i) = gmin * (gmax / gmin)**(dble(i - 1) / dble(numg - 1))
      end do build_g

      !---> Fraction of accreted kinetic energy injected into non-thermal electrons
      L_e = eps_e * 4d0 * pi * r(0)**2 * n_ext * energy_p * cLight * Bbulk(0) * Gbulk(0) * (Gbulk(0) - 1d0)
      Q0 = L_e * pwl_norm(4d0 * pi * r(0)**3 * energy_e / (12d0 * Gbulk(0)), pind - 1d0, g1, g2)
      ! Q0 = L_e / ( ( g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) &
      !       - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6) ) &
      !       * (4d0 * pi * r(0)**3 / 3d0) * energy_e )
      ! L_e = 4d0 * Gbulk(0) * n_ext * eps_e
      ! Q0 = L_e * pwl_norm(4d0 * pi * r(0)**3 / 3d0, pind, g1, g2)
      Qinj = injection_pwl(t(0), 1d200, gamma_e, g1, g2, pind, Q0)
      n_e(:, 0) = Qinj
      dotg(:, 0) = (b_const * B**2 * gamma_e + cLight * Bbulk(0) * Gbulk(0) / r(0)) * gamma_e

      write(*, *) ''
      write(*, "('---> Blast-wave parameters')")
      write(*, "('Gamma_0  =', ES15.7)") G0
      write(*, "('n_ext    =', ES15.7)") n_ext
      write(*, "('pind     =', ES15.7)") pind
      write(*, "('E_0      =', ES15.7)") E0
      write(*, "('eps_e    =', ES15.7)") eps_e
      write(*, "('eps_B    =', ES15.7)") eps_B
      ! write(*, *) ''
      write(*, "('---> Calculating the emission')")
      write(*, *) screan_head

      ! ###### #    #  ####  #      #    # ##### #  ####  #    #
      ! #      #    # #    # #      #    #   #   # #    # ##   #
      ! #####  #    # #    # #      #    #   #   # #    # # #  #
      ! #      #    # #    # #      #    #   #   # #    # #  # #
      ! #       #  #  #    # #      #    #   #   # #    # #   ##
      ! ######   ##    ####  ######  ####    #   #  ####  #    #
      do i=1, numt
         r(i) = r0 * ( (rmax / r0)**(dble(i) / dble(numt)) )
         Gbulk(i) = adiab_bw(R(i), G0, E0, n_ext)
         dr = r(i) - r(i - 1)
         Bbulk(i) = bofg(Gbulk(i))
         call rk2_arr(t(i-1), 1d0 / (Bbulk(i-1:i) * Gbulk(i-1:i) * cLight), dr, t(i))
         call rk2_arr(t_lab(i-1), 1d0 / (Bbulk(i-1:i) * cLight), dr, t_lab(i))
         dt = t(i) - t(i - 1)
         rb = r(i) / (12d0 * Gbulk(i))
         tlc = rb / cLight

         call FP_FinDif_difu(dt, &
               &             gamma_e, &
               &             n_e(:, i - 1), &
               &             n_e(:, i), &
               &             dotg(:, i - 1), &
               &             zeros, &
               &             Qinj, &
               &             tlc, &
               &             tlc)

         !--->  Magnetic field
         B = dsqrt(32d0 * pi * energy_p * eps_B * n_ext * (Gbulk(i) - 1d0) * Gbulk(i))

         !---> Minimum and maximum Lorentz factors of the particles distribution
         g1 = eps_e * mass_p * (pind - 2d0) * (Gbulk(i) - 1d0) / ((pind - 1d0) * mass_e)
         g2 = g2_const / dsqrt(B)

         !---> Fraction of accreted kinetic energy injected into non-thermal electrons
         L_e = eps_e * 4d0 * pi * r(i)**2 * n_ext * energy_p * cLight * Bbulk(i) * Gbulk(i) * (Gbulk(i) - 1d0)
         Q0 = L_e * pwl_norm(4d0 * pi * r(i)**2 * energy_e * rb, pind - 1d0, g1, g2)
         ! Q0 = L_e / ((g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) &
         !       - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6)) &
         !       * (4d0 * pi * r(i)**3 / (12d0 * Gbulk(i))) * energy_e)
         ! L_e = 4d0 * Gbulk(i) * n_ext * eps_e
         ! Q0 = L_e * pwl_norm(1d0, pind, g1, g2)
         Qinj = injection_pwl(t(i), 1d200, gamma_e, g1, g2, pind, Q0)
         dotg(:, i) = (b_const * B**2 * gamma_e + cLight * Bbulk(i) * Gbulk(i) / r(i)) * gamma_e

         !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j)
         do j=1, numf
            ! nu(j, i) = nu_com_f(nu_obs(j), z, D(i))
            call syn_emissivity(jnut(j, i), nu(j), gamma_e, n_e(:, i), B)
            call syn_absorption(anut(j, i), nu(j), gamma_e, n_e(:, i), B)
            call syn_broken(nu(j), t_lab(i), Gbulk(i), g1, pind, B, n_ext, jnut_b(j, i))
         end do
         !$OMP END PARALLEL DO

         !   ####  #    #     ####   ####  #####  ###### ###### #    #
         !  #    # ##   #    #      #    # #    # #      #      ##   #
         !  #    # # #  #     ####  #      #    # #####  #####  # #  #
         !  #    # #  # #         # #      #####  #      #      #  # #
         !  #    # #   ##    #    # #    # #   #  #      #      #   ##
         !   ####  #    #     ####   ####  #    # ###### ###### #    #
         if ( mod(i, nmod) == 0 .or. i == 1 ) &
               write(*, on_screen) i, t_lab(i), r(i), Gbulk(i), B

      end do

      write(*, *) horiz_line
      write(*, "('---> Calculating fluxes')")

      !  ####    ##   #    # # #    #  ####
      ! #       #  #  #    # # ##   # #    #
      !  ####  #    # #    # # # #  # #
      !      # ###### #    # # #  # # #  ###
      ! #    # #    #  #  #  # #   ## #    #
      !  ####  #    #   ##   # #    #  ####
#ifdef HDF5
      ! write(*, *) ''
      write(*, "('---> Creating HDF5')")
      ! ------  Opening output file  ------
      call h5open_f(herror)
      call h5io_createf(output_file, file_id, herror)
      ! ------  Saving initial parameters  ------
      call h5io_createg(file_id, "Parameters", group_id, herror)
      call h5io_wint0 (group_id, 'numt',        numt, herror)
      call h5io_wint0 (group_id, 'numf',        numf, herror)
      call h5io_wint0 (group_id, 'numg',        numg, herror)
      call h5io_wdble0(group_id, 't_max',       tmax, herror)
      call h5io_wdble0(group_id, 'tstep',       tstep, herror)
      call h5io_wdble0(group_id, 'Gamma_bulk0', G0, herror)
      call h5io_wdble0(group_id, 'gamma_min',   gmin, herror)
      call h5io_wdble0(group_id, 'gamma_max',   gmax, herror)
      call h5io_wdble0(group_id, 'gamma_1',     g1, herror)
      call h5io_wdble0(group_id, 'gamma_2',     g2, herror)
      call h5io_wdble0(group_id, 'pwl-index',   pind, herror)
      call h5io_wdble0(group_id, 'epsilon_e',   eps_e, herror)
      call h5io_wdble0(group_id, 'epsilon_B',   eps_B, herror)
      call h5io_wdble0(group_id, 'nu_min',      numin, herror)
      call h5io_wdble0(group_id, 'nu_max',      numax, herror)
      call h5io_wdble0(group_id, 'E0',          E0, herror)
      call h5io_wdble0(group_id, 'R0',          r0, herror)
      call h5io_wdble0(group_id, 'n_ext',       n_ext, herror)
      call h5io_closeg(group_id, herror)
      ! ------  Saving Numerical data  ------
      call h5io_createg(file_id, "Numeric", group_id, herror)
      call h5io_wdble1(group_id, 't_com',    t(1:), herror)
      call h5io_wdble1(group_id, 't_lab',   t_lab(1:), herror)
      call h5io_wdble1(group_id, 'r',       r(1:), herror)
      call h5io_wdble1(group_id, 'Gamma',   Gbulk(1:), herror)
      call h5io_wdble1(group_id, 'nu',      nu, herror)
      call h5io_wdble1(group_id, 'gamma_e', gamma_e, herror)
      call h5io_wdble1(group_id, 'nu',      nu, herror)
      call h5io_wdble2(group_id, 'jnut',    jnut, herror)
      call h5io_wdble2(group_id, 'anut',    anut, herror)
      call h5io_wdble2(group_id, 'n_e',     n_e(:, 1:), herror)
      call h5io_wdble2(group_id, 'dotg',    dotg(:, 1:), herror)
      call h5io_closeg(group_id, herror)
      ! ------  Saving SPN98 data  ------
      call h5io_createg(file_id, "Analytic", group_id, herror)
      call h5io_wdble2(group_id, "emiss_broken", jnut_b, herror)
      call h5io_closeg(group_id, herror)
      ! ------  Closing output file  ------
      call h5io_closef(file_id, herror)
      call h5close_f(herror)
#endif
      ! write(*,*) ''
      write(*, "('==========  FINISHED  ==========')")
   end subroutine syn_afterglow

end program tests
