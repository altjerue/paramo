module benchmarks
!  Description:
!    Benchmarking Paramo
! 
use data_types
use constants
use transformers
use misc
use pwl_integ
use hdf5
use h5_inout
use SRtoolkit
use Aglow_models
use anaFormulae
use dist_evol
use radiation
use specialf
implicit none

! TODO: Steaty state
! TODO: Non-steady state
! TODO: FPE Convergence
! TODO: Measure KN errore
! TODO: Conservation of particles
! TODO: SPN98
! TODO: PM09 & PVP14
! TODO: PK00
private
integer :: Ng, Nt, Nv, i, j, k
real(dp) :: theta_obs, z, d_Lum, pind, g_min, g_max, g1, g2, tmax, tstep, &
      eps_B, eps_e, E0, G0, n_ext, B, b_const, beta_bulk, dr, dt, g2_const, &
      g1_const, eps_g2, gamma_bulk0, mu_obs, uB, urad_const
real(dp), allocatable, dimension(:,:) :: nn, dotg
public afterglow_syn_lc!, steady_state, BB_RadCool

contains

! subroutine steady_state
!    !
!    !  Description:
!    !    This subroutine generates steaty state (SS) solutions 
!    !
!    implicit none
!    integer(HID_T) :: file_id
!    integer :: i, k, numg, numt, herror, numf, j
!    real(dp) :: g1, g2, gmin, gmax, tmax, tstep, qind, tacc, tesc, R, B, numax, numin, uB
!    real(dp), allocatable, dimension(:) :: t, g, Q0, D0, C0, aux0, dt, dg, &
!       Ntot1, Ntot2, Ntot3, Ntot4, Ntot5, Ntot6, zero1, zero2, freqs, &
!       Inu1, Inu4, Inu5, Inu6
!    real(dp), allocatable, dimension(:, :) :: n1, n2, n3, n4, n5, n6
!    real(dp), allocatable, dimension(:, :) :: jmbs1, ambs1, jssc1, jmbs4, ambs4, jssc4, jmbs5, ambs5, jssc5, jmbs6, ambs6, jssc6
!
!    numg = 128
!    numt = 300
!    numf = 192
!    g1 = 1e2
!    g2 = 1e6
!    gmin = 1.01d0
!    gmax = 1.5d0 * g2
!    numin = 1d10
!    numax = 1d27
!    tstep = 1e0
!    tmax = 1e7
!    qind = 0d0
!    R = 1e16
!    B = 1d0
!    uB = B**2 / (8d0 * pi)
!
!    allocate(g(numg), Q0(numg), D0(numg), t(0:numt), dt(numt), dg(numg), &
!       Ntot1(numt), Ntot2(numt), Ntot3(numt), Ntot4(numt), Ntot5(numt), &
!       Ntot6(numt), C0(numg), aux0(numg), zero1(numg), zero2(numg), &
!       freqs(numf), Inu1(numf), Inu4(numf), Inu5(numf), Inu6(numf))
!    allocate(n1(0:numt, numg), n2(0:numt, numg), n3(0:numt, numg), &
!       n4(0:numt, numg), n5(0:numt, numg), n6(0:numt, numg), &
!       jmbs1(numf, numt), ambs1(numf, numt), jssc1(numf, numt), &
!       jmbs4(numf, numt), ambs4(numf, numt), jssc4(numf, numt), &
!       jmbs5(numf, numt), ambs5(numf, numt), jssc5(numf, numt), &
!       jmbs6(numf, numt), ambs6(numf, numt), jssc6(numf, numt))
!
!    build_f: do j=1,numf
!       freqs(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numf - 1)) )
!    end do build_f
!
!    build_g: do k = 1, numg
!       g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numg - 1))
!       if ( k > 1 ) dg(k) = g(k) - g(k - 1)
!    end do build_g
!    dg(1) = dg(2)
!
!    t(0) = 0d0
!    zero1 = 1d-200
!    zero2 = 0d0
!    C0 = 3.48d-11 ! 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)
!    tacc = 1d0 / (C0(1) * 10d0**4.5d0) !tesc
!    tesc = tacc ! 1d200 ! 1.5d0 * R / cLight !
!    D0 = 0.5d0 * pofg(g)**2 / tacc
!    n1(0, :) = injection_pwl(1d0, tacc, g, g1, g2, qind, 1d0)
!    n2(0, :) = n1(0, :)
!    n3(0, :) = n1(0, :)
!    n4(0, :) = n1(0, :)
!    n5(0, :) = n1(0, :)
!    n6(0, :) = n1(0, :)
!
!    write(*, "(4ES15.7)") C0(1), D0(1), tacc, tmax
!
!    time_loop: do i = 1, numt
!
!       t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numt - 1)) )
!       dt(i) = t(i) - t(i - 1)
!
!       Q0 = injection_pwl(t(i), tacc, g, g1, g2, qind, 1d0)
!       ! call FP_FinDif_cool(dt(i), g, n1(i - 1, :), n1(i, :), C0, zero2, 1d200)
!       ! call FP_FinDif_cool(dt(i), g, n2(i - 1, :), n2(i, :), C0, Q0, 1d200)
!       ! call FP_FinDif_cool(dt(i), g, n3(i - 1, :), n3(i, :), C0, Q0, tesc)
!       call FP_FinDif_difu(dt(i), g, n1(i - 1, :), n1(i, :), C0 * pofg(g)**2, zero1, zero2, 1d200, R / cLight)
!       call FP_FinDif_difu(dt(i), g, n2(i - 1, :), n2(i, :), C0 * pofg(g)**2, zero1, Q0,    1d200, R / cLight)
!       call FP_FinDif_difu(dt(i), g, n3(i - 1, :), n3(i, :), C0 * pofg(g)**2, zero1, Q0,    tesc,  R / cLight)
!       call FP_FinDif_difu(dt(i), g, n4(i - 1, :), n4(i, :), C0 * pofg(g)**2, D0,    zero2, 1d200, R / cLight)
!       call FP_FinDif_difu(dt(i), g, n5(i - 1, :), n5(i, :), C0 * pofg(g)**2, D0,    Q0,    1d200, R / cLight)
!       call FP_FinDif_difu(dt(i), g, n6(i - 1, :), n6(i, :), C0 * pofg(g)**2, D0,    Q0,    tesc,  R / cLight)
!
!       Ntot1(i) = sum(n1(i, :) * dg)
!       Ntot2(i) = sum(n2(i, :) * dg)
!       Ntot3(i) = sum(n3(i, :) * dg)
!       Ntot4(i) = sum(n4(i, :) * dg)
!       Ntot5(i) = sum(n5(i, :) * dg)
!       Ntot6(i) = sum(n6(i, :) * dg)
!
!    end do time_loop
!
!    call h5open_f(herror)
!    call h5io_createf("SSsol_dbg.h5", file_id, herror)
!    call h5io_wint0(file_id, 'numdt', numt, herror)
!    call h5io_wint0(file_id, 'numbins', numg, herror)
!    call h5io_wdble0(file_id, 't_max', tmax, herror)
!    call h5io_wdble0(file_id, 'tstep', tstep, herror)
!    call h5io_wdble0(file_id, 'radius', R, herror)
!    call h5io_wdble0(file_id, 'gamma_min', gmin, herror)
!    call h5io_wdble0(file_id, 'gamma_max', gmax, herror)
!    call h5io_wdble0(file_id, 'gamma_1', g1, herror)
!    call h5io_wdble0(file_id, 'gamma_2', g2, herror)
!    call h5io_wdble0(file_id, 'pwl-index', qind, herror)
!    call h5io_wdble0(file_id, 't_acc', tacc, herror)
!    call h5io_wdble0(file_id, 't_esc', tesc, herror)
!    call h5io_wdble1(file_id, 'time', t(1:), herror)
!    call h5io_wdble1(file_id, 'gamma', g, herror)
!    call h5io_wdble1(file_id, 'freqs', freqs, herror)
!    call h5io_wdble2(file_id, 'dist1', n1(1:, :), herror)
!    call h5io_wdble2(file_id, 'dist2', n2(1:, :), herror)
!    call h5io_wdble2(file_id, 'dist3', n3(1:, :), herror)
!    call h5io_wdble2(file_id, 'dist4', n4(1:, :), herror)
!    call h5io_wdble2(file_id, 'dist5', n5(1:, :), herror)
!    call h5io_wdble2(file_id, 'dist6', n6(1:, :), herror)
!    call h5io_wdble2(file_id, 'jsyn1', jmbs1, herror)
!    call h5io_wdble2(file_id, 'jsyn4', jmbs4, herror)
!    call h5io_wdble2(file_id, 'jsyn5', jmbs5, herror)
!    call h5io_wdble2(file_id, 'jsyn6', jmbs6, herror)
!    call h5io_wdble2(file_id, 'asyn1', ambs1, herror)
!    call h5io_wdble2(file_id, 'asyn4', ambs4, herror)
!    call h5io_wdble2(file_id, 'asyn5', ambs5, herror)
!    call h5io_wdble2(file_id, 'asyn6', ambs6, herror)
!    call h5io_wdble2(file_id, 'jssc1', jssc1, herror)
!    call h5io_wdble2(file_id, 'jssc4', jssc4, herror)
!    call h5io_wdble2(file_id, 'jssc5', jssc5, herror)
!    call h5io_wdble2(file_id, 'jssc6', jssc6, herror)
!    call h5io_wdble1(file_id, 'Ntot1', Ntot1, herror)
!    call h5io_wdble1(file_id, 'Ntot2', Ntot2, herror)
!    call h5io_wdble1(file_id, 'Ntot3', Ntot3, herror)
!    call h5io_wdble1(file_id, 'Ntot4', Ntot4, herror)
!    call h5io_wdble1(file_id, 'Ntot5', Ntot5, herror)
!    call h5io_wdble1(file_id, 'Ntot6', Ntot6, herror)
!    call h5io_closef(file_id, herror)
!    call h5close_f(herror)
!
!    write(*, "('--> Fokker-Planck solver test')")
!
! end subroutine steady_state

! subroutine BB_RadCool
!    !
!    !  Description:
!    !    Blackbody radiative cooling
!    !
!    implicit none
!    integer :: Ng, Nf, j, k
!    real(dp) :: T, Theta, lC, xi_c, gmin, gmax, fmin, fmax
!    real(dp), allocatable, dimension(:) :: Ibb, dotg, nu, g
!    Ng = 384
!    Nf = 512
!    allocate(nu(Nf), g(Ng), dotg(Ng), Ibb(Nf))
!    gmin = 1d7
!    gmax = 1e13
!    fmin = 1d5
!    fmax = 1.26d13
!    T = 2.72d0
!    Theta = kBoltz * T / energy_e
!    lC = hPlanck / (mass_e * cLight)
!    xi_c = 4d0 * hPlanck / energy_e
!
!    do k = 1, Ng
!       g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(Ng - 1))
!    end do
!
!    do j = 1, Nf
!       nu(j) = fmin * (fmax / fmin)**(dble(j - 1) / dble(Nf - 1))
!       Ibb(j) = BBintensity(nu(j), T)
!    end do
!
!    call rad_cool(dotg, g, nu, 4 * pi * Ibb / cLight, .true.)
!
!    do k = 1, Ng
!       write(*, *) g(k), dotg(k)
!    end do
!
! end subroutine BB_RadCool

subroutine afterglow_syn_lc
   !
   !  Description:
   !    Producing GRB afterglwo light curves of GRB190114C using models: SPN98,
   !    PM09, PVP14. Full FP evolution with radiation is also produced to
   !    compare with analytic models.
   !
   !  Note:
   !    All quantities are setup as described in SPN98
   !
   implicit none

   real(dp) :: Rd, R0, Q0, Rd2, tinj, tlc, N_e, v_min, v_max
   real(dp), allocatable, dimension(:) :: t, t_obs, v, v_obs, g, Ddiff, R, D, &
         gamma_bulk, volume
   real(dp), allocatable, dimension(:,:) :: j_v, a_v, Qinj, Fv_SPN98, P_syn_max
   
   !---> Parameters
   Nt = 400
   Ng = 384
   Nv = 512
   d_Lum = pc2cm(2.3d9)
   z = 0.4245d0
   pind = 2.6d0
   eps_B = 8d-5
   eps_e = 0.07d0
   E0 = 8d53
   gamma_bulk0 = 700d0
   n_ext = 0.5d0
   theta_obs = 2d0
   v_min = 1d9
   v_max = 1d18

   allocate(t(0:Nt), t_obs(0:Nt), v(Nv), v_obs(Nv), g(Ng), Ddiff(Ng), &
         gamma_bulk(0:Nt), R(0:Nt), D(0:Nt), volume(0:Nt))
   allocate(j_v(Nv, Nt), a_v(Nv, Nt), Qinj(Ng, 0:Nt), dotg(Ng, 0:Nt), &
         nn(Ng, 0:Nt), Fv_SPN98(Nv, 0:Nt), P_syn_max(Nv, Nt))

   !-----> About the observer
   theta_obs = theta_obs * pi / 180d0
   mu_obs = dcos(theta_obs)

   !---> Spherical shock radius and bulk Lorentz factor
   call deceleration_radius(Rd, Rd2, E0, gamma_bulk0, n_ext, .false., 0d0)
   t_obs(0) = 0.9d0 * Rd / (4d0 * gamma_bulk0**2 * cLight)
   call blastwave_approx_SPN98(gamma_bulk0, E0, n_ext, t_obs(0), gamma_bulk(0), R(0), .true.)
   volume(0) = 4d0 * pi * R(0)**3 / 3d0
   D(0) = Doppler(gamma_bulk(0), mu_obs)
   t(0) = t_com_f(t_obs(0), z, gamma_bulk(0), 0d0, mu_obs)

   !---> Magnetic field following eq. (2) in SPN98
   b_const = dsqrt(32d0 * pi * mass_p * eps_B * n_ext) * cLight
   B = b_const * gamma_bulk(0)
   uB = B**2 / (8d0 * pi)

   !---> Minimum and maximum Lorentz factors of the particles distribution
   g2_const = dsqrt(6d0 * pi * eCharge * eps_g2 / sigmaT)
   g1_const = eps_e * mass_p * (pind - 2d0) / ((pind - 1d0) * mass_e)
   g2 = g2_const / dsqrt(B)
   g1 = g1_const * (gamma_bulk(0) - 1d0)

   !---> Total number of swept-up electrons in the postshock SPN98
   N_e = volume(0) * n_ext
   Q0 = N_e / ( g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) * volume(0) )
   ! Q0 = N_e / ((g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) &
   !       - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6)) * volume(0) * energy_e)

   build_freqs: do j = 1, Nv
      v_obs(j) = v_min * ( (v_max / v_min)**(dble(j - 1) / dble(Nv - 1)) )
   end do build_freqs

   build_gammas: do k = 1, Ng
      g(k) = g_min * (g_max / g_min)**(dble(k - 1) / dble(Ng - 1))
   end do build_gammas

   dotg(:, 0) = urad_const * uB * pofg(g)**2
   Ddiff = 1d-200
   Qinj(:, 0) = injection_pwl(t(0), tinj, g, g1, g2, pind, Q0)
   nn(:, 0) = Qinj(:, 0)

   !---> Evolution loop
   evolution: do i = 1, Nt

      !---> Localizing the emission region
      t_obs(i) = t_obs(0) * (tmax / t_obs(0))**(dble(i) / dble(Nt))
      call blastwave_approx_SPN98(gamma_bulk0, E0, n_ext, t_obs(i), gamma_bulk(i), R(i), .true.)
      volume(i) = 4d0 * pi * R(i)**3 / 3d0
      dr = R(i) - R(i - 1)
      D(i) = Doppler(gamma_bulk(i), mu_obs)
      beta_bulk = bofg(gamma_bulk(i))
      t(i) = t(i - 1) + 0.5d0 * dr * ( (1d0 / (beta_bulk * gamma_bulk(i))) &
            + (1d0 / (bofg(gamma_bulk(i - 1)) * gamma_bulk(i - 1))) ) / cLight
      dt = t(i) - t(i - 1)

      !---> FPE solver
      call FP_FinDif_difu(dt, &
            &             g, &
            &             nn(:, i - 1), &
            &             nn(:, i), &
            &             dotg(:, i - 1), &
            &             Ddiff, &
            &             Qinj(:, i - 1), &
            &             1d200, &
            &             tlc)

      !---> Magnetic field following eq. (2) in SPN98
      B = b_const * gamma_bulk(i)
      uB = B**2 / (8d0 * pi)

      !---> Minimum and maximum Lorentz factors of the particles distribution
      g2 = g2_const / dsqrt(B)
      g1 = g1_const * gamma_bulk(i)

      !---> Total number of swept-up electrons in the postshock SPN98
      N_e = volume(i) * n_ext
      Q0 = N_e / ( g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) * volume(i) )
      ! Q0 = N_e / ((g1**(2d0 - pind) * Pinteg(g2 / g1, pind - 1d0, 1d-6) &
      !       - g1**(1d0 - pind) * Pinteg(g2 / g1, pind, 1d-6)) * volume(i))
      Qinj(:, i) = injection_pwl(t(i), tinj, g, g1, g2, pind, Q0)

      !---> Numerical synchrotron
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j, k)
      do j = 1, Nv
         v(j) = nu_com_f(v_obs(j), z, D(i))
         call mbs_emissivity(j_v(j, i), v(j), g, nn(:, i), B)
         call mbs_absorption(a_v(j, i), v(j), g, nn(:, i), B)
      end do
      !$OMP END PARALLEL DO

      !---> Peak spectral power as in eq. (5) in SPN98
      P_syn_max(j, i) = mass_e * cLight**2 * sigmaT * gamma_bulk(i) * B

      !---> Cooling coefficient
      !!!---> Numeric using finite differences
      dotg(:, i) = urad_const * uB * pofg(g)**2 + pofg(g) * dlog(volume(i) / volume(i - 1)) / (3d0 * dt)
      !!!---> MSB00
      !dotg(:, i) = dotg(:, i) + cLight * beta_bulk * gamma_bulk(i) * pofg(gg) / R(i)
      !!!---> Eq. (11) in Hao's paper
      !dotg(:, i) = dotg(:, i) + 1.6d0 * cLight * gamma_bulk(i) * pofg(gg) / R(i)

      !---> On-screen printout
      ! if ( mod(i, nmod) == 0 .or. i == 1 ) &
      !       write(*, on_screen) i, t_obs(i), R(i), gamma_bulk(i), B

   end do evolution

   call syn_afterglow_SPN98(v_obs, t_obs, z, E0, eps_e, eps_B, G0, pind, n_ext, d_Lum, .true., Fv_SPN98)

   !---> Saving data
   call h5open_f(herror)
   call h5io_createf("Afterglow_bench.h5", file_id, herror)
   call h5io_wint0(file_id, 'Nt', Nt, herror)
   call h5io_wint0(file_id, 'Ng', Ng, herror)
   call h5io_wdble0(file_id, 't_max', tmax, herror)
   call h5io_wdble0(file_id, 'tstep', tstep, herror)
   call h5io_wdble0(file_id, 'radius', R, herror)
   call h5io_wdble0(file_id, 'gamma_min', g_min, herror)
   call h5io_wdble0(file_id, 'gamma_max', g_max, herror)
   call h5io_wdble0(file_id, 'gamma_1', g1, herror)
   call h5io_wdble0(file_id, 'gamma_2', g2, herror)
   call h5io_wdble0(file_id, 'pwl-index', pind, herror)
   call h5io_wdble1(file_id, 'time', t(1:), herror)
   call h5io_wdble1(file_id, 'gamma', g, herror)
   call h5io_wdble1(file_id, 'freqs', v_obs, herror)
   call h5io_wdble2(file_id, 'dist', n(1:, :), herror)
   call h5io_wdble2(file_id, 'j_v', j_v, herror)
   call h5io_wdble2(file_id, 'a_v', a_v, herror)
   call h5io_wdble1(file_id, 'n_e', , herror)
   call h5io_wdble1(file_id, 'Ntot2', Ntot2, herror)
   call h5io_wdble1(file_id, 'Ntot3', Ntot3, herror)
   call h5io_wdble1(file_id, 'Ntot4', Ntot4, herror)
   call h5io_wdble1(file_id, 'Ntot5', Ntot5, herror)
   call h5io_wdble1(file_id, 'Ntot6', Ntot6, herror)
   call h5io_closef(file_id, herror)
   call h5close_f(herror)

end subroutine afterglow_syn_lc

end module benchmarks
