!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Benchmarking Paramo
! 
! author: Jesus M. Rueda-Becerril
! last update: 22/Sep/2020
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module benchmarks
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
use K1
use K2
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
integer :: Ng, Nt, Nv
real(dp) :: theta_obs, z, dLum
real(dp) :: pind, gmin, gmax, g1, g2, tmax, tstep
real(dp) :: numin, numax
real(dp) :: eps_B, eps_e, E0, G0, n_ext, u_ext, nu_ext
real(dp), allocatable, dimension(:) :: t, t_obs, nu, nu_obs
real(dp), allocatable, dimension(:,:) :: jnu, anu, 
public steady_state, 
contains

subroutine steady_state
!
!  Description:
!    This subroutine generates steaty state (SS) solutions 
!
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

subroutine BB_RadCool
!
!  Description:
!    Blackbody radiative cooling
!
implicit none
integer :: Ng, Nf, j, k
real(dp) :: T, Theta, lC, xi_c, gmin, gmax, fmin, fmax
real(dp), allocatable, dimension(:) :: Ibb, dotg, nu, g
Ng = 384
Nf = 512
allocate(nu(Nf), g(Ng), dotg(Ng), Ibb(Nf))
gmin = 1d7
gmax = 1e13
fmin = 1d5
fmax = 1.26d13
T = 2.72d0
Theta = kBoltz * T / energy_e
lC = hPlanck / (mass_e * cLight)
xi_c = 4d0 * hPlanck / energy_e

do k = 1, Ng
   g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(Ng - 1))
end do

do j = 1, Nf
   nu(j) = fmin * (fmax / fmin)**(dble(j - 1) / dble(Nf - 1))
   Ibb(j) = BBintensity(nu(j), T)
end do

call rad_cool(dotg, g, nu, 4 * pi * Ibb / cLight, .true.)

do k = 1, Ng
   write(*, *) g(k), dotg(k)
end do

end subroutine BB_RadCool

subroutine MaxwellDist
implicit none
integer :: i
real(dp), dimension(100) :: g
do i=1,100
   g(i) = 1.001d0 * (1d4 / 1.001d0)**(dble(i - 1) / 100d0)
   write(*,*) g(i), RMaxwell(g(i), 100d0)
end do
end subroutine MaxwellDist

subroutine afterglow_lc
!
!  Description:
!    Producing GRB afterglwo light curves of GRB190114C using models: SPN98,
!    PM09, PVP14. Full FP evolution with radiation is also produced to
!    compare with analytic models.
!
implicit none

real(dp) :: Rd, R0

Nt = 400
Ng = 384
Nv = 512
dL = pc2cm(2.3e9)
z = 0.4245
pind = 2.6
epsB = 8e-5
epse = 0.07
E0 = 8e53
G0 = 700.
n0 = 0.5
uext = 7.5e-9
nuext = sapy.eV2Hz(0.02)   

call blastwave_approx(tobs, z, G0, E0, n, Gshk, Rshk, adiabatic)
call syn_aglow_SPN98(nuo, to, z, E0, epse, epsB, G0, pind, n, d_lum, adiab, flux)


end subroutine afterglow_lc

end module benchmarks


program benchmarking
!
!  Description:
!    Running all benchmarks to test code performance, and reproducing all the test
!
use benchmarking
implicit none

call SPN98_afterflow
write(*, *) '=======  FINISHED  ======='
end program benchmarking
