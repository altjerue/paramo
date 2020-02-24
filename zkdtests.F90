!some sudo code
!gram looks like a name to call the program by so make this.
!--looked it up still don't know time to ask jesus
program testszkd
  !add all dependencies
  use data_types
  use constants
  use misc
  use pwl_integ
  use SRtoolkit
  use dist_evol
  use K2
  use hdf5
  use h5_inout
  !not sure what implicit is
  !---implicit none is to keep i,j,k,l,m,n to be assumed to be integeres
  !--- implicit none should always be used
  implicit none

  !add subroutine to do are calculation:
  call steadyish_state
  !write finised ie learn how write function works
  write(*,*) '===== Finished ====== '
  !create the functrion you called earlier to do the calculation
contains
  !contains is some way of organizing subroutines ask jesus
  subroutine steadyish_state
    !reserve memory for all needed variable and outputs? so learn how variable are made in fortran
    !--outputs are managed by using the intent(out)
    !--intent doesn't appear to be needed if it is not a parameter or variable that leaves
    implicit none
    integer(HID_T) :: file_id
    integer :: i, k, numg, numt, herror, l
    real(dp) :: g1, g2, gmin, gmax, tacc, qind, R, tmax, tstep
    real(dp), allocatable, dimension(:) :: t, g, dt, C0,D0, Q0, zero1, zero2,Diff
    real(dp), allocatable, dimension(:, :) :: n1
    !then initiallize all inputs taht are scalar
    !next initiallize all matrix and vectors
    numg=128
    numt=300
    qind = 0d0
    tstep = 1e0
    tmax = 1e8
    g1=1e1
    g2=1e6
    gmin = 1.01d0
    gmax = 1.5d0 * g2
    R = 1e16

    allocate(g(numg),t(0:numt),dt(numt),Diff(numt),D0(numg), Q0(numg),C0(numg), zero1(numg), zero2(numg))
    allocate(n1(0:numt, numg))

    !-- not sure what build_g is doing
    build_g: do k = 1, numg
      g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numg - 1))
    end do build_g




    t(0) = 0d0
    zero1 = 1d-200
    zero2 = 0d0
    C0 = 3.48d-11 ! 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)
    tacc = 1d0 / (C0(1) * 10d0**4.5d0) !tesc
    D0 = 0.5d0 * pofg(g)**2 / tacc
    n1(0, :) = injection_pwl(1d0, tacc, g, g1, g2, qind, 1d0)


    !D_t(0)=D0
    time_loop: do i = 1, numt

      t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numt - 1)) )
      dt(i) = t(i) - t(i - 1)

      Diff = D0*t(i)/tmax
      !Q0 = injection_pwl(t(i), tacc, g, g1, g2, qind, 1d0)
      !--- learn how for loops and if statements work in fortran
      !--looks like you can name for loops adn if you do you end with end do NAME
      call FP_FinDif_difu(dt(i), g, n1(i - 1, :), n1(i, :), C0 * pofg(g)**2, Diff, zero2, 1d200, R / cLight)

    end do time_loop

    call h5open_f(herror)
    call h5io_createf("SSsol_zkd.h5", file_id, herror)
    call h5io_wdble1(file_id, 'time', t(1:), herror)
    call h5io_wdble1(file_id, 'gamma', g, herror)
    call h5io_wdble2(file_id, 'dist1', n1(1:, :), herror)
    call h5io_closef(file_id, herror)
    call h5close_f(herror)

    !end function and in theory feel proud
  end subroutine steadyish_state

end program testszkd
