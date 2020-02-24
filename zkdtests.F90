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
     integer :: i, k, numg, numt
     real(dp) :: g1, g2, gmin, gmax, tacc, qind, R, tmax, tstep
     real(dp), allocatable, dimension(:) :: t, g, dt, C0, Q0, zero1, zero2
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
     R = 1e16

     allocate(g(numg),t(0:numt),dt(numt), zero1(numg), zero2(numg))
     allocate(n1(0:numt, numg))

 !-- not sure what build_g is doing
     build_g: do k = 1, numg
        g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numg - 1))
     end do build_g

     t(0) = 0d0
     zero1 = 1d-200
     zero2 = 0d0
     C0(1) = 3.48d-11 ! 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)
     tacc = 1d0 / (C0(1) * 10d0**4.5d0) !tesc
     n1(0, :) = injection_pwl(1d0, tacc, g, g1, g2, qind, 1d0)

     write(*,*) "compile tests"
     time_loop: do i = 1, numt

        t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numt - 1)) )
        dt(i) = t(i) - t(i - 1)

        Q0 = injection_pwl(t(i), tacc, g, g1, g2, qind, 1d0)
        !--- learn how for loops and if statements work in fortran
         !--looks like you can name for loops adn if you do you end with end do NAME
        call FP_FinDif_difu(dt(i), g, n1(i-1,:), n1(i,:), C0 * pofg(g)**2, zero1, zero2, 1d200, R / cLight)

        write(*,*) "loop iteration"
      end do time_loop
     !end function and in theory feel proud
   end subroutine steadyish_state

end program testszkd
