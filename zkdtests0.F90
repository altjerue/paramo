!some sudo code
!gram looks like a name to call the program by so make this.
!--looked it up still don't know time to ask jesus
program testszkd
  !add all dependencies
  use dist_evol
  use h5_inout
  use radiation
  use anaFormulae
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
    integer :: i, k, numg, numt, herror, l,numf,j
    real(dp) :: g1, g2, gmin, gmax, tacc, qind, R, tmax, tstep, DynT, tesc, th, sig, lva, tc, ll, va, rho ,Gammah,Gammaa,&
     Gamma0,Gamma2, gam0, u_e,B,p,numin,numax
    real(dp), allocatable, dimension(:) :: t, g, dt, C0,D0, Q0, zero1, zero2,Diff, gdotty, dg, total, Ap, Dpp,&
      n,ke,nuj,tempg,Mgam
    real(dp), allocatable, dimension(:, :) :: n1,n2,n3,jmbs
    !then initiallize all inputs taht are scalar
    !next initiallize all matrix and vectors

    !read in data


    numg=128
    numt =300
    numf = 256
    qind=0d0

    g1=1e1
    g2=1e10
    gmin=1.01d0
    gmax =1.5d0 * g2
    R=1e16
    sig=9d-1
    gam0=3d2
    Gamma0=1.8d0
    Gamma2=Gamma0
    Gammah=-3d0
    Gammaa=-8d0
    rho=1d0
    ll=4d2*rho
    va=cLight*((sig/(sig+1d0))**0.5d0)
    lva=ll/va
    !tc=(1d0)*4d0*lva/sig
    tc=1.33d0
    tmax=1.5d0
    tstep=tmax/numt
    B = 1d0
    p = 2.2d0
    numin = 1d8
    numax = 1d28

!
    write(*,*) tc,cLight,sig
    !write(*,*) numg,numt,qind,tstep,tmax,g1,g2,gmin,gmax,R
    gmax = gmax*g2

    allocate(n(numg), g(numg),t(0:numt),dt(numt),Diff(numg),D0(numg), Q0(numg),C0(numg), zero1(numg), zero2(numg),&
     gdotty(numg), dg(numg),total(numg), Ap(numg), Dpp(numg),nuj(numf),tempg(numg),Mgam(0:numt))

    allocate(n1(0:numt, numg),n2(0:numt, numg),n3(0:numt, numg),jmbs(0:numt,numf))

    u_e = 1d0
    write(*,*) "test"
    !ke = u_e * pwl_norm(mass_e * cLight**2, p - 1d0, g1, g2)
    !-- not sure what build_g is doing
    build_g: do k = 1, numg

      g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numg - 1))
      !n(k) = ke * powlaw_dis(g(k), g1, g2, p)
      if(k>1) then
        dg(k-1)=g(k)-g(k-1)

      end if
      !write(*,*) dg(k-1)
    end do build_g


    Ap=(Gammah*gam0 + Gammaa*g)/tc
    !Ap=(Gammaa*g)/tc
    Dpp=(Gamma0*(gam0**2) + Gamma2*(g**2))/tc

    !build_gdotty: do k = 1, numg
      !gdotty(k) = Ap + 2*Gamma2*g(k) + D
    !end do build_gdotty
    gdotty= (-1d0)*(Ap + (1)*2d0*Gamma2*g/tc + 2d0*Dpp/g - (g**2)/(gam0*tc))


    t(0) = 0d0
    !!!!!QUESTION: What's all this? Try to keep your code clean by deleting or
    !!!!!          commenting out those variables that you don't use. That will
    !!!!!          save you time debugging. It will also help others to read/understand
    !!!!!          your code
    zero1 = 1d-200
    !zero2 = 0d0
    C0 = 3.48d-11 ! 4d0 * sigmaT * uB / (3d0 * mass_e * cLight)
  !  tacc = 1d0 / (C0(1) * 10d0**4.5d0) !tesc
  !  tesc=tacc
  !  D0 = 0.5d0 * pofg(g)**2 / tacc
    th = 1d2
    n1(0, :) = RMaxwell_v(g,th)

  !  n2(0, :) = n1(0, :)
    !n3(0, :) =n1(0, :)
    !dynT=R/cLight
    Diff = 2*Dpp
    !total_part=0
    !D_t(0)=D0

    time_loop: do i = 1, numt

      t(i) = tstep * dble(i)
      !dt(i) = t(i) - t(i - 1) !!!NOTE: In this kind of time-step you're using, dt = tstep for all i
      !write(*,*) dt(i)
      !zero2=-1.8*gam0*n1(i-1,:)/tc

      !Q0 = injection_pwl(t(i), tacc, g, g1, g2, qind, 1d0)
      !--- learn how for loops and if statements work in fortran
      !--looks like you can name for loops adn if you do you end with end do NAME
      call FP_FinDif_difu(tstep, g, n1(i - 1, :), n1(i, :), gdotty, Diff, zero2, 1d200, ll / cLight)!what is r/clight???  tlc was added since old
      !call FP_FinDif_difu(dt(i), g, n1(i - 1, :), n1(i, :), zero2, zero2, zero2, 1d200, R / cLight)
      !call FP_FinDif_difu(dt(i), gamma, distroin, distrout, gammadot, diffusion, Injection, escape, R / cLight)
      !call FP_FinDif_difu(dt(i), g, n2(i - 1, :), n2(i, :), (t(i)/dynT)*C0 * pofg(g)**2, zero1, zero2, 1d200, R / cLight)
      !call FP_FinDif_difu(dt(i), g, n3(i - 1, :), n3(i, :), (t(i)/dynT)*C0 * pofg(g)**2, Diff, zero2, 1d200, R / cLight)
      do l=2, numg
        tempg(l-1) = (g(l-1)*n1(i,l-1)+n1(i,l)*g(l))*dg(l-1)/2d0
      end do
      Mgam(i)=sum(tempg)
      tempg=0
      B=Mgam(i)*((16*Pi*sig*(1)*Mgam(i)*mass_e*(cLight**2)/3))**0.5

      write(*,*) B
      do j = 1, numf

         nuj(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numf - 1)) )

         call mbs_emissivity(jmbs(i,j),nuj(j), g, n1(i,:), B)
      end do

      do l=2, numg
        !write(*,*) (n1(i,l-1)+n1(i,l))*dg(l-1)
        !write(*,*) dg(l)

        total(l-1) = (n1(i,l-1)+n1(i,l))*dg(l-1)/2d0
        !write(*,*) total(l-1)
      end do
      write(*,*) sum(total),"iteration: ",i
      !total_part=0
    end do time_loop

    call h5open_f(herror)
    call h5io_createf("/media/sf_vmshare/SSsol_zkd31.h5", file_id, herror)
    call h5io_wdble1(file_id, 'time', t(1:), herror)
    call h5io_wdble1(file_id, 'gamma', g, herror)
    call h5io_wdble1(file_id, 'nu', nuj, herror)
    call h5io_wdble2(file_id, 'dist1', n1(:, :), herror)
    call h5io_wdble2(file_id, 'sync1', jmbs(:, :), herror)
    !call h5io_wdble2(file_id, 'dist0', n2(1:, :), herror)
  !  call h5io_wdble2(file_id, 'dist3', n3(1:, :), herror)
    call h5io_closef(file_id, herror)
    call h5close_f(herror)



    !end function and in theory feel proud
  end subroutine steadyish_state

end program testszkd
