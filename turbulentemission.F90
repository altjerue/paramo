program turbulentemission

  use dist_evol
  use h5_inout
  use radiation
  use anaFormulae

  implicit none

  call distributionEmission

  write(*,*) '===== Finished ====== '

contains

  subroutine distributionEmission

    implicit none

    integer(HID_T) :: file_id
    integer :: numg, numt, numf, i, k, j, l, herror
    real(dp) :: tstep, tmax, numin, numax, gmin, gmax, sig, gam0, Gamma0, Gamma2, Gammah, Gammaa, va, Rva, kk, tc, R, B_lab, B_co, th, thss, n0, Bulk_lorentz, l_rho, rho
    real(dp), allocatable, dimension(:) :: t, dt, Ap, Dpp, Diff, gdotty, g, dg, total, nuj, tempg, Rarray, Mgam, zero1, zero2
    real(dp), allocatable, dimension(:, :) :: n1, jmbs




    numg=128
    numt =300
    numf = 256

    gmin=1.01d0
    gmax =1.5d10
    tmax = 1.5d0
    tstep = 1d-2


    allocate(g(numg),t(0:numt),dt(numt), zero1(numg), zero2(numg), Diff(numg), gdotty(numg), dg(numg),total(numg), Ap(numg), Dpp(numg),nuj(numf),tempg(numg),Mgam(0:numt),Rarray(1))

    allocate(n1(0:numt, numg), jmbs(0:numt,numf))

    build_g: do k = 1, numg

      g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numg - 1))
      if(k>1) then
        dg(k-1)=g(k)-g(k-1)

      end if
    end do build_g


    zero1 = 1d-200
    zero2 = 0d0
    t(0) = 0
    n0=1
    Bulk_lorentz=10

    !zhdankin parameters
    l_rho = 60.4
    th = 1d2
    thss=75d0
    kk=0.033
    sig=0.9d0
    gam0=3d2
    Gamma0=1.8d0
    Gamma2=Gamma0
    Gammah=-3d0
    Gammaa=-8d0
    va=cLight*((sig/(sig+1d0))**0.5d0)
    tc=(1/LOG(kk)) +1
    !R=(1d0)*tc*sig*va/4
    B_lab = dsqrt(sig*16*Pi*n0*gam0*mass_e*(cLight**2)/3)
    rho=gam0*mass_e*(cLight**2)/(eCharge*B_lab)
    R=l_rho*rho*2*Pi
    Rva=R/va
    Rarray(1)=R
    Ap=(Gammah*gam0 + Gammaa*g)/tc
    Dpp=(Gamma0*(gam0**2) + Gamma2*(g**2))/tc


    !emission parameters
    numin = 1d5
    numax = 1d20

    !distribution parameters
    gdotty= (-1d0)*(Ap + (1)*2d0*Gamma2*g/tc + 2d0*Dpp/g - (g**2)/(gam0*tc))
    Diff = 2*Dpp
    n1(0, :) = RMaxwell_v(g,th)




    write(*,*) "tmax: ", tmax, "sig: ", sig,"tc: ",tc,"R: ",R,"B_lab: ", B_lab, "Theta: ", Th

    time_loop: do i = 1, numt

      t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numt - 1)) )
      dt(i) = t(i) - t(i - 1)

      call FP_FinDif_difu(dt(i), g, n1(i - 1, :), n1(i, :), gdotty, Diff, zero2, 1d200, R / cLight)

      do l=2, numg
        total(l-1) = (n1(i,l-1)+n1(i,l))*dg(l-1)/2d0
        tempg(l-1) = (g(l-1)*n1(i,l-1)+n1(i,l)*g(l))*dg(l-1)/2d0
      end do

      Mgam(i)=sum(tempg)
      tempg=0
      B_co=B_lab/Bulk_lorentz

      do j = 1, numf

         nuj(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numf - 1)) )

         call mbs_emissivity(jmbs(i,j),nuj(j), g, n1(i,:), B_co)

      end do

      write(*,*) "# of particles ",sum(total),"iteration: ",i,"B_co ", B_co, "Mgam ",Mgam(i)


    end do time_loop

    call h5open_f(herror)
    call h5io_createf("/media/sf_vmshare/n0_1_turbulentemission_fig17mT_100.h5", file_id, herror)
    call h5io_wdble1(file_id, 'R', Rarray, herror)
    call h5io_wdble1(file_id, 'time', t(1:), herror)
    call h5io_wdble1(file_id, 'Mgam', Mgam, herror)
    call h5io_wdble1(file_id, 'gamma', g, herror)
    call h5io_wdble1(file_id, 'nu', nuj, herror)
    call h5io_wdble2(file_id, 'dist1', n1(:, :), herror)
    call h5io_wdble2(file_id, 'sync1', jmbs(:, :), herror)
    call h5io_closef(file_id, herror)
    call h5close_f(herror)


  end subroutine distributionEmission
end program turbulentemission
