program commissorecreation

  use dist_evol
  use h5_inout
  use radiation
  use anaFormulae
  use K2
  use K1

  implicit none

  call distributionEmission

  write(*,*) '===== Finished ====== '

contains

  subroutine distributionEmission

    implicit none

    integer(HID_T) :: file_id
    integer :: numg, numt, numf, i, k, j, l, herror, case
    real(dp) :: tstep, tmax, numin, numax, gmin, gmax, sig, gam0, va, Lva, LL, B_lab, B_co, th, thss, n0, Bulk_lorentz, L_de0, nuext,Uph,Uph_co,tcm,&
      tcorg,mfp,de0,delB0_B0,sigma0,th0,gamth0,w0,wp0,K3,N
    real(dp), allocatable, dimension(:) :: t, dt, Dpp, Diff, gdotty, g, dg, total, nuj, dnuj, tempg,tempnu, Rarray, Mgam, zero1, zero2,U, Inu,&
      wpprime,vparprime, wp, vpar, rho, t2
    real(dp), allocatable, dimension(:, :) :: n1, jmbs, jic, ambs,jssc


    !mfp=1.5d0!1.85d0
    mfp=1d0

    numg=128
    numt =300
    numf = 300


    tcm=1d0
    gmin=1.0001d0
    gmax =1.5d15


    allocate(g(numg),t(0:numt),dt(numt), zero1(numg), zero2(numg), Diff(numg), gdotty(numg), dg(numg),total(numg), Dpp(numg),nuj(numf), dnuj(numf), tempnu(numf), tempg(numg),Mgam(0:numt),Rarray(1),U(numt),&
     Inu(numf),wp(numg),vpar(numg),wpprime(numg),vparprime(numg),rho(numg),t2(numg))

    allocate(n1(0:numt, numg), jmbs(0:numt,numf),jic(0:numt,numf),jssc(0:numt,numf),ambs(0:numt,numf))

    build_g: do k = 1, numg

      g(k) = gmin * (gmax / gmin)**(dble(k - 1) / dble(numg - 1))
      if(k>1) then
        dg(k-1)=g(k)-g(k-1)

      end if
    end do build_g


    !emission parameters
    numin = 1d5
    numax = 1d25
    nuext=1d15
    Bulk_lorentz=10
    build_nuj: do j=1, numf
      nuj(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numf - 1)) )
      if(j>=1) then
        dnuj(j)=nuj(j)-nuj(j-1)

      end if
    end do build_nuj

    !comisso parameters
    L_de0=820
    sigma0=10
    delB0_B0=1
    th0=3d-1
    N=4

    zero1 = 1d-200
    zero2 = 0d0
    t(0) = 0
    n0=1

    sig=sigma0*(delB0_B0**2)
    va=cLight*((sig/(sig+1d0))**0.5d0)
    call K2_init
    K3=K1_func(1/th0) + (4d0/(1/th0))*K2_func(1/th0)
    w0=(K3/K2_func(1/th0))
    gamth0=w0-th0
    wp0=dsqrt(4*Pi*n0*(eCharge**2)/(gamth0*mass_e))
    B_lab = dsqrt(sig*4*Pi*w0*mass_e*(cLight**2))
    rho=g*mass_e*(cLight**2)/(eCharge*B_lab)
    de0=cLight/wp0
    LL=L_de0*de0
    Lva=LL/va
    Rarray(1)=LL
    vpar=cLight*dsqrt(1-(1/(g**2)))
    wp=dsqrt(4*Pi*n0*(eCharge**2)/(g*mass_e))
    wpprime=dsqrt(4*Pi*n0*(eCharge**2)/(mass_e))*(-0.5d0*(g**(-1.5d0)))
    vparprime=2*(cLight**2)/(vpar*(g**3))


    tmax = tcm*12d0*LL/cLight
    tstep = (1/tcm)*1d-2

    !distribution parameters
    t2=(((1d0/3d0)*((1d0-((va/cLight)**2d0))**-1d0)*((va/cLight)**2d0)*(cLight*N/(mfp*rho)))**(-1))
    !Dpp=(g**2)*wp/(t2*vpar)
    !Dpp=0.1d0*(sig*(cLight/R)*(g**2))
    Dpp=(g**2)/t2
    !gdotty= (-1d0)*(2d0*g*wp/(t2*vpar) + ((g**2)/t2)*((wpprime/vpar) - wp*vparprime/(vpar**2)) + 2d0*Dpp/g)
    !gdotty=(-1d0)*(2d0*Dpp/g + 0.1d0*(sig*(cLight/R)*(g*2d0)))/mfp
    gdotty= (-1d0)*(1d0/t2 + 2d0*Dpp/g)
    Diff = 2*Dpp!/mfp
    n1(0, :) = n0*RMaxwell_v(g,th0)




    write(*,*) "tmax: ", tmax, "sig: ", sig,"L: ",LL,"B_lab: ", B_lab, "Theta: ", Th,"va: ",va,"sigmaT",sigmaT

    time_loop: do i = 1, numt

      ambs(i,:)=0d0 !!no absorption

      t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numt - 1)) )
      dt(i) = t(i) - t(i - 1)

    !  if(t(i)>= tc*1.5d0)  then
    !    write(*,*) "COOLING"
    !    Diff=zero1
    !    gdotty=((g**2)/(gam0*tc))+((2/3)*g/Rva)
    !  end if

      call FP_FinDif_difu(dt(i), g, n1(i - 1, :), n1(i, :), gdotty, Diff, zero2, 1d200, LL / cLight)

      do l=2, numg
        total(l-1) = (n1(i,l-1)+n1(i,l))*dg(l-1)/2d0
        tempg(l-1) = (g(l-1)*n1(i,l-1)+n1(i,l)*g(l))*dg(l-1)/2d0
      end do

      Mgam(i)=sum(tempg)
      tempg=0
      B_co=B_lab!*Bulk_lorentz
      Uph_co=Uph!*(Bulk_lorentz**2)
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)

      do j = 1, numf


         !call mbs_emissivity(jmbs(i,j),nuj(j), g, n1(i,:), B_co)
         !!ssc needs to be in a seperate loop
         !call mbs_absorption(ambs(i,j),nuj(j), g, n1(i,:), B_co)



      end do
      !$OMP END PARALLEL DO


      !radiation transfer
    !call RadTrans_blob(Inu, R, jmbs(i, :), ambs(i, :))


      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)
      do j = 1, numf
        !call IC_iso_powlaw(jssc(i,j),nuj(j), nuj,Inu, n1(i,:), g)
        !call IC_iso_monochrom(jic(i,j), nuj(j), Uph, nuext, n1(i,:), g)
      end do
    ! $OMP END PARALLEL DO

      do j = 2, numf
          tempnu(j-1) = (jmbs(i,j-1)+jmbs(i,j))*dnuj(j-1)/2d0
      end do
        !!!!!energy density
      U(i)=(sum(tempnu)/(nuj(numf)-nuj(1)))*(4/3)*Pi*(LL**3)*(LL/cLight)


      write(*,*) "# of particles ",sum(total),"iteration: ",i,"B_co ", B_co, "Mgam ",Mgam(i), " Energy Density: ",U(i),"Magnetic Energy Density: ",((B_co)**2)/(8*Pi),"Uph ",Uph_co


    end do time_loop


    call h5open_f(herror)
    call h5io_createf("/media/sf_vmshare/test.h5", file_id, herror)
    call h5io_wdble1(file_id, 'R', Rarray, herror)
    call h5io_wdble1(file_id, 'time', t(1:), herror)
    call h5io_wdble1(file_id, 'Mgam', Mgam, herror)
    call h5io_wdble1(file_id, 'gamma', g, herror)
    call h5io_wdble1(file_id, 'nu', nuj, herror)
    call h5io_wdble1(file_id, 'energydensity', U, herror)
    call h5io_wdble2(file_id, 'dist1', n1(:, :), herror)
    call h5io_wdble2(file_id, 'sync1', jmbs(:, :), herror)
    call h5io_wdble2(file_id, 'IC1', jic(:, :), herror)
    call h5io_wdble2(file_id, 'ssc', jssc(:, :), herror)
    call h5io_wdble2(file_id, 'abs', ambs(:, :), herror)
    call h5io_closef(file_id, herror)
    call h5close_f(herror)


  end subroutine distributionEmission
end program commissorecreation
