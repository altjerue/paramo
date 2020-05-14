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
    integer :: numg, numt, numf, i, k, j, l, herror, case
    real(dp) :: tstep, tmax, numin, numax, gmin, gmax, sig, gam0, Gamma0, Gamma2, Gammah, Gammaa, va, Rva, kk, tc, R, B_lab, B_co, th, thss, n0, Bulk_lorentz, l_rho, rho, nuext,Uph,Uph_co,tcm,t2,&
      tcorg,mfp
    real(dp), allocatable, dimension(:) :: t, dt, Ap, Dpp, Diff, gdotty, g, dg, total, nuj, dnuj, tempg,tempnu, Rarray, Mgam, zero1, zero2,U, Inu
    real(dp), allocatable, dimension(:, :) :: n1, jmbs, jic, ambs,jssc


    mfp=1.5d0!1.85d0

    numg=128
    numt =300
    numf = 256


    tcm=1d0
    gmin=1.01d0
    gmax =1.5d10*(1d10)
    tmax = tcm*1.5d0 * (9d0)
    tstep = (1/tcm)*1d-2!*(1d-4)


    allocate(g(numg),t(0:numt),dt(numt), zero1(numg), zero2(numg), Diff(numg), gdotty(numg), dg(numg),total(numg), Ap(numg), Dpp(numg),nuj(numf), dnuj(numf), tempnu(numf), tempg(numg),Mgam(0:numt),Rarray(1),U(numt), Inu(numf))

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
    build_nuj: do j=1, numf
      nuj(j) = numin * ( (numax / numin)**(dble(j - 1) / dble(numf - 1)) )
      if(j>=1) then
        dnuj(j)=nuj(j)-nuj(j-1)

      end if
    end do build_nuj


    zero1 = 1d-200
    zero2 = 0d0
    t(0) = 0
    n0=1
    Bulk_lorentz=10

    !zhdankin parameters
  case=2
    if ( case==1 ) then
      l_rho = 28.3d0!m
      kk=0.033d0
      sig=0.9d0
      Gamma0=1.8d0
      Gammah=-3d0
      Gammaa=-8d0
      tcorg=0.45512721235884102
    end if

    if ( case==2 ) then
      l_rho = 29.6d0!t
      !l_rho=21.4
      kk=0.014d0
      sig=0.2d0
      Gamma0=3d0
      Gammah=-1d0
      Gammaa=-20d0
      tcorg=7.6608494666808964
    end if

    if ( case==3 ) then
      !l_rho = 38.9d0!b
      l_rho=39.1d0
      kk=0.058d0
      sig=3.4d0
      Gamma0=1.4d0
      Gammah=-5d0
      Gammaa=-3.5d0
      tcorg=6.7050164467544707E-002
    end if

    th = 1d2
    thss=75d0
    gam0=3d2
    Gamma2=Gamma0
    va=cLight*((sig/(sig+1d0))**0.5d0)

    !R=(1d0)*tc*sig*va/4

    B_lab = dsqrt(sig*16*Pi*n0*gam0*mass_e*(cLight**2)/(2d0*3d0))
    rho=gam0*mass_e*(cLight**2)/(eCharge*B_lab)

    R=l_rho*rho*2*Pi !!stay ar

    Rva=R/va

    tc=4d0*R/(sig*va)






    t2=(((dsqrt(2d0)/3d0)*((1d0-((va/cLight)**2d0))**-1d0)*((va/cLight)**2d0)*(cLight/(mfp*R)))**(-1))
    !tc=(1/LOG(kk)) +1

    t2=t2

    Rarray(1)=R
    Ap=(Gammah*gam0 + Gammaa*g)/tc
    Dpp=(Gamma0*(gam0**2) + Gamma2*(g**2))/tc
    Uph=mass_e*cLight*sig*va/(16d0*sigmaT*(thss)*R)

    !Dpp=((g**2d0)/t2)+((gam0**2d0)/t2)

    !distribution parameters
    gdotty= (-1d0)*(Ap + (1)*2d0*Gamma2*g/tc + 2d0*Dpp/g - (g**2)/(gam0*tc))
    !gdotty= (-1d0)*(Ap + (1)*2d0*g/t2 + 2d0*Dpp/g - (g**2)/(gam0*tc))
    Diff = 2*Dpp
    n1(0, :) = RMaxwell_v(g,th)




    write(*,*) "tmax: ", tmax, "sig: ", sig,"tc: ",tc,"R: ",R,"B_lab: ", B_lab, "Theta: ", Th,"va: ",va,"sigmaT",sigmaT

    time_loop: do i = 1, numt

      ambs(i,:)=0d0 !!no absorption

      t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numt - 1)) )
      dt(i) = t(i) - t(i - 1)
      write(*,*) "test1"
      call FP_FinDif_difu(dt(i), g, n1(i - 1, :), n1(i, :), gdotty, Diff, zero2, 1d200, R / cLight)

      do l=2, numg
        total(l-1) = (n1(i,l-1)+n1(i,l))*dg(l-1)/2d0
        tempg(l-1) = (g(l-1)*n1(i,l-1)+n1(i,l)*g(l))*dg(l-1)/2d0
      end do

      Mgam(i)=sum(tempg)
      tempg=0
      B_co=B_lab!*Bulk_lorentz
      Uph_co=Uph!*(Bulk_lorentz**2)
      write(*,*) "test2"
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j)

      do j = 1, numf


         !call mbs_emissivity(jmbs(i,j),nuj(j), g, n1(i,:), B_co)
         !!ssc needs to be in a seperate loop
         !call mbs_absorption(ambs(i,j),nuj(j), g, n1(i,:), B_co)



      end do
      !$OMP END PARALLEL DO
      write(*,*) "test3"


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
      U(i)=(sum(tempnu)/(nuj(numf)-nuj(1)))*(4/3)*Pi*(R**3)*(R/cLight)


      write(*,*) "# of particles ",sum(total),"iteration: ",i,"B_co ", B_co, "Mgam ",Mgam(i), " Energy Density: ",U(i),"Magnetic Energy Density: ",((B_co)**2)/(8*Pi),"Uph ",Uph_co


    end do time_loop
    write(*,*)"TC: ",tc
    write(*,*)"T2: ",t2

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
end program turbulentemission
