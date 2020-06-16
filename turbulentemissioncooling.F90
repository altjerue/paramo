    program turbulentemission

  use dist_evol
  use h5_inout
  use radiation
  use anaFormulae

  implicit none

  integer :: m

  do m=1,216
    call distributionEmission(m)
  end do

  write(*,*) '===== Finished ====== '

contains

  subroutine distributionEmission(m)

    implicit none
    integer, intent(in):: m
    integer(HID_T) :: file_id
    integer :: numg, numt, numf, i, k, j, l, herror, case
    real(dp) :: del_ph,L_jet,R_turb,p,tstep, tmax, numin, numax, gmin, gmax,&
        sig, gam0, Gamma0, Gamma2, Gammah, Gammaa, va, Rva, kk, tc, R, B_rms,&
        B0, th, thss, n0, Bulk_lorentz, l_rho, rho, nuext,Uph,Uph_co,tcm,t2,&
        tcorg,mfp,ninj,cool,L_jetm,Rm,R_turbm,omega_j
    real(dp), allocatable, dimension(:) :: t,dt,Ap,Dpp,Diff,gdotty,g,dg,total,&
        nuj, dnuj, tempg,tempnu, Rarray, Mgam, zero1, zero2,U, Inu
    real(dp), allocatable, dimension(:, :) :: n1, jmbs, jic, ambs,jssc
    character(len=50) :: file_name
    character(len=10) :: file_idd
    integer, parameter :: last=7
    real :: dat(last)




      mfp=1.5d0!1.85d0


      numg=128*2
      numt =300
      numf = 700


      tcm=1d0
      gmin=1.000001d0
      gmax =1.5d10*(1d25)
      tmax = tcm*1.5d0 * (1d18)
      tstep = (1/tcm)*1d18


      allocate(g(numg),t(0:numt),dt(numt), zero1(numg), zero2(numg), Diff(numg),&
          gdotty(numg), dg(numg),total(numg), Ap(numg), Dpp(numg),nuj(numf),&
          dnuj(numf), tempnu(numf), tempg(numg),Mgam(0:numt),Rarray(1),U(numt),&
          Inu(numf))

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

        !!get parameters from file:
    write(file_idd, '(i0)') m
    open(1,file='/mnt/c/Users/zachk/Documents/Ubuntu/Inputs/input'// trim(adjustl(file_idd)) //'.txt')
    do i=1,last
        read(1,*) dat(i)
        write(*,*) dat(i)
    end do
    close(1)

      zero1 = 1d-200
      zero2 = 0d0
      t(0) = 0
      !n0=1
      Bulk_lorentz=10

    case=int(dat(1))
    !!!!!WARNING:
    !!!!!   * Try not to use variables with same name as intrinsic functions/subroutines/modules; e.g., case here.
    !!!!!   * Try not to mix types of variables
    !!!!!       * int:  converts a floating point/double into integer
    !!!!!       * dble: converts a real/integer into double precision
    cool=dat(2)
    L_jetm=dat(3)
    Rm=dat(4)
    del_ph=dat(5)
    ninj=dat(6)
    R_turbm=dat(7)




      !zhdankin parameters
    !case=3
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
        tmax=tmax*1d1
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

      L_jet=(L_jetm)*1d46!
      R=(Rm)*1d18
      omega_j=4*Pi/100
      !del_ph=1d-2
      !ninj=1d-1
      R_turb=(R_turbm)*R/(Bulk_lorentz)
      va=cLight*((sig/(sig+1d0))**0.5d0)
      gam0=(3*Pi/4)*sig*(ninj/del_ph)*(r/(Bulk_lorentz*R_turb))*mass_e*(cLight**3)*R/(sigmaT*L_jet*Bulk_lorentz)
      B0=dsqrt(L_jet*(4*Pi)/(omega_j*(R**2)*(Bulk_lorentz**2)*cLight))
      B_rms=dsqrt(2d0)*B0
      n0=(B_rms*3d0)/(sig*16*Pi*gam0*mass_e*(cLight**2))


      Uph=del_ph*(Bulk_lorentz**2)*L_jet/(4*Pi*cLight*(R**2))
      thss=ninj*mass_e*(cLight**2)*sig*va/(16*sigmaT*Uph*R_turb*cLight)
      th = 1d2

    !  thss=75d0
      !gam0=3d2
      Gamma2=Gamma0
      !va=cLight*((sig/(sig+1d0))**0.5d0)
      !R=(1d0)*tc*sig*va/4

      !R=mass_e*cLight*sig*va/(16d0*sigmaT*(thss)*Uph)
    !  R=1d15
    !  tcm=R/((1d0)*tcorg*sig*va/4)
      tc=4d0*R_turb/(ninj*sig*va)
      !thss=thss*tc/tcorg
      !gam0=gam0*tc/tcorg
      !R=(1d0)*tc*sig*va/4

      !B0 = dsqrt(sig*16*Pi*n0*gam0*mass_e*(cLight**2)/(2d0*3d0))
      !rho=gam0*mass_e*(cLight**2)/(eCharge*B_lab)

      !R=l_rho*rho*2*Pi !!stay ar

      rho=R_turb/(2*Pi*l_rho)

      Rva=R/va



      tmax = (4d0)*2d0 * tc
      tstep = tmax/1d3




      t2=(((dsqrt(2d0)/3d0)*((1d0-((va/cLight)**2d0))**-1d0)*((va/cLight)**2d0)*(cLight/(mfp*R_turb)))**(-1))
      !tc=(1/LOG(kk)) +1

      !t2=318874d0

      Rarray(1)=R
      Ap=(Gammah*gam0 + Gammaa*g)/tc
      Dpp=(Gamma0*(gam0**2) + Gamma2*(g**2))/tc
      !Uph=mass_e*cLight*sig*va/(16d0*sigmaT*(thss)*R)

      !Dpp=((g**2d0)/t2)+((gam0**2d0)/t2)


      !distribution parameters
      gdotty= (-1d0)*(1d0)*(Ap + (1)*2d0*Gamma2*g/tc + 2d0*Dpp/g - (g**2)/(gam0*tc))
      !gdotty= (-1d0)*(Ap + (1)*2d0*g/t2 + 2d0*Dpp/g - (g**2)/(gam0*tc))

      Diff = 2*Dpp
      n1(0, :) = n0*RMaxwell_v(g,th)


      p=1d0

      file_name="TEfull"//trim(file_idd)//""
      open(1, file = '/mnt/c/Users/zachk/Documents/Ubuntu/Outputs/output'//trim(file_idd)//'.txt', status = 'new')
      write(*,*) "n0: ", n0,"tmax: ", tmax, "sig: ", sig,"tc: ",tc,"R: ",R,&
          "R_turb",R_turb,"B0: ", B0, "Thetass: ", thss,"gam0: ",gam0,&
          "va: ",va,"Uph",Uph, "Th:",th,"L_jet:",L_jet,"del_ph:",del_ph,&
          "ninj",ninj
      write(1,*) "n0: ", n0,"tmax: ", tmax, "sig: ", sig,"tc: ",tc,"R: ",R,&
          "R_turb",R_turb,"B0: ", B0, "Thetass: ", thss,"gam0: ",gam0,&
          "va: ",va,"Uph",Uph, "Th:",th,"L_jet:",L_jet,"del_ph:",del_ph,&
          "ninj",ninj
      time_loop: do i = 1, numt

        ambs(i,:)=0d0 !!no absorption

        t(i) = tstep * ( (tmax / tstep)**(dble(i - 1) / dble(numt - 1)) )
        dt(i) = t(i) - t(i - 1)
        write(*,*) "test1"
        if(cool>0) then
          if(t(i)>= (R_turb/cLight)*1.5d0)  then
            write(*,*) "COOLING"
            Diff=2*((Gamma0*(gam0**2) + Gamma2*(g**2))/tc)*(((R_turb/cLight)*1.5d0/t(i))**(p))
            gdotty=(-1d0)*(Ap + 2*g*(((R_turb/cLight)*1.5d0/t(i))**(p))/tc + Diff/(g) - (g**2)/(gam0*tc))
            if(sum(n1(i-1,:))<1d-10) then
              Diff=1d-200
              gdotty=1d-200
            end if
          end if
        end if

        call FP_FinDif_difu(dt(i), g, n1(i - 1, :), n1(i, :), gdotty, Diff, zero2, 1d200, R / cLight)

        do l=2, numg
          total(l-1) = (n1(i,l-1)+n1(i,l))*dg(l-1)/2d0
          tempg(l-1) = (g(l-1)*n1(i,l-1)+n1(i,l)*g(l))*dg(l-1)/2d0
        end do

        Mgam(i)=sum(tempg)
        tempg=0
        !B_co=B_lab!*Bulk_lorentz
        Uph_co=Uph!*(Bulk_lorentz**2)
        write(*,*) "test2"
        !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
        !$OMP& PRIVATE(j)

        do j = 1, numf

          if(i==1 .or. i>295) then
           call mbs_emissivity(jmbs(i,j),nuj(j), g, n1(i,:), B0)
          end if
           !!ssc needs to be in a seperate loop
           !call mbs_absorption(ambs(i,j),nuj(j), g, n1(i,:), B_co)



        end do
        !$OMP END PARALLEL DO
        write(*,*) "test3"


        !radiation transfer
       if(i==1 .or. i>295) then
         call RadTrans_blob(Inu, R, jmbs(i, :), ambs(i, :))
       end if

        !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
        !$OMP& PRIVATE(j)
        do j = 1, numf
          if(i==1 .or. i>295 .or. cool>0) then
            call IC_iso_powlaw(jssc(i,j),nuj(j), nuj,Inu, n1(i,:), g)
            call IC_iso_monochrom(jic(i,j), nuj(j), Uph, nuext, n1(i,:), g)
          end if
        end do
       !$OMP END PARALLEL DO

        do j = 2, numf
            tempnu(j-1) = (jmbs(i,j-1)+jmbs(i,j))*dnuj(j-1)/2d0
        end do
          !!!!!energy density
        U(i)=(sum(tempnu)/(nuj(numf)-nuj(1)))*(4d0/3d0)*Pi*(R**3)*(R/cLight)


        write(*,*) "# of particles ",sum(total)/n0,"iteration: ",i,"B0 ", B0,&
            "Mgam ",Mgam(i), " Energy Density: ",U(i),&
            "Magnetic Energy Density: ",((B0)**2)/(8*Pi),"Uph ",Uph_co


      end do time_loop
      write(*,*)"TC: ",tc
      write(*,*)"T2: ",t2

      write(1,*) "# of particles ",sum(total)/n0,"iteration: ",i,&
          "coolingint:",int(tc*cool*1.5d0/tstep),"B0 ", B0, "Mgam ",Mgam(i),&
          " Energy Density: ",U(i),"Magnetic Energy Density: ",((B0)**2)/(8*Pi),&
          "Uph ",Uph_co,"TC: ",tc,"T2: ",t2
      close(1)
      call h5open_f(herror)
      call h5io_createf("/mnt/c/Users/zachk/Documents/Ubuntu/DataFiles/"//trim(file_name)//".h5", file_id, herror)
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
