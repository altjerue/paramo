subroutine turBlaz(params_file, output_file, cool_withKN, with_abs)
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
   use anaFormulae
   use distribs
   use radiation
   implicit none

   character(len=*),intent(in) :: params_file
   logical,intent(in) :: cool_withKN,with_abs
   character(len=*),intent(inout) :: output_file
   integer,parameter :: nmod=50
   character(len=*),parameter :: screan_head=&
      '| Iteration |        Time |   Time step |    nu_0(g2) |       N_tot |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------',&
      on_screen="(' | ',I9,' | ',ES11.4,' | ',ES11.4,' | ',ES11.4,' | ',ES11.4,' |')"
#ifdef HDF5
   integer(HID_T) :: file_id,group_id
#endif
   integer :: i,j,k,numbins,numdf,numdt,time_grid,herror
   integer :: l,mtb_case
   real(dp) :: uB,R,gmin,gmax,numin,numax,pind,D,g1,g2,tstep,Qnth,tmax,&
         d_lum,z,tinj,gamma_bulk,theta_obs,Rdis,mu_obs,nu_ext,tesc,tlc,&
         L_jet,volume,beta_bulk,urad_const,usyn
   real(dp) :: B_0,B_rms,gam0,Gamma0,Gamma2,Gammaa,Gammah,kk,l_rho,mfp,n0,p,&
         rho,Rva,sig,t2,tcm,tcorg,th,thss,uph,va,zero1,zero2,tc
   real(dp),allocatable,dimension(:) :: freqs,t,Ntot,Inu,g,dt,dg,urad,dfreq
   real(dp),allocatable,dimension(:) :: tempg,tempnu,Ap,Dpp,total,Mgam,ubol
   real(dp),allocatable,dimension(:,:) :: gdotty,n1,jnut,jmbs,jssc,jeic,&
         ambs,anut,Qinj,Diff
   real(dp),allocatable,dimension(:,:) :: dotg_temp,dotg_temp2
   character(len=256) :: mtb_label
   logical :: with_cool


   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   call read_params(params_file)
   d_lum=par_d_lum
   z=par_z
   ! tstep=par_tstep
   ! tmax=par_tmax
   L_jet=par_L_j
   gamma_bulk=par_gamma_bulk
   gmin=par_gmin
   gmax=par_gmax
   numin=par_numin
   numax=par_numax
   numbins=par_numbins
   numdt=par_numdt
   numdf=par_numdf
   time_grid=par_time_grid


   allocate(t(0:numdt),freqs(numdf),Ntot(0:numdt),g(numbins),dfreq(numdf),&
         dt(numdt),Inu(numdf),dg(numbins),urad(numbins))
   allocate(tempg(numbins),tempnu(numdf),Ap(numbins),Dpp(numbins),ubol(numdt),&
         total(numbins),Mgam(numdt))
   allocate(n1(numbins,0:numdt),gdotty(numbins,0:numdt),&
         ambs(numdf,numdt),jmbs(numdf,numdt),jnut(numdf,numdt),&
         jssc(numdf,numdt),anut(numdf,numdt),jeic(numdf,numdt),&
         Qinj(numbins,numdt),Diff(numbins,0:numdt))
   allocate(dotg_temp(numbins,numdt),dotg_temp2(numbins,numdt))

   build_f: do j=1,numdf
      freqs(j)=numin*((numax/numin)**(dble(j-1)/dble(numdf-1)))
      if(j>1) dfreq(j)=freqs(j)-freqs(j-1)
   end do build_f
   dfreq(1)=dfreq(2)

   build_g: do k=1,numbins
      g(k)=gmin*(gmax/gmin)**(dble(k-1)/dble(numbins-1))
      ! g(k)=(gmin-1d0)*((gmax-1d0)/(gmin-1d0))**(dble(k-1)/dble(numbins-1))+1d0
      if (k>1) dg(k)=g(k)-g(k-1)
   end do build_g
   dg(1)=dg(2)

   !   # #    # # #####     ####   ####  #    # #####
   !   # ##   # #   #      #    # #    # ##   # #    #
   !   # # #  # #   #      #      #    # # #  # #    #
   !   # #  # # #   #      #      #    # #  # # #    #
   !   # #   ## #   #      #    # #    # #   ## #    #
   !   # #    # #   #       ####   ####  #    # #####
   mfp=1.5d0!1.85d0 !<--- lambda_mfp multiplier
   tcm=1d0

   with_cool=.true.

   zero1=1d-200
   zero2=0d0
   t(0)=0d0
   n0=1d0

   beta_bulk=bofg(gamma_bulk)           ! Beta bulk
   theta_obs=par_theta_obs*pi/180d0 ! Observer viewing angle
   mu_obs=dcos(theta_obs)
   D=Doppler(gamma_bulk,mu_obs)        ! Doppler factor

   ! ------------   zhdankin parameters   ------------
   !    Middle (1),top (2),bottom (3) panels in Fig. 17 from Zhdankin et al. (2020)
   mtb_case=3
   select case(mtb_case)
   case(1)
      l_rho=28.3d0!m
      kk=0.033d0
      sig=0.9d0
      Gamma0=1.8d0
      Gammah=-3d0
      Gammaa=-8d0
      tcorg=0.45512721235884102d0
      mtb_label="m"
   case(2)
      l_rho=29.6d0!t
      !l_rho=21.4
      kk=0.014d0
      sig=0.2d0
      Gamma0=3d0
      Gammah=-1d0
      Gammaa=-20d0
      tcorg=7.6608494666808964d0
      ! tmax=tmax*1d1
      mtb_label="t"
   case(3)
      !l_rho=38.9d0
      l_rho=39.1d0!b
      kk=0.058d0
      sig=3.4d0
      Gamma0=1.4d0
      Gammah=-5d0
      Gammaa=-3.5d0
      tcorg=6.7050164467544707d-2
      mtb_label="b"
   case default
      stop
   end select
   output_file=trim(mtb_label)//trim(output_file)
   ! -------------------------------------------------

   th=1d2
   thss=75d0
   gam0=3d2
   Gamma2=Gamma0
   va=cLight*dsqrt(sig/(sig+1d0)) !<--- Alfven speed

   !R=(1d0)*tc*sig*va/4

   !   ---> Mean magnetic field,B_0,in Zhdankin et al. (2020)
   B_0=dsqrt(sig*16d0*pi*n0*gam0*mass_e*cLight**2/6d0)
   ! B_lab=dsqrt(sig*16d0*pi*n0*gam0*mass_e*cLight**2/(2d0*3d0))
   B_rms=dsqrt(B_0**2*2d0)
   uB=B_0**2/(8d0*pi)

   ! ---> Larmor radius
   rho=gam0*mass_e*cLight**2/(eCharge*B_rms)

   R=l_rho*rho*2d0*pi !!stay ar

   volume=4d0*pi*R**3/3d0
   tlc=R/cLight
   tesc=1d200
   tinj=tlc

   Rva=R/va ! <--- Alfven corssing time

   !   ---> Eq. (A10) in Zhdankin et al. (2020).
   !!!!!NOTE: they assume eta_inj=1
   tc=4d0*R/(sig*va)
   tmax=tc*1d0
   tstep=tc*1e-3

   !   ---> Following Eq. (39) from Comisso & Sironi (2019)
   !!!!!FIXME
   t2=3d0*(mfp*cLight/R)/(gofb(va/cLight)*va/cLight)**2
   ! t2=(((dsqrt(2d0)/3d0)*((1d0-((va/cLight)**2d0))**-1d0)*((va/cLight)**2d0)*(cLight/(mfp*R)))**(-1))

   ! tc=(1d0/dlog(kk)) +1
   ! t2=t2

   !   ---> Followeing Eqs. (A7) in Zhdankin et al. (2020)
   Ap=(Gammah*gam0+Gammaa*g)/tc
   Dpp=(Gamma0*gam0**2+Gamma2*g**2)/tc
   ! Dpp=((g**2d0)/t2)+((gam0**2d0)/t2)

   ! ------------- External radiation field -------------
   ! ---> Following Eq. (7) in Zhdankin et al. (2020)
   uph=mass_e*cLight*sig*va/(16d0*sigmaT*thss *R)
   nu_ext=par_nu_ext
   ! uext=par_uext*gamma_bulk**2*(1d0+beta_bulk**2/3d0) ! Eq. (5.25) Dermer & Menon (2009)
   urad_const=4d0*sigmaT*cLight/(3d0*energy_e)
   ! -------------
   !distribution parameters
   ! gdotty(:,0)=urad_const*(uB+uext)*pofg(g)**2
   call rad_cool_mono(dotg_temp(:,0),g,nu_ext,uph,cool_withKN)
   ! dotg_temp(:,0) = dotg_temp(:,0)
   gdotty(:,0)=-(Ap + (1d0*2d0*Gamma2*g/tc) + (2d0*Dpp/g) - dotg_temp(:,0))! - (g**2/(gam0*tc)))
   !gdotty= (-1d0)*(Ap+(1)*2d0*g/t2+2d0*Dpp/g-(g**2)/(gam0*tc))
   Diff(:,0)=2d0*Dpp

   n1(:,0)=RMaxwell(g,th)
   !!!!!WARNING: The following expression is not normalized
   ! n1(:,0)=n0*dexp(-g/th)/(8d0*pi*(th*mass_e*cLight)**3)

   p=0.5d0


   ! -----> Output with initial setup
   write(*,"('--> Simulation setup')")
   write(*,"('theta_obs =',ES15.7)") par_theta_obs
   write(*,"('Doppler   =',ES15.7)") D
   write(*,"('L_jet     =',ES15.7)") L_jet
   write(*,"('Q_nth     =',ES15.7)") Qnth
   write(*,"('u_B       =',ES15.7)") uB
   write(*,"('B         =',ES15.7)") B_0
   write(*,"('nu_ext    =',ES15.7)") par_nu_ext
   write(*,"('Gamma     =',ES15.7)") gamma_bulk
   write(*,"('t_dyn     =',ES15.7)") tlc
   write(*,"('t_esc     =',ES15.7)") tesc
   write(*,"('t_inj     =',ES15.7)") tinj
   ! -------------------------------
   write(*,"('tc        =',ES15.7)") tc
   write(*,"('tstep     =',ES15.7)") tstep
   write(*,"('tmax      =',ES15.7)") tmax
   write(*,"('sig       =',ES15.7)") sig
   write(*,"('R         =',ES15.7)") R
   write(*,"('Theta     =',ES15.7)") Th
   write(*,"('va        =',ES15.7)") va
   write(*,"('uph       =',ES15.7)") uph
   Ntot(0)=sum(n1(:,0)*dg)!,mask=n1(:,i)>1d-200)
   write(*,"('Ntot      =',ES15.7)") Ntot(0)


   write(*,"('--> Calculating the emission')")
   if (cool_withKN) then
      write(*, "('--> Radiative cooling: Klein-Nishina')")
      output_file="KNcool-"//trim(output_file)
   else
      write(*, "('--> Radiative cooling: Thomson')")
      output_file="Thcool-"//trim(output_file)
   end if
   write(*,*) ''
   write(*,"('Using tstep =',ES10.3)") tstep
   write(*,"('Wrting data in: ',A)") trim(output_file)
   write(*,*) ''
   write(*,*) screan_head


   ! ###### #    #  ####  #      #    # ##### #  ####  #    #
   ! #      #    # #    # #      #    #   #   # #    # ##   #
   ! #####  #    # #    # #      #    #   #   # #    # # #  #
   ! #      #    # #    # #      #    #   #   # #    # #  # #
   ! #       #  #  #    # #      #    #   #   # #    # #   ##
   ! ######   ##    ####  ######  ####    #   #  ####  #    #
   time_loop: do i=1,numdt

      t(i)=tstep*((tmax/tstep)**(dble(i-1)/dble(numdt-1)))
      dt(i)=t(i)-t(i-1)


      !  ###### ###### #####
      !  #      #      #    #
      !  #####  #####  #    #
      !  #      #      #    #
      !  #      #      #    #
      !  ###### ###### #####
      call FP_FinDif_difu(dt(i),&
            &             g,&
            &             n1(:,i-1),&
            &             n1(:,i),&
            &             gdotty(:,0),&
            &             Diff(:,0),&
            &             zeros1D(numbins,.true.),&
            &             tesc,&
            &             tlc)

      do l=2,numbins
         total(l-1)=(n1(l-1,i)+n1(l,i))*dg(l-1)/2d0
         tempg(l-1)=(g(l-1)*n1(l-1,i)+n1(l,i)*g(l))*dg(l-1)/2d0
      end do

      Mgam(i)=sum(tempg)
      tempg=0
      ! B_co=B_lab!*Bulk_lorentz
      ! Uph_co=Uph!*(Bulk_lorentz**2)


      !  #####    ##   #####  #   ##   ##### #  ####  #    #
      !  #    #  #  #  #    # #  #  #    #   # #    # ##   #
      !  #    # #    # #    # # #    #   #   # #    # # #  #
      !  #####  ###### #    # # ######   #   # #    # #  # #
      !  #   #  #    # #    # # #    #   #   # #    # #   ##
      !  #    # #    # #####  # #    #   #   #  ####  #    #
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j)
      do j=1,numdf
         call mbs_emissivity(jmbs(j,i),freqs(j),g,n1(:,i),B_0)
         if (with_abs) call mbs_absorption(ambs(j,i),freqs(j),g,n1(:,i),B_0)
      end do
      !$OMP END PARALLEL DO

      call RadTrans_blob(Inu,R,jmbs(:,i),ambs(:,i))
      !---> energy density
      call bolometric_integ(freqs,4d0*pi*Inu/cLight,ubol(i))
      ! U(i)=(sum(tempnu)/(nuj(numdf)-nuj(1)))*(4/3)*Pi*(R**3)*(R/cLight)

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j)
      do j=1,numdf
         call IC_iso_powlaw(jssc(j,i),freqs(j),freqs,Inu,n1(:,i),g)
         call IC_iso_monochrom(jeic(j,i),freqs(j),uph,nu_ext,n1(:,i),g)
         jnut(j,i)=jmbs(j,i)+jssc(j,i)+jeic(j,i)
         anut(j,i)=ambs(j,i)
      end do
      !$OMP END PARALLEL DO

      do j=2,numdf
         tempnu(j-1)=(jmbs(j-1,i)+jmbs(j,i))*dfreq(j-1)/2d0
      end do

   
      !   ####   ####   ####  #      # #    #  ####
      !  #    # #    # #    # #      # ##   # #    #
      !  #      #    # #    # #      # # #  # #
      !  #      #    # #    # #      # #  # # #  ###
      !  #    # #    # #    # #      # #   ## #    #
      !   ####   ####   ####  ###### # #    #  ####
      ! gdotty(:,i)=0d0
      if (with_cool) then
         call bolometric_integ(freqs,4d0*pi*Inu/cLight,usyn)
      !    ! call RadTrans_blob(Inu,R,jssc(:,i)+jeic(:,i),anut(:,i))
      !    ! call RadTrans_blob(Inu,R,jmbs(:,i),ambs(:,i))
         call rad_cool_pwl(dotg_temp2(:,i),g,freqs,4d0*pi*Inu/cLight,cool_withKN)
      !    call rad_cool_mono(dotg_temp(:,i),g,nu_ext,uph,cool_withKN)
      end if
      ! gdotty(:,i)=gdotty(:,0)!+dotg_temp(:,i)
      dotg_temp(:,i)=dotg_temp(:,0)

      ! ----->   N_tot
      Ntot(i)=sum((n1(:,i)+n1(:,i-1))*dg*0.5d0)!,mask=n1(:,i)>1d-200)

      if (mod(i,nmod)==0.or.i==1) &
            ! write(*,on_screen) i,t(i),ubol(i),Ntot(i)/Ntot(0),Ntot(i)
            write(*,on_screen) i,t(i),tstep,maxval(dotg_temp2(:, i)),Ntot(i)

   end do time_loop

   !!! use the final distribution at t=tc and use it as the injection term for a new evolution

   ! Qinj = n1(:,numdt)
   ! time_loop2: do i=1,numt2
   !    !  ###### ###### #####
   !    !  #      #      #    #
   !    !  #####  #####  #    #
   !    !  #      #      #    #
   !    !  #      #      #    #
   !    !  ###### ###### #####
   !    call FP_FinDif_difu(dt(i),&
   !    &             g,&
   !    &             n1(:,i-1),&
   !    &             n1(:,i),&
   !    &             gdotty(:,0),&
   !    &             zeros1D()!Diff(:,0),&
   !    &             Qinj,&
   !    &             tesc,&
   !    &             tlc)
   ! end do time_loop2


   !  ####    ##   #    # # #    #  ####
   ! #       #  #  #    # # ##   # #    #
   !  ####  #    # #    # # # #  # #
   !      # ###### #    # # #  # # #  ###
   ! #    # #    #  #  #  # #   ## #    #
   !  ####  #    #   ##   # #    #  ####
   write(*,*) "---> Saving"
#ifdef HDF5
   !----->   Opening output file
   call h5open_f(herror)
   call h5io_createf(trim(output_file),file_id,herror)
   !----->   Saving initial parameters
   call h5io_createg(file_id,"Parameters",group_id,herror)
   call h5io_wint0(group_id,'numdt',numdt,herror)
   call h5io_wint0(group_id,'numdf',numdf,herror)
   call h5io_wint0(group_id,'numbins',numbins,herror)
   call h5io_wdble0(group_id,'t_max',tmax,herror)
   call h5io_wdble0(group_id,'tstep',tstep,herror)
   call h5io_wdble0(group_id,'R_b',R,herror)
   call h5io_wdble0(group_id,'R_em',Rdis,herror)
   call h5io_wdble0(group_id,'d_lum',d_lum,herror)
   call h5io_wdble0(group_id,'redshift',z,herror)
   call h5io_wdble0(group_id,'Gamma_bulk',gamma_bulk,herror)
   call h5io_wdble0(group_id,'sigma',sig,herror)
   call h5io_wdble0(group_id,'theta_obs_deg',par_theta_obs,herror)
   call h5io_wdble0(group_id,'gamma_min',gmin,herror)
   call h5io_wdble0(group_id,'gamma_max',gmax,herror)
   call h5io_wdble0(group_id,'gamma_1',g1,herror)
   call h5io_wdble0(group_id,'gamma_2',g2,herror)
   call h5io_wdble0(group_id,'pwl-index',pind,herror)
   call h5io_wdble0(group_id,'u_ext',uph,herror)
   call h5io_wdble0(group_id,'nu_ext',par_nu_ext,herror)
   call h5io_wdble0(group_id,'L_jet',L_jet,herror)
   call h5io_wdble0(group_id,'nu_min',numin,herror)
   call h5io_wdble0(group_id,'nu_max',numax,herror)
   call h5io_closeg(group_id,herror)
   !----->   Saving data
   call h5io_wdble0(file_id,'tc',tc,herror)
   call h5io_wdble0(file_id,'t_inj',tinj,herror)
   call h5io_wdble0(file_id,'t_esc',tesc,herror)
   call h5io_wdble0(file_id,'Bfield',B_0,herror)
   call h5io_wdble0(file_id,'uB',uB,herror)
   call h5io_wdble0(file_id,'uph',uph,herror)
   call h5io_wdble0(file_id,'usyn',usyn,herror)
   call h5io_wdble1(file_id,'time',t(1:),herror)
   call h5io_wdble1(file_id,'nu',freqs,herror)
   call h5io_wdble1(file_id,'gamma',g,herror)
   call h5io_wdble2(file_id,'jnut',jnut,herror)
   call h5io_wdble2(file_id,'jmbs',jmbs,herror)
   call h5io_wdble2(file_id,'jssc',jssc,herror)
   call h5io_wdble2(file_id,'jeic',jeic,herror)
   call h5io_wdble2(file_id,'anut',anut,herror)
   call h5io_wdble2(file_id,'ambs',ambs,herror)
   call h5io_wdble2(file_id,'n_e',n1(:,1:),herror)
   call h5io_wdble2(file_id,'dgdt',gdotty(:,1:),herror)
   call h5io_wdble2(file_id,'cool_coef_eic',dotg_temp,herror)
   call h5io_wdble2(file_id,'cool_coef_ssc',dotg_temp2,herror)
   !----->   Closing output file
   call h5io_closef(file_id,herror)
   call h5close_f(herror)
#endif
   write(*,*) '=======  FINISHED  ======='
   write(*,*) ''

end subroutine turBlaz
