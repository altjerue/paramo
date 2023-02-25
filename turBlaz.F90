subroutine turBlaz(params_file,output_file,cool_withKN,with_abs)
   use data_types
   use constants
   use params
   use misc
   use pwl_integ
   use hdf5
   use h5_inout
   use SRtoolkit
   use radiation
   use distribs
   implicit none

   character(len=*),intent(in) :: output_file,params_file
   character(len=1) :: no_rad_case
   logical,intent(in) :: cool_withKN,with_abs
   integer,parameter :: nmod=5

   character(len=*),parameter :: screan_head=&
      '| Iteration |        Time |   Time step |       N_tot |'&
      //new_line('A')//&
      ' ---------------------------------------------------------------------',&
      on_screen="(' | ',I9,' | ',ES11.4,' | ',ES11.4,' | ',ES11.4,' |')"

   integer(HID_T) :: file_id,group_id
   integer :: i,j,k,numbins,numdf,numdt,herror,l,ml
   real(dp) :: uB,R,gmin,gmax,numin,numax,tstep,tmax,&
         gamma_bulk,nu_ext,tesc,tlc,L_jet,ninj,omega_j,R_turb,&
         R_turbm,del_ph,gavg,tc0,tc1,tc2
   real(dp) :: B_0,B_rms,gam0,Gamma0,Gamma2,Gammaa,Gammah,n0,p,&
         sig,th,uph,va,zero1,zero2,tc,A0,a,b,c,c1,hplank,&
         P0,e0,hbar2
   real(dp),allocatable,dimension(:,:) :: gdotty,n1,jnut,jmbs,jssc,jeic,&
            ambs,anut,Diff
   real(dp),allocatable,dimension(:) :: freqs,t,Ntot,Inu,g,dt,dg,urad,dfreq
   real(dp),allocatable,dimension(:) :: tempg,tempnu,Ap,Dpp,total,Mgam,ubol,&
         dotgKN,gdot_db,dotgKN2
   logical :: with_cool,cool_after,do_full_radiation,norad,printout,nossc,no_photon_cooling


   !  ####  ###### ##### #    # #####
   ! #      #        #   #    # #    #
   !  ####  #####    #   #    # #    #
   !      # #        #   #    # #####
   ! #    # #        #   #    # #
   !  ####  ######   #    ####  #
   printout=.true.
   call read_params(params_file)
   L_jet=par_L_j
   gamma_bulk=par_gamma_bulk
   gmin=par_gmin
   gmax=par_gmax
   numin=par_numin
   numax=par_numax

   numbins=par_NG
   numdt=par_NT
   numdf=par_NF
   ninj=par_ed2 !efficienty of driving frequency to plasma
   R=par_R !distance from black hole
   R_turbm=par_R0 !tubulence scale multiplier
   del_ph=par_ed1 !fraction of jet energy converted into turbulence
   sig = par_sigma
   no_rad_case=par_lg1!only caclulate radiation on last iteration
   norad = .true.
   nossc = .true.
   no_photon_cooling = .true.
   if(no_rad_case=='F') then
     norad = .false.
   end if

   allocate(t(0:numdt),freqs(numdf),Ntot(0:numdt),g(numbins),dfreq(numdf),&
         dt(numdt),Inu(numdf),dg(numbins),urad(numbins))
   allocate(tempg(numbins),tempnu(numdf),Ap(numbins),Dpp(numbins),ubol(numdt),&
         total(numbins),Mgam(numdt),dotgKN(numbins),dotgKN2(numbins),gdot_db(numbins))
   allocate(n1(numbins,0:numdt),gdotty(numbins,0:numdt),&
         ambs(numdf,numdt),jmbs(numdf,numdt),jnut(numdf,numdt),&
         jssc(numdf,numdt),anut(numdf,numdt),jeic(numdf,numdt),&
         Diff(numbins,0:numdt))

   build_f: do j=1,numdf
      freqs(j)=numin*((numax/numin)**(dble(j-1)/dble(numdf-1)))
      if(j>1) dfreq(j)=freqs(j)-freqs(j-1)
   end do build_f
   dfreq(1)=dfreq(2)

   build_g: do k=1,numbins
      g(k)=gmin*(gmax/gmin)**(dble(k-1)/dble(numbins-1))
      if (k>1) dg(k)=g(k)-g(k-1)
   end do build_g
   dg(1)=dg(2)

   !   # #    # # #####     ####   ####  #    # #####
   !   # ##   # #   #      #    # #    # ##   # #    #
   !   # # #  # #   #      #      #    # # #  # #    #
   !   # #  # # #   #      #      #    # #  # # #    #
   !   # #   ## #   #      #    # #    # #   ## #    #
   !   # #    # #   #       ####   ####  #    # #####
   do_full_radiation = .true.
   with_cool=.true.
   cool_after=.false.
   zero1=1d-200
   zero2=0d0
   t(0)=0d0
   Gamma0=0.05d0*sig + 2.09d0
   Gammaa=0.124d0*sig - 9.5d0
   Gammah = -0.42*sig - 2.46d0
   Gamma2=Gamma0
   va=cLight*dsqrt(sig/(sig+1d0)) !<--- Alfven speed
   omega_j=2d0*Pi*(1d0-dcos(1d0/gamma_bulk))
   R_turb=(R_turbm)*R/(gamma_bulk)
   B_0=dsqrt(L_jet*(4d0*Pi)/(omega_j*(R**2)*(gamma_bulk**2)*cLight))
   B_rms=dsqrt((B_0**2)*2d0)
   uB=(B_rms**2d0)/(8d0*pi)
   uph=del_ph*(gamma_bulk**2)*(1d11)*(1 + ((1-(1/(gamma_bulk**2)))/3))/(4d0*Pi*cLight)
   n0 = 2*uB/(sig* mass_p* (cLight**2))
   c1 = (32d0/(81d0))*dsqrt(3d0)
   hplank=6.62607015d-27
   hbar2=hplank/(2d0*pi)
   P0 = (eCharge**2d0)/((hbar2**2d0)*cLight*2d0*dsqrt(3d0))
   A0 = 3d0 * c1 * sigmaT * P0 *(hbar2**2d0)*(eCharge ** 2d0)*(8*pi) / ((mass_e**3d0) * (cLight ** 4d0))
   if(printout .eqv. .true.) then
     write(*,*) "P0: ",P0
   end if
   if(printout .eqv. .true.) then
     write(*,*) "A0: ",A0
   end if
   a = A0*R_turb*n0 * uB*(mass_e*(cLight**2d0))
   b = (16d0/9d0)*sigmaT*cLight*(uB + uph)
   c = -ninj*va*uB/(2*n0*R_turb)
   !Cooling type
   if((dsqrt((b**2) - 4*a*c) - b)/(2*a)< 1d-20 .or. nossc .eqv. .true. .or. no_photon_cooling .eqv. .true.)then
     !synchrotron and eic
     gam0=(-c/b)**0.5d0
     write(*,*) "NO SSC"
     if(no_photon_cooling .eqv. .true.)then
       write(*,*) "NO PHOTON COOLING"
       !just synchrotron
       b = (16d0/9d0)*sigmaT*cLight*(uB)
       gam0=(-c/b)**0.5d0
     end if
   else
     !synch, ssc, and eic
     gam0=((dsqrt((b**2) - 4*a*c) - b)/(2*a))**0.5d0
   end if
   tc1 = 1/(((16/9)*sigmaT*cLight*(uB + uph)*gam0/(mass_e*(cLight**2))))
   tc2 = 1/((A0*n0*uB*R_turb*(gam0**3))/(mass_e*(cLight**2)))
   !tc0 = tc1 !without ssc
  if((dsqrt((b**2) - 4*a*c) - b)/(2*a)< 1d-20 .or. nossc .eqv. .true. .or. no_photon_cooling .eqv. .true.)then
    !no ssc cooling
    tc0=tc1
    if(no_photon_cooling .eqv. .true.) then
      !just synchrotron
      tc0 =1/((16/9)*sigmaT*cLight*(uB)*gam0/(mass_e*(cLight**2)))
    end if
  else
    !ssc eic and synch coolingtime
    tc0 =1/(((16/9)*sigmaT*cLight*(uB + uph)*gam0/(mass_e*(cLight**2))) + ((A0*n0*uB*R_turb*(gam0**3))/(mass_e*(cLight**2))))
  end if
  if(printout .eqv. .true.) then
    write(*,*) "A0: ", A0, "a: ",a,"b: ",b,"c:",c,"gam0:",gam0,"tc0:",tc0,"tc1:",tc1,"tc2:",tc2,"R_turb",R_turb
    write(*,*) "tc1*gam0",tc1*gam0 ,"tc2*gam0**3",tc2*(gam0)
  end if
   tlc=1d0!R_turb/cLight
   tesc=1d200
   gam0=gam0
   tc=tc0
   tmax=tc*3d0!par_tmax
   th=gam0/3d2 !injection temperature
   tstep=tmax/numdt

   !   ---> Followeing Eqs. (A7) in Zhdankin et al. (2020)
   Ap=(Gammah*gam0 + Gammaa*g)/tc
   Dpp=(Gamma0*(gam0**2d0) + Gamma2*(g**2d0))/tc
   ! Dpp=((g**2d0)/t2)+((gam0**2d0)/t2)


   ! ----->    External radiation field
   nu_ext=par_nu_ext!*gamma_bulk

   !distribution parameters
   gdot_db = (-1d0)*(Ap+ 2d0*Gamma2*g/tc+ 2d0*Dpp/g)
   gdotty(:,0)=((g**2d0)/(gam0*tc))
   Diff(:,0)=2d0*Dpp
   Diff(:2,0)=1d-100
   gdotty(:2,0)=1d-100
   gdot_db(:2) = 1d-100
   n1(:,0)=n0*RMaxwell(g,th)
   Ntot(0)=sum((n1(:,0)+n1(:,0))*dg*0.5d0)
   n1(:,0)=n0*n1(:,0)/Ntot(0)

   ! -----> Output with initial setup
   if(printout .eqv. .true.) then
     write(*,"('--> Simulation setup')")
     write(*,"('L_jet     =',ES15.7)") L_jet
     write(*,"('ninj       =',ES15.7)") ninj
     write(*,"('del_ph         =',ES15.7)") del_ph
     write(*,"('nu_ext    =',ES15.7)") par_nu_ext
     write(*,"('GB     =',ES15.7)") gamma_bulk
     write(*,"('t_dyn     =',ES15.7)") tlc
     write(*,"('t_esc     =',ES15.7)") tesc
     write(*,"('R_turbm     =',ES15.7)") R_turbm
     write(*,"('sig       =',ES15.7)") sig
     write(*,"('R         =',ES15.7)") R
     ! -------------------------------
     write(*,"('tc        =',ES15.7)") tc
     write(*,"('tstep     =',ES15.7)") tstep
     write(*,"('tmax      =',ES15.7)") tmax
     write(*,"('Theta     =',ES15.7)") Th
     write(*,"('va        =',ES15.7)") va
     write(*,"('uph       =',ES15.7)") uph
     Ntot(0)=sum(n1(:,0)*dg)!,mask=n1(:,i)>1d-200)
     write(*,"('Ntot      =',ES15.7)") Ntot(0)
     write(*,"('n0        =',ES15.7)") n0
     write(*,"('gam0        =',ES15.7)") gam0
      write(*,"('sigmat        =',ES15.7)") sigmaT

     write(*,"('--> Calculating the emission')")
     write(*,*) ''
     write(*,"('Using tstep =   ',F5.3)") tstep
     write(*,*) ''
     write(*,*) screan_head
   end if


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
            &             gdotty(:,i-1)+ gdot_db,& !change 1 to i if you want iteration
            &             Diff(:,1-1),&
            &             zeros1D(numbins, .true.),&
            &             tesc,&
            &             tlc,.false.)

      do l=2,numbins
        ! if(g(l)>0 .and. g(l)<1e6) then
        total(l-1)=(n1(l-1,i)+n1(l,i))*dg(l-1)/2d0
        tempg(l-1)=(g(l-1))*total(l-1)
        ! end if
      end do



      Mgam(i)=sum(tempg)/n0
      Mgam(i) = Mgam(i)
      tempg=0
      if(i == numdt) then
        norad = .false.
      end if

      if(norad .eqv. .false.) then
        !  #####    ##   #####  #   ##   ##### #  ####  #    #
        !  #    #  #  #  #    # #  #  #    #   # #    # ##   #
        !  #    # #    # #    # # #    #   #   # #    # # #  #
        !  #####  ###### #    # # ######   #   # #    # #  # #
        !  #   #  #    # #    # # #    #   #   # #    # #   ##
        !  #    # #    # #####  # #    #   #   #  ####  #    #
        !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j)
        do j=1,numdf
           call syn_emissivity(jmbs(j,i),freqs(j),g,n1(:,i),B_0)
           if (with_abs) call syn_absorption(ambs(j,i),freqs(j),g,n1(:,i),B_0)
        end do
        !$OMP END PARALLEL DO

        call RadTrans_blob(Inu,R_turb,jmbs(:,i),ambs(:,i))

        !---> energy density
        call bolometric_integ(freqs,4d0*pi*Inu/cLight,ubol(i))
        if (do_full_radiation) then
          !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(j)
          do j=1,numdf
             call IC_iso_powlaw(jssc(j,i),freqs(j),freqs,Inu,n1(:,i),g)
             call IC_iso_monochrom(jeic(j,i),freqs(j),uph,nu_ext,n1(:,i),g)
             jnut(j,i)=jmbs(j,i)+jssc(j,i)+jeic(j,i)
             anut(j,i)=ambs(j,i)
          end do
          !$OMP END PARALLEL DO
        end if
        !
        ! do j=2,numdf
        !    tempnu(j-1)=(jmbs(j-1,i)+jmbs(j,i))*dfreq(j-1)/2d0
        ! end do
      end if

      !   ####   ####   ####  #      # #    #  ####
      !  #    # #    # #    # #      # ##   # #    #
      !  #      #    # #    # #      # # #  # #
      !  #      #    # #    # #      # #  # # #  ###
      !  #    # #    # #    # #      # #   ## #    #
      !   ####   ####   ####  ###### # #    #  ####
      ! gdotty(:,i)=0d0
      if (with_cool) then
          ! call bolometric_integ(freqs,4d0*pi*Inu/cLight,ubol(i))
          ! call RadTrans_blob(Inu,R_turb,jssc(:,i)+jeic(:,i),anut(:,i))
          ! call bolometric_integ(freqs,4d0*pi*Inu/cLight,ubol(i))
          ! call RadTrans_blob(Inu,R_turb,jmbs(:,i),ambs(:,i))
          call rad_cool_pwl(dotgKN,g,freqs,4d0*pi*Inu/cLight,cool_withKN)
          call rad_cool_mono(dotgKN2,g,nu_ext,uph,cool_withKN)
      end if

      gdotty(:,i)=dotgKN + (4d0/3d0)*sigmaT*cLight*uB*(g**2)/(mass_e*(cLight**2d0)) + dotgKN2
      ! gdotty(:,i)=dotgKN + (g**2)/(gam0*tc)
      ! ! gdotty(:,i)=(4d0/3d0)*sigmaT*cLight*uB*(g**2)
      !
      !
      ! ! gdotty(:,i)=(4d0/3d0)*sigmaT*cLight*(ubol(i) + uph)*(g**2)
      !
      !redefine parameters
      ! if(mod(i, 50) == 0) then
      !   ml = minloc(abs(g - Mgam(i)),DIM=1)
      !   tc = (gdotty(ml,i)/g(ml))**-1d0
      !   Ap=(Gammah*Mgam(i) + Gammaa*g)/tc
      !   Dpp=(Gamma0*(Mgam(i)**2d0) + Gamma2*(g**2d0))/tc
      !   gdot_db = (-1d0)*(Ap + (1d0)*2d0*Gamma2*g/tc + 2d0*Dpp/g)
      !   Diff(:,0) = 2d0*Dpp
      !   Diff(:2,0)=1d-100
      !   gdotty(:2,i)=1d-100
      !   gdot_db(:2) = 1d-100
      ! end if


      ! ----->   N_tot
      Ntot(i)=sum((n1(:,i)+n1(:,i-1))*dg*0.5d0)!,mask=n1(:,i)>1d-200)

      ! write(*,*) 'delta n', Ntot(i) - Ntot(i-1), 'delta other', sum(((2d0*g/(gam0*tc)) - (Gammaa/tc) - (4d0*Gamma2/tc) + (2d0*Gamma2*(gam0**2d0)/((g**2d0)*tc)) + (2d0*Gamma2/tc))*n1(:,i)*dg)*dt(i)
      if ( i == numdt .or. i==1 .or. mod(i,10) ==0 .and. printout .eqv. .true.) &
            write(*,*) i,t(i),"n_0",n0,"percent_N",Ntot(i)/n0,"Ntotal: ",Ntot(i),"gavg: ",Mgam(i),"predicted_gam0: ",gam0
            ! write(*,*) "gdotty(:,i)/gdotty(:,0): ", max(gdotty(:,i)/gdotty(:,0))," min: ", min(gdotty(:,i)/gdotty(:,0))

   end do time_loop


   !  ####    ##   #    # # #    #  ####
   ! #       #  #  #    # # ##   # #    #
   !  ####  #    # #    # # # #  # #
   !      # ###### #    # # #  # # #  ###
   ! #    # #    #  #  #  # #   ## #    #
   !  ####  #    #   ##   # #    #  ####
   if(printout .eqv. .true.) then
     write(*,*) "---> Saving",trim(output_file)
   end if
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
   call h5io_wdble0(group_id,'Gamma_bulk',gamma_bulk,herror)
   call h5io_wdble0(group_id,'sigma',sig,herror)
   call h5io_wdble0(group_id,'gamma_min',gmin,herror)
   call h5io_wdble0(group_id,'gamma_max',gmax,herror)
   call h5io_wdble0(group_id,'u_ext',uph,herror)
   call h5io_wdble0(group_id,'nu_ext',par_nu_ext,herror)
   call h5io_wdble0(group_id,'nu_min',numin,herror)
   call h5io_wdble0(group_id,'nu_max',numax,herror)
   call h5io_closeg(group_id,herror)
   !----->   Saving data
   call h5io_wdble0(file_id,'R_turb',R_turb,herror)
   call h5io_wdble0(file_id,'R',R,herror)
   call h5io_wdble0(file_id,'Lj',L_jet,herror)
   call h5io_wdble0(file_id,'gam0',gam0,herror)
   call h5io_wdble0(file_id,'tc',tc,herror)
   call h5io_wdble0(file_id,'t_esc',tesc,herror)
   call h5io_wdble0(file_id,'Bfield',B_0,herror)
   call h5io_wdble0(file_id,'n0',n0,herror)
   call h5io_wdble0(file_id,'uB',uB,herror)
   call h5io_wdble0(file_id,'uph',uph,herror)
   call h5io_wdble1(file_id,'time',t(1:),herror)
   call h5io_wdble1(file_id,'nu',freqs,herror)
   call h5io_wdble1(file_id,'ubol',ubol(1:),herror)
   call h5io_wdble1(file_id,'dotgkn',dotgKN,herror)
   call h5io_wdble1(file_id,'gamma',g,herror)
   call h5io_wdble2(file_id,'jnut',jnut,herror)
   call h5io_wdble2(file_id,'jmbs',jmbs,herror)
   call h5io_wdble2(file_id,'jssc',jssc,herror)
   call h5io_wdble2(file_id,'jeic',jeic,herror)
   call h5io_wdble2(file_id,'anut',anut,herror)
   call h5io_wdble2(file_id,'ambs',ambs,herror)
   call h5io_wdble2(file_id,'n_e',n1(:,1:),herror)
   call h5io_wdble2(file_id,'dgdt',gdotty(:,1:),herror)

   !----->   Closing output file
   call h5io_closef(file_id,herror)
   call h5close_f(herror)
   if(printout .eqv. .true.) then
     write(*,*) '=======  FINISHED  ======='
     write(*,*) ''
   end if

end subroutine turBlaz
