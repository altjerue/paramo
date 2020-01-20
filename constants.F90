module constants
   use data_types
   implicit none

   !
   !  #    #   ##   ##### #    #
   !  ##  ##  #  #    #   #    #
   !  # ## # #    #   #   ######
   !  #    # ######   #   #    #
   !  #    # #    #   #   #    #
   !  #    # #    #   #   #    #
   !
   real(dp), parameter :: &
      pi = 3.141592653589793238462643d0, &
      twopi = 6.283185307179586476925287d0, &
      sqrtpi = 1.772453850905516027298167d0, &
      sqrt2 = 1.414213562373095048801689d0, &
      isqrt2 = 0.707106781186547524400844d0, &
      lzero = -460.517018598809136803598291d0  ! log(1e-200)


   !
   ! #####  #    # #   #  ####
   ! #    # #    #  # #  #
   ! #    # ######   #    ####
   ! #####  #    #   #        #
   ! #      #    #   #   #    #
   ! #      #    #   #    ####
   !
   real(dp), parameter :: &
      cLight  = 2.99792458d10, &                 ! cm / s
      mass_p = 1.672621898d-24, &                ! g
      mass_e = 9.10938356d-28, &                 ! g
      energy_e = 8.187105776823886d-07, &        ! erg
      energy_p = 1.5032776159851257d-3, &        ! erg
      eCharge = 4.80320467299766d-10, &          ! cm^(3/2) g^(1/2) / s
      sigmaT = 6.6524587158d-25, &               ! 1 / cm^2
      hPlanck = 6.62607004d-27, &                ! erg s
      hbar = 1.0545718d-27, &                    ! erg s
      kBoltz = 1.38064852d-16, &                 ! erg / K
      sigmaSB = 5.670367d-5, &                   ! erg / cm^2 / K^4 / s
      Ggrav = 6.67408d-8, &                      ! c^3 / g / s^2
      eVolt = 1.60218d-12, &                     ! erg
      nuconst = 2.7992491077281560779657886d6, & ! eCharge / (2 * pi * m_e * cLight)
      jmbconst = 6.66456981963510022816d-30, &   ! sqrt(3) * eCharge**2 / (2 * cLight)
      ambconst = 3.65807942558050123993d-3, &    ! sqrt(3) * eCharge**2 / (4 * m_e * cLight)
      Bcritical = 4.414e13, &                    ! mass_e**2 * cLight**3 / (eCharge * hbar)
      mec2_h = 1.235589965126603d20, &           ! m_e c^2 / h
      h_mec2 = 8.093299785722493d-21             ! h / m_e c^2


   !
   !   ##    ####  ##### #####   ####
   !  #  #  #        #   #    # #    #
   ! #    #  ####    #   #    # #    #
   ! ######      #   #   #####  #    #
   ! #    # #    #   #   #   #  #    #
   ! #    #  ####    #   #    #  ####
   !
   real(dp), parameter :: &
      asto_unit = 1495978707d13, &          ! cm
      solar_mass = 1.9884754153381438d33, & ! g
      solar_radius = 6.957d10, &            ! cm
      solar_lum = 3.828d33, &               ! erg / s
      parsec = 3.085677581467192d18, &      ! cm
      jansky = 1d-23, &                     ! erg / s cm^2 Hz
      lightyr = 9.46073d17                  ! cm


   !
   ! #    # #    # #    # #####  ###### #####   ####
   ! ##   # #    # ##  ## #    # #      #    # #
   ! # #  # #    # # ## # #####  #####  #    #  ####
   ! #  # # #    # #    # #    # #      #####       #
   ! #   ## #    # #    # #    # #      #   #  #    #
   ! #    #  ####  #    # #####  ###### #    #  ####
   !
   real(dp), parameter :: &
      chunche_c100g100 = 2.2619939050180366385d-6, &
      chunche_c100g20  = 2.1157699720918349273d-1

end module constants
