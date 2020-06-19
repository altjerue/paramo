module anaFormulae
   ! ************************************************************************
   !  Analytic expresions and formulas
   !  ================================
   !
   !  - CyclotronLimit     : function
   !  - RJfAGN_Eq340_pow   : function
   !  - RJfAGN_Eq340_emiss : function
   !  - P81_eq8            : function
   !  - P81_trapzd         : subroutine
   !  - P81_qromb          : function
   !  - RMA                : function
   !  - SL2007             : function
   !  - SL2007_table       : subroutine
   !  - SL2007_alt         : function
   !  - RMA_trapzd         : subroutine
   !  - RMA_qromb          : function
   !  - AMA_trapzd         : subroutine
   !  - AMA_qromb          : function
   !
   ! ************************************************************************
   use data_types
   use constants
   use misc
#ifdef BRWN
   use ISO_C_BINDING
#endif
   use pwl_integ
   implicit none

#ifdef BRWN
   interface
      function tgamma(y) bind(c)
        use ISO_C_BINDING
        real(c_double), value :: y
        real(c_double) :: tgamma
      end function tgamma
   end interface
#endif


contains
   !
   !
   !  Equation 4 in Marcowith & Malzac (2003)
   !
   function CyclotronLimit(beta, m, nu_b) result(cyclo)
      implicit none
      integer, intent(in) :: m
      real(dp), intent(in) :: beta,nu_b
      real(dp) :: cyclo
#ifdef BRWN
      cyclo = 8d0 * pi**2 * nu_b * dble( (m + 1) * ( m**(2 * m + 1) ) ) * &
      beta**(2 * m) / tgamma(dble(2 * m + 2))
#else
      cyclo = 8d0 * pi**2 * nu_b * dble( (m + 1) * ( m**(2 * m + 1) ) ) * &
      beta**(2 * m) / dgamma(dble(2 * m + 2))
#endif
      if (cyclo.lt.1d-200) cyclo = 1d-200
   end function CyclotronLimit

   !
   !     Eq. 3.40 from Relativistic Jets from Active Galactic Nuclei
   !
   function RJfAGN_eq340_pow(gam, beta, chi) result(P_nu)
      implicit none
      real(dp), intent(in) :: chi, gam, beta
      real(dp) :: P_nu, chi_new
      chi_new = 2d0 * chi / (3d0 * gam**2)
      !chi_new = 0.8d0 * chi / gam**2
      !P_nu = 1.6d1 * pi * (4d0 * pi * chi / (3d0 * gam**2))**(1d0/3d0) * &
      !     dexp(- 4d0 * pi * chi / (3d0 * gam**2)) &
      !P_nu = 1.6d1 * pi * (chi/gam**2)**(1d0/3d0) * &
      !     dexp(- chi/gam**2) &
      !     / (2.7d1 * dgamma(4d0 / 3d0))
      P_nu = 8d0 * beta**2 * chi_new**(1d0 / 3d0) * dexp(- chi_new) &
#ifdef BRWN
           / (27d0 * tgamma(4d0 / 3d0))
#else
           / (27d0 * dgamma(4d0 / 3d0))
#endif
      if (P_nu .lt. 1d-200) P_nu = 1d-200
   end function RJfAGN_eq340_pow


   function RJfAGN_eq340_emiss(B, n0, nu, g1, g2, q) result(j_nu)
      implicit none
      real(dp), intent(in) :: nu, n0, B, g1, g2, q
      real(dp) :: j_nu, uB, nu0
      uB = B**2 / 8d0 / pi
      nu0 = 3d0 * B * nuconst / 2d0
      if (nu.ge.nu0*g1**2.and.nu.le.nu0*g2**2) then
         j_nu = 4d0 * cLight * (eCharge**2 / (mass_e * cLight**2))**2 * uB * &
               n0 * nu0**((q - 3d0) / 2d0) * nu**((1d0 - q) / 2d0) / 9d0
      else
         j_nu = 1d-200
      end if
   end function RJfAGN_eq340_emiss
   ! =========================================================================


   ! =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
   !
   !  Results from Petrosian (1981)
   !  =============================
   !
   ! :::: Equation 8 in Petrosian (1981) ::::
   function P81_eq8(beta, gam, theta, chi) result(jnu_of_theta)
      implicit none
      real(dp), intent(in) :: beta, theta, chi, gam
      real(dp) :: scrZ_max, m, t, jnu_of_theta
      t = beta * gam * dsin(theta)
      m = chi * (1d0 + t**2) / gam
      scrZ_max = t * dexp( 1d0 / dsqrt(1d0 + t**2) ) / ( 1d0 + dsqrt(1d0 + t**2) )
      jnu_of_theta = dsqrt(pi * chi) * ( (1d0 + 2d0 / (dtan(theta) * gam)**2) * &
      (1d0 - (beta * dcos(theta))**2)**(2.5d-1) ) * scrZ_max**(2d0 * m) / gam
   end function P81_eq8


   ! :::: Trapezoidal integrator over viewing angles ::::
   subroutine P81_trapzd(theta_a, theta_b, beta, gam, chi, s, n)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: theta_a, theta_b, beta, gam, chi
      real(dp), intent(inout) :: s
      integer :: it,i
      real(dp) :: del,fsum,fa,fb,th

      if (n == 1) then
         fa = P81_eq8(beta, gam, theta_a, chi)
         fb = P81_eq8(beta, gam, theta_b, chi)
         s = 5d-1 * (theta_b - theta_a) * (fa + fb)
      else
         it = 2**(n - 2)
         del = (theta_b - theta_a) / dble(it)
         th = theta_a + 5d-1 * del
         fsum = 0d0
         do i=1,it
            fsum = fsum + P81_eq8(beta, gam, th, chi)
            th = th + del
         end do
         s = 5d-1 * (s + del * fsum)
      end if

   end subroutine P81_trapzd


   ! :::: Romberg integrator over viewing angles ::::
   function P81_qromb(theta_a,theta_b,beta,gam,chi)
      implicit none
      real(dp), intent(in) :: theta_a, theta_b, beta, gam, chi
      integer, parameter :: JMAX = 20, JMAXP = JMAX + 1, K = 5, KM = K - 1
      real(dp), parameter :: EPS = 1d-12
      integer :: j
      real(dp) :: P81_qromb, dqromb
      real(dp), dimension(JMAXP) :: h, s
      h(1) = 1d0
      do j=1,JMAX
         call P81_trapzd(theta_a, theta_b, beta, gam, chi, s(j), j)
         if (j >= K) then
            call polint(h(j - KM:j), s(j - KM:j), 0d0, P81_qromb, dqromb)
            if (abs(dqromb).le.EPS*abs(P81_qromb)) return
         end if
         s(j+1) = s(j)
         h(j+1) = 0.25d0 * h(j)
      end do
      call an_error('P81_qromb: too many steps')
   end function P81_qromb
   ! =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+



   ! =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
   !
   !  Rueda-Mimica-Aloy (RMA) function
   !  ================================
   !

   ! ::::: The RMA function :::::
   ! RMAfit(x) = Rsync(x) = 0.5 pi x CS(x)
   function RMA_new(chi, g) result(res)
      implicit none
      real(dp), intent(in) :: chi, g
      real(dp), parameter :: c1 = 3.2180900500625734d-4, &
            c2 = 6.50532122717873d-1, c3 = 1.5579904689804556d1
      real(dp) :: res, x
      if (chi > 0.75d0 / g) then
         x = 2d0 * chi / (3d0 * g**2)
         if ( x < c1 ) then
            res = 1.8084180211028020864d0 * x**(1d0 / 3d0)
         else if (x >= c1 .and. x <= c2) then
            res = dexp( -0.7871626401625178d0 &
                  - 0.7050933708504841d0 * dlog(x) &
                  - 0.35531869295610624d0 * dlog(x)**2 &
                  - 0.06503312461868385d0 * dlog(x)**3 &
                  - 0.0060901233982264096d0 * dlog(x)**4 &
                  - 0.00022764616638053332d0 * dlog(x)**5 )
            ! res = 10d0**( -0.35564612225908254d0 &
            !       - 0.3421635631371654d0 * dlog10(x) &
            !       - 0.18290602166517914d0 * dlog10(x)**2 &
            !       - 0.03776013298031654d0 * dlog10(x)**3 &
            !       - 0.004040039762244288d0 * dlog10(x)**4 &
            !       - 0.0001732560180040394d0 * dlog10(x)**5 )
         else if (x > c2 .and. x <= c3) then
            res = dexp( -0.8236455154570651d0 &
               - 0.831668613094906d0 * dlog(x) &
               - 0.525630345887699d0 * dlog(x)**2 &
               - 0.22039314697105414d0 * dlog(x)**3 &
               + 0.01669179529512499d0 * dlog(x)**4 &
               - 0.028650695862677572d0 * dlog(x)**5 )
            ! res = 10d0**( -0.357998258501421d0 &
            !       - 0.36360117602083497d0 * dlog10(x) &
            !       - 0.21939774566168257d0 * dlog10(x)**2 &
            !       - 0.10439150658509294d0 * dlog10(x)**3 &
            !       + 0.010445217874656604d0 * dlog10(x)**4 &
            !       - 0.012831897130695337d0 * dlog10(x)**5 )
         else
            res = 0.5d0 * pi * dexp(-x) * (1d0 - 11d0 / (18d0 * x))
         end if
      else
         res = 0d0
      end if
   end function RMA_new


   function RMA(chi, g) result(res)
      implicit none
      real(dp), intent(in) :: chi, g
      real(dp) :: res, cs, x
      if (chi > 0.8d0 / g) then
         x = 2d0 * chi / (3d0 * g**2)
         if (x >= 1d2 .and. x <= 7d2) then
            cs = dexp(-x) * x**(-2d0 / 3d0) / ( 0.869d0 * dexp(-x) + x**(1d0 / 3d0) )
         else if ( x > 7d2 ) then
            cs = 0d0
         else
            cs = x**(-2d0 / 3d0) / ( 0.869d0 + x**(1d0 / 3d0) * dexp(x) )
         end if
         res = x * cs !* pi * sqrt(3d0) / 8d0
      else
         res = 0d0
      end if
   end function RMA


   ! ::::: Eq. 16 in Schlickeiser & Lerch (2007) :::::
   function SL07(chi, g) result(res)
      implicit none
      real(dp), intent(in) :: chi, g
      real(dp) :: res, cs, x
      x = 2d0 * chi / (3d0 * g**2)
      if (x >= 1d2 .and. x <= 5d2) then
         cs = dexp(-x) * x**(-2d0 / 3d0) / ( 0.869d0 * dexp(-x) + x**(1d0 / 3d0) )
      else if (x > 5d2) then
         cs = 0d0
      else
         cs = x**(-2d0 / 3d0) / ( 0.869d0 + x**(1d0 / 3d0) * dexp(x) )
      end if
      res = cs * x !* pi * sqrt(3d0) / 8d0
   end function SL07


   function SL07_alt(chi, g) result(res)
      implicit none
      real(dp), intent(in) :: chi, g
      real(dp) :: res, cs, x, a, b, c, infpow
      a = 1.0d0
      b = 1.0d0
      c = 0.5d0
      infpow = 5d0 / 6d0
      x = 2d0 * chi / (3d0 * g**2)
      !!$print*,"a=",a,"  b=",b,"  c=",c
      !!$x = chi / g**2
      if ( x.ge.1d2 .and. x.le.5d2) then
         cs = a * dexp(-x) * x**(-2d0 / 3d0) / ( b * dexp(-x) + c * x**infpow )
      else if ( x.gt.5d2 ) then
         cs = 0d0
      else
         cs = a * x**(-2d0/3d0) / ( b + c * x**infpow * dexp(x) )
      end if
      res = cs * x !* pi * sqrt(3d0) / 8d0
   end function SL07_alt


   !============================================================================
   !  ###### #    # #  ####   ####  # #    # # ##### #   #
   !  #      ##  ## # #      #      # #    # #   #    # #
   !  #####  # ## # #  ####   ####  # #    # #   #     #
   !  #      #    # #      #      # # #    # #   #     #
   !  #      #    # # #    # #    # #  #  #  #   #     #
   !  ###### #    # #  ####   ####  #   ##   #   #     #
   function j_mb(nu, B, n0, gmin, gmax, qq, RMAfunc) result(emiss)
      implicit none
      real(dp), intent(in) :: nu, B, gmin, gmax, qq, n0
      interface
         function RMAfunc(c, g)
            use data_types
            real(dp), intent(in) :: c, g
            real(dp) :: RMAfunc
         end function RMAfunc
      end interface
      real(dp) :: emiss, chi, nu_b, I2
      nu_b = nuconst * B
      chi = nu / nu_b
      I2 = RMA_qromb(chi, qq, dlog(gmin), dlog(gmax), RMAfunc)
      emiss = pi * eCharge**2 * nu_b * n0 * I2 * gmin**qq / (2d0 * cLight)
      ! emiss = jmbconst * nu_b * n0 * I2 * gmin**qq
   end function j_mb

   subroutine RMA_trapzd(chi, q, lga, lgb, s, n, RMAfunc)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: chi, q, lga, lgb
      real(dp), intent(inout) :: s
      interface
         function RMAfunc(c, g)
            use data_types
            real(dp), intent(in) :: c, g
            real(dp) :: RMAfunc
         end function RMAfunc
      end interface
      integer :: it,i
      real(dp) :: del, fsum, lg, fa, fb, ega, egb, eg
      if ( n == 1 ) then
         ega = dexp(lga)
         egb = dexp(lgb)
         fa = ega**(1d0 - q) * RMAfunc(chi, ega)
         fb = egb**(1d0 - q) * RMAfunc(chi, egb)
         s = 0.5d0 * (lgb - lga) * (fa + fb)
      else
         it = 2**(n - 2)
         del = (lgb - lga) / dble(it)
         lg = lga + 0.5d0 * del
         eg = dexp(lg)
         fsum = 0d0
         itloop: do i = 1, it
            fsum = fsum + eg**(1d0 - q) * RMAfunc(chi, eg)
            lg = lg + del
            eg = dexp(lg)
         end do itloop
         s = 0.5d0 * (s + del * fsum)
      end if
   end subroutine RMA_trapzd

   function RMA_qromb(chi, q, lga, lgb, RMAfunc) result(qromb)
      implicit none
      real(dp), intent(in) :: chi, q, lga, lgb
      interface
         function RMAfunc(c, g)
            use data_types
            real(dp), intent(in) :: c, g
            real(dp) :: RMAfunc
         end function RMAfunc
      end interface
      integer, parameter :: JMAX = 60, JMAXP = JMAX + 1, K = 10, KM = K - 1
      real(dp), parameter :: EPS = 1d-5
      real(dp), dimension(JMAXP) :: h, s
      real(dp) :: dqromb, qromb
      integer :: j
      h(1) = 1d0
      do j=1, JMAX
         call RMA_trapzd(chi, q, lga, lgb, s(j), j, RMAfunc)
         if ( j >= K ) then
            call polint(h(j - KM:j), s(j - KM:j), 0d0, qromb, dqromb)
            if ( dabs(dqromb) <= EPS * dabs(qromb) ) return
         end if
         s(j + 1) = s(j)
         h(j + 1) = 0.25d0 * h(j)
      end do
      print*,'RMA_qromb error'
      print*,'chi    =', chi
      print*,'q      =', q
      print*,'ga     =', dexp(lga)
      print*,'gb     =', dexp(lgb)
      print*,'qromb  =', qromb
      print*,'dqromb =', dqromb
      call an_error('RMA_qromb: too many steps')
   end function RMA_qromb
   !============================================================================



   !============================================================================
   !    ##   #####   ####   ####  #####  #####  ##### #  ####  #    #
   !   #  #  #    # #      #    # #    # #    #   #   # #    # ##   #
   !  #    # #####   ####  #    # #    # #    #   #   # #    # # #  #
   !  ###### #    #      # #    # #####  #####    #   # #    # #  # #
   !  #    # #    # #    # #    # #   #  #        #   # #    # #   ##
   !  #    # #####   ####   ####  #    # #        #   #  ####  #    #
   function a_mb(nu, B, n0, gmin, gmax, qq, RMAfunc) result(absor)
      implicit none
      real(dp), intent(in) :: nu, B, gmin, gmax, qq, n0
      interface
         function RMAfunc(c, g)
            use data_types
            real(dp), intent(in) :: c, g
            real(dp) :: RMAfunc
         end function RMAfunc
      end interface
      real(dp) :: absor, chi, nu_b, A2
      nu_b = nuconst * B
      chi = nu / nu_b
      A2 = ARMA_qromb(chi, qq, dlog(gmin), dlog(gmax), RMAfunc)
      absor = pi * eCharge**2 * nu_b * n0 * A2 * gmin**qq / (4d0 * mass_e * cLight * nu**2)
      ! absor = ambconst * nu_b * n0 * A2 * gmin**qq / nu**2
   end function a_mb

   subroutine ARMA_trapzd(chi, q, lga, lgb, s, n, RMAfunc)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: chi, q, lga, lgb
      real(dp), intent(inout) :: s
      interface
         function RMAfunc(c, g)
            use data_types
            real(dp), intent(in) :: c, g
            real(dp) :: RMAfunc
         end function RMAfunc
      end interface
      integer :: it, i
      real(dp) :: del, fsum, lg, fa, fb, ega, egb, eg
      if (n == 1) then
         ega = dexp(lga)
         egb = dexp(lgb)
         fa = ega**(-q) * RMAfunc(chi, ega) * (q + 1d0 + ega**2 / (ega**2 - 1d0))
         fb = egb**(-q) * RMAfunc(chi, egb) * (q + 1d0 + egb**2 / (egb**2 - 1d0))
         s = 0.5d0 * (lgb - lga) * (fa + fb)
      else
         it = 2**(n - 2)
         del = (lgb - lga) / dble(it)
         lg = lga + 0.5d0 * del
         eg = dexp(lg)
         fsum = 0d0
         itloop: do i=1,it
            fsum = fsum + eg**(-q) * RMAfunc(chi, eg) * (q + 1d0 + eg**2 / (eg**2 - 1d0))
            lg = lg + del
            eg = dexp(lg)
         end do itloop
         s = 0.5d0 * (s + del * fsum) ! where del = (lgb - lga) / it
      end if
   end subroutine ARMA_trapzd

   function ARMA_qromb(chi, q, lga, lgb, RMAfunc) result(qromb)
      implicit none
      real(dp), intent(in) :: chi, q, lga, lgb
      interface
         function RMAfunc(c, g)
            use data_types
            real(dp), intent(in) :: c, g
            real(dp) :: RMAfunc
         end function RMAfunc
      end interface
      integer, parameter :: JMAX = 60, JMAXP = JMAX + 1, K = 10, KM = K - 1
      real(dp), parameter :: EPS = 1d-5
      real(dp), dimension(JMAXP) :: h, s
      real(dp) :: dqromb, qromb
      integer :: j
      h(1) = 1d0
      do j = 1, JMAX
         call ARMA_trapzd(chi, q, lga, lgb, s(j), j, RMAfunc)
         if (j >= K) then
            call polint(h(j-KM:j), s(j-KM:j), 0d0, qromb, dqromb)
            if (dabs(dqromb) .le. EPS * dabs(qromb)) return
         end if
         s(j + 1) = s(j)
         h(j + 1) = 0.25d0 * h(j)
      end do
      print*,'ARMA_qromb error'
      print*,'chi    =', chi
      print*,'q      =', q
      print*,'ga     =', dexp(lga)
      print*,'gb     =', dexp(lgb)
      print*,'qromb  =', qromb
      print*,'dqromb =', dqromb
      call an_error('ARMA_qromb: too many steps')
   end function ARMA_qromb
   !============================================================================

end module anaFormulae
