module radiation
   use data_types
   use constants
   use misc
#ifdef INTEL
   use ISO_C_BINDING
#endif
   use pwl_integ
   use SRtoolkit
   !$ use omp_lib
   implicit none

#ifdef INTEL
   interface
      function tgamma(y) bind(c)
        use ISO_C_BINDING
        real(c_double), value :: y
        real(c_double) :: tgamma
      end function tgamma
   end interface
#endif

   interface RadTrans
      module procedure RadTrans_s
      module procedure RadTrans_v
   end interface RadTrans

   interface opt_depth
      module procedure opt_depth_s
      module procedure opt_depth_v
   end interface opt_depth

   interface BBintensity
      module procedure BBintensity_s
      module procedure BBintensity_v
   end interface BBintensity

   interface IC_emis_full
      module procedure IC_emis_full_s
      module procedure IC_emis_full_v
   end interface IC_emis_full

contains

   !============================================================================
   !>  Equation 4 in Marcowith & Malzac (2003)
   function CyclotronLimit(beta, m, nu_b) result(cyclo)
      implicit none
      integer, intent(in) :: m
      real(dp), intent(in) :: beta,nu_b
      real(dp) :: cyclo
#ifdef INTEL
      cyclo = 8d0 * pi**2 * nu_b * dble( (m + 1) * ( m**(2 * m + 1) ) ) * &
      beta**(2 * m) / tgamma(dble(2 * m + 2))
#else
      cyclo = 8d0 * pi**2 * nu_b * dble( (m + 1) * ( m**(2 * m + 1) ) ) * &
      beta**(2 * m) / dgamma(dble(2 * m + 2))
#endif
      if (cyclo.lt.1d-200) cyclo = 1d-200
   end function CyclotronLimit
   !============================================================================


   !============================================================================
   !> Eq. (3.40) from Relativistic Jets from Active Galactic Nuclei
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
#ifdef INTEL
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
      if ( (nu >= nu0 * g1**2) .and. (nu <= nu0* g2**2) ) then
         j_nu = 4d0 * cLight * (eCharge**2 / (mass_e * cLight**2))**2 * uB * &
               n0 * nu0**((q - 3d0) / 2d0) * nu**((1d0 - q) / 2d0) / 9d0
      else
         j_nu = 1d-200
      end if
   end function RJfAGN_eq340_emiss
   !============================================================================


   !============================================================================
   !> Equation (8) in Petrosian (1981)
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


   !> Trapezoidal integrator over viewing angles for P81
   subroutine P81_trapzd(theta_a, theta_b, beta, gam, chi, s, n)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: theta_a, theta_b, beta, gam, chi
      real(dp), intent(inout) :: s
      integer :: it, i
      real(dp) :: del, fsum, fa, fb, th
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
      integer, parameter :: JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1
      real(dp), parameter :: EPS=1d-12
      integer :: j
      real(dp) :: P81_qromb, dqromb
      real(dp), dimension(JMAXP) :: h, s
      h(1) = 1d0
      do j=1,JMAX
         call P81_trapzd(theta_a, theta_b, beta, gam, chi, s(j), j)
         if (j >= K) then
            call polint(h(j - KM:j), s(j - KM:j), 0d0, P81_qromb, dqromb)
            if (abs(dqromb) <= EPS * abs(P81_qromb)) return
         end if
         s(j+1) = s(j)
         h(j+1) = 0.25d0 * h(j)
      end do
      call an_error('P81_qromb: too many steps')
   end function P81_qromb
   !============================================================================


   !============================================================================
   !> The RMA function: RMAfit(x) = Rsync(x) = 0.5 pi x CS(x)
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


   function RMA_nobreak(chi, g) result(res)
      implicit none
      real(dp), intent(in) :: chi, g
      real(dp), parameter :: c1 = 3.2180900500625734d-4, &
            c2 = 6.50532122717873d-1, c3 = 1.5579904689804556d1
      real(dp) :: res, x
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
   end function RMA_nobreak


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
      ! emiss = pi * eCharge**2 * nu_b * n0 * I2 * gmin**qq / (2d0 * cLight)
      emiss = jmbconst * nu_b * n0 * I2 * gmin**qq
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


   subroutine syn_emissivity(jnu, freq, gg, nn, B)
      implicit none
      real(dp), intent(in) :: freq, B
      real(dp), intent(in), dimension(:) :: gg, nn
      real(dp), intent(out) :: jnu
      integer :: k, Ng
      real(dp) :: qq
      Ng = size(gg, dim=1)
      jnu = 0d0
      calc_jnu: do k = 1, Ng - 1
         if ( nn(k) > 1d-100 .and. nn(k + 1) > 1d-100) then
            qq = -dlog(nn(k + 1) / nn(k)) / dlog(gg(k + 1) / gg(k))
            if ( qq > 8d0 ) qq = 8d0
            if ( qq < -8d0 ) qq = -8d0
            !!!TODO: RMA_new argument is hardcoded. Change this
            jnu = jnu + j_mb(freq, B, nn(k), gg(k), gg(k + 1), qq, RMA_new)
         end if
      end do calc_jnu
      if ( jnu < 1d-200 ) jnu = 0d0
   end subroutine syn_emissivity


   !> Synchrotron emissivity following eq. (22) in Ryan G., van Eerten H., Piro L., Troja E., 2020, ApJ, 896, 166
   subroutine syn_broken(nu, tc, Gbulk, gm, p, B, n0, emiss)
      implicit none
      real(dp), intent(in) :: nu, p, Gbulk, gm, n0, tc, B
      real(dp), intent(out) :: emiss
      real(dp) :: num, nuc, epsP, xiN, gc, n
      xiN = 1d0
      n = 4d0 * n0 * Gbulk
      gc = 6d0 * pi * Gbulk * mass_e * cLight / (sigmaT * B**2 * tc)
      num = 3d0 * eCharge * B * gm**2 / (fourpi * mass_e * cLight)
      nuc = 3d0 * eCharge * B * gc**2 / (fourpi * mass_e * cLight)
      if ( nu < num .and. num < nuc ) then
         emiss = (nu / num)**(1d0 / 3d0)
      else if ( num <= nu .and. nu <= nuc ) then
         emiss = (nu / num)**(-0.5d0 * (p - 1d0))
      else if ( num < nuc .and. nuc < nu ) then
         emiss = (nuc / num)**(-0.5d0 * (p - 1d0)) * (nu / nuc)**(-0.5d0 * p)
      else if ( nu < nuc .and. nuc < num ) then
         emiss = (nu / nuc)**(1d0 / 3d0)
      else if ( nuc <= nu .and. nu <= num ) then
         emiss = (nu / nuc)**(-0.5d0)
      else
         emiss = (num / nuc)**(-0.5d0) * (nu / num)**(-0.5d0 * p)
      end if
      epsP = 0.5d0 * (p - 1d0) * sqrt3 * eCharge**3 * xiN * n * B / energy_e
      emiss = emiss * epsP
   end subroutine syn_broken


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
    !  absor = pi * eCharge**2 * nu_b * n0 * A2 * gmin**qq / (4d0 * mass_e * cLight * nu**2)
      absor = ambconst * nu_b * n0 * A2 * gmin**qq / nu**2
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


   subroutine syn_absorption(anu, freq, gg, nn, B)
      implicit none
      real(dp), intent(in) :: freq, B
      real(dp), intent(in), dimension(:) :: gg, nn
      real(dp), intent(out) :: anu
      integer :: k, Ng
      real(dp) :: qq
      Ng = size(gg, dim=1)
      anu = 0d0
      calc_anu: do k = 1, Ng - 1
         if ( nn(k) > 1d-100 .and. nn(k + 1) > 1d-100) then
            qq = -dlog(nn(k + 1) / nn(k)) / dlog(gg(k + 1) / gg(k))
            if ( qq > 8d0 ) qq = 8d0
            if ( qq < -8d0 ) qq = -8d0
            !!!TODO: RMA_new argument is hardcoded. Change this
            anu = anu + a_mb(freq, B, nn(k), gg(k), gg(k + 1), qq, RMA_new)
         end if
      end do calc_anu
      if ( anu < 1d-200 ) anu = 0d0
   end subroutine syn_absorption

  !!todo comment
   subroutine gg_absorption_trapzd(nu1, qq, lnu_a, lnu_b, s, n, func)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: nu1, qq, lnu_a, lnu_b
      real(dp), intent(inout) :: s
      interface
         function func(nu1_t, nu_t, qq_t)
            use data_types
            real(dp), intent(in) :: nu1_t, nu_t, qq_t
            real(dp) :: func
         end function func
      end interface
      integer :: it, i
      real(dp) :: del, fsum, lnu, fa, fb, enu_a, enu_b, enu
      if (n == 1) then
         enu_a = dexp(lnu_a)
         enu_b = dexp(lnu_b)
         fa = func(nu1,enu_a,qq)
         fb = func(nu1,enu_b,qq)
         s = 0.5d0 * (lnu_b - lnu_a) * (fa + fb)
      else
         it = 2**(n - 2)
         del = (lnu_b - lnu_a) / dble(it)
         lnu = lnu_a + 0.5d0 * del
         enu = dexp(lnu)
         fsum = 0d0
         itloop: do i=1,it
            fsum = fsum + func(nu1,enu,qq) !eg**(-q) * RMAfunc(chi, eg) * (q + 1d0 + eg**2 / (eg**2 - 1d0))
            lnu = lnu + del
            enu = dexp(lnu)
            ! write(*,*) fsum
         end do itloop
         s = 0.5d0 * (s + del * fsum) ! where del = (lgb - lga) / it
      end if
   end subroutine gg_absorption_trapzd

   !todo comment
   function gg_absorption_qromb(lnu_a, lnu_b, func,nu1, qq) result(qromb)
      implicit none
      real(dp), intent(in) :: lnu_a, lnu_b, nu1, qq
      interface
         function func(nu1_t, nu_t, qq_t)
            use data_types
            real(dp), intent(in) :: nu1_t, nu_t, qq_t
            real(dp) :: func
         end function func
      end interface
      integer, parameter :: JMAX = 30, JMAXP = JMAX + 1, K = 10, KM = K - 1
      real(dp), parameter :: EPS = 1d-3
      real(dp), dimension(JMAXP) :: h, s
      real(dp) :: dqromb, qromb
      integer :: j
      h(1) = 1d0
      do j = 1, JMAX
         call gg_absorption_trapzd(nu1, qq, lnu_a, lnu_b, s(j), j, func)
         if (j >= K) then
            call polint(h(j-KM:j), s(j-KM:j), 0d0, qromb, dqromb)
            if (dabs(dqromb) .le. EPS * dabs(qromb)) return
         end if
         s(j + 1) = s(j)
         h(j + 1) = 0.25d0 * h(j)
      end do
      print*,'gg_absorption_qromb error'
      print*,'nu1    =', nu1
      print*,'qq     =', qq
      print*,'nu_a     =', dexp(lnu_a)
      print*,'nu_b     =', dexp(lnu_b)
      print*,'qromb  =', qromb
      print*,'dqromb =', dqromb
      call an_error('gg_absorption_qromb: too many steps')
   end function gg_absorption_qromb

   !!calculates 10.9 from dermer and menon 2009
   !!eps is the convergence value for teh series
   function phi_bar(s0_d, eps) result(res)
     implicit none
     real(dp), intent(in) :: s0_d, eps
     real(dp) :: res, beta_0, w_0,a,b,c,d,sum,sum_temp,s0
     integer :: k,Nmax

     if(s0_d<1.0001d0) then
        s0=1.0001d0
     else
       s0 = s0_d
     end if
     ! write(*,*) s0_d
     if(s0_d < 1d0) then
       write(*,*) s0_d
     end if
     if(s0 > 100) then
       s0 = 100
     end if
     Nmax = 10 !usually converges ~10
     beta_0 = dsqrt(1d0 - (s0**(-1d0)))
     ! if(beta_0 > 1d0 - 1d-20) then
     !   beta_0 = 1d0 - 1d-20
     ! end if
     w_0 = (1d0 + beta_0)/(1d0 - beta_0)
     !
     ! if(w_0 > 1d100) then
     !   w_0 = 1d100
     ! end if
     a = ((1d0 + (beta_0**2d0))/(1d0 - (beta_0**2d0)))*LN1(w_0,1d-6)
     b = -(beta_0**2d0)*LN1(w_0,1d-6) - LN2(w_0,1d-6) - (4d0 * beta_0/(1d0 - (beta_0**2d0)))
     c = (2d0 * beta_0) + (4d0 * LN1(w_0,1d-6) * LN1(w_0 + 1d0,1d-6))
     sum = 0d0
     do k = 1, Nmax
       sum_temp = sum
       sum = sum + ((-1d0)**dble(k-1))*(dble(k)**(-2d0))*(w_0**(-dble(k)))
       if(isnan(sum))then
          write(*,*) 'sum: ',sum
          write(*,*) 'sumtemp', sum_temp
          write(*,*) 'w0', w_0
       end if
       ! write(*,*) 'entering loop3' ,abs(sum-sum_temp), w_0, beta_0, s0
       if(abs(sum-sum_temp) < eps) then
          ! write(*,*) 'exit loop3'
          exit
       end if
     end do
     d = -(4d0)*( (0.5d0)*LN2(w_0,1d-6)  + (pi**2d0)/(12d0) - sum)
     res = a + b + c + d
     ! write(*,*) 'a: ',a,' b: ',b, ' c: ', c ,' d:', d
     if(isnan(a)) write(*,*) 'a'
     if(isnan(b)) write(*,*) 'b'
     if(isnan(c)) write(*,*) 'c'
     if(isnan(d)) write(*,*) 'd: ',d
     if(isnan(sum)) write(*,*) 'sum: ',sum
     if(isnan(w_0)) write(*,*) 'w0: ',w_0


   end function phi_bar

   function gg_integration_function(nu1,nu,qq) result(res)
     implicit none
     real(dp), intent(in) :: nu1, nu,qq
     real(dp) :: res,phi_b,s0

     s0 = ((h_mec2)**2d0)*nu1*nu
     phi_b = phi_bar(s0,1d-1)
     if(s0<1d0) then
       write(*,*) s0,h_mec2, nu1, nu
     end if
     res = (nu**(-(qq + 3d0)))*phi_b
     ! if(isnan(phi_b)) write(*,*) 'phi_b'
     ! if(isnan(res)) write(*,*) 'res_gg_integration_function'
   end function gg_integration_function

   !!todo comment
   subroutine gg_absorption_i(nu_a,nu_b,anu_i,Inu_a,Inu_b,nu1)
      implicit none
      real(dp), intent(in) :: nu_a, nu_b,Inu_a,Inu_b,nu1
      real(dp), intent(out) :: anu_i
      real(dp), dimension(10):: x,f
      real(dp) :: qq,nu_min
      integer :: i

      ! nu_min = min(nu_a,(mec2_h**2d0)/nu1)
      nu_min = nu_a
      anu_i = 0d0
      if ( Inu_a > 1d-100 .and. Inu_b > 1d-100 .and. nu_b>nu_min) then
        qq = -dlog(Inu_b / Inu_a) / dlog(nu_b /nu_a)
        if ( qq > 8d0 ) qq = 8d0
        if ( qq < -8d0 ) qq = -8d0
        ! build_x: do i=1,10
        !   x(i) = nu_min * (nu_b / nu_min)**(dble(i - 1) / dble(10 - 1))
        !   ! write(*,*) x(i)
        !   f(i) = gg_integration_function(nu1,x(i),qq)
        ! end do build_x
        anu_i = gg_abs_con*(Inu_a/(nu1))*(nu_a**(qq))*gg_absorption_qromb(dlog(nu_a), dlog(nu_b), gg_integration_function,nu1, qq)!qromb_arr(x, nu_min, nu_b, f)
      end if
      if ( anu_i < 1d-200 ) anu_i = 0d0
   end subroutine gg_absorption_i

   subroutine gg_absorption(nu,anu,Inu,nu1)
     implicit none
     real(dp), intent(in) :: nu1
     real(dp), intent(in), dimension(:) :: Inu, nu
     real(dp), intent(out) :: anu
     integer :: k, Nf
     real(dp) :: qq, anu_tmp
     Nf = size(nu, dim=1)
     anu = 0d0
     calc_anu: do k = 1, Nf - 1
       ! write(*,*) 'loop 2'
       if(nu1*nu(k) > 1/(h_mec2**2d0) .and. nu1*nu(k+1) > 1/(h_mec2**2d0)) then
         ! write(*,*) 'loop 2'
         call gg_absorption_i(nu(k),nu(k+1),anu_tmp,Inu(k),Inu(k+1),nu1)
         anu = anu + anu_tmp
       end if
     end do calc_anu
     if ( anu < 1d-200 ) anu = 0d0
   end subroutine gg_absorption

   !============================================================================


   !============================================================================
   !  #####  #        ##    ####  #    #    #####   ####  #####  #   #
   !  #    # #       #  #  #    # #   #     #    # #    # #    #  # #
   !  #####  #      #    # #      ####      #####  #    # #    #   #
   !  #    # #      ###### #      #  #      #    # #    # #    #   #
   !  #    # #      #    # #    # #   #     #    # #    # #    #   #
   !  #####  ###### #    #  ####  #    #    #####   ####  #####    #
   function BBintensity_s(nu, T) result(B)
      implicit none
      real(dp), intent(in) :: nu, T
      real(dp) :: B
      B = 2d0 * hPlanck * nu**3 / (cLight**2 * (dexp(hPlanck * nu / (kBoltz * T)) - 1d0))
   end function BBintensity_s


   function BBintensity_v(nu, T) result(B)
      implicit none
      real(dp), intent(in) :: T
      real(dp), intent(in), dimension(:) :: nu
      real(dp), dimension(size(nu)) :: B
      B = 2d0 * hPlanck * nu**3 / (cLight**2 * (dexp(hPlanck * nu / (kBoltz * T)) - 1d0))
   end function BBintensity_v


   function BBenergy_dens(T) result(u)
      implicit none
      real(dp), intent(in) :: T
      real(dp) :: u
      u = 4d0 * sigmaSB * T**4 / cLight
   end function BBenergy_dens
   !============================================================================


   !   ####  #####  #####        #####  ###### #####  ##### #    #
   !  #    # #    #   #          #    # #      #    #   #   #    #
   !  #    # #    #   #          #    # #####  #    #   #   ######
   !  #    # #####    #   ###    #    # #      #####    #   #    #
   !  #    # #        #   ###    #    # #      #        #   #    #
   !   ####  #        #   ###    #####  ###### #        #   #    #
   function opt_depth_s(absor, s) result(tau)
      implicit none
      real(dp), intent(in) :: s
      real(dp), intent(in), dimension(:) :: absor
      real(dp), dimension(size(absor)) :: tau
      tau = absor * s
   end function opt_depth_s

   function opt_depth_v(absor, s) result(tau)
      implicit none
      real(dp), intent(in), dimension(:) :: s
      real(dp), intent(in), dimension(:, :) :: absor
      integer :: Ns, j, Nf
      real(dp), dimension(size(absor, dim=1)) :: tau
      Ns = size(s, dim=1)
      Nf = size(absor, dim=1)
      if ( Ns == 1 ) then
         tau = absor(:, 1) * s(1)
      else
         do j = 1, Nf
            tau = 0.5d0 * sum( (absor(j, :Ns - 1) + absor(j, 2:)) * (s(2:) - s(:Ns - 1)) )
         end do
      end if
   end function opt_depth_v

   function opt_depth_blob(absor, R) result(u)
      implicit none
      real(dp), intent(in) :: absor, R
      real(dp) :: u, tau
      tau = dmax1(1d-100, 2d0 * R * absor)
      if ( tau <= 1d-10 ) then
         u = 1d0
      else
         if ( tau > 100d0 ) then
            u = 0.5d0 - 1d0 / tau**2d0
         else if ( tau >= 0.01d0 .and. tau <= 100d0 ) then
            u = 0.5d0 * (1d0 - 2d0 * (1d0 - (1d0 + tau) * dexp(-tau)) / tau**2)
         else
            u = (tau / 3d0) - 0.125d0 * tau**2d0
         end if
         u = 3d0 * u / tau
      end if
   end function opt_depth_blob

   function opt_depth_slab(absor, r) result(u)
      implicit none
      real(dp), intent(in) :: absor, r
      real(dp) :: u, tau
      tau = dmax1(1d-100, r * absor)
      if ( tau <= 1d-10 ) then
         u = 1d0
      else
         u = (1d0 - dexp(-tau)) / tau
      end if
   end function opt_depth_slab


   ! # #    # ##### ###### #    #  ####  # ##### #   #
   ! # ##   #   #   #      ##   # #      #   #    # #
   ! # # #  #   #   #####  # #  #  ####  #   #     #
   ! # #  # #   #   #      #  # #      # #   #     #
   ! # #   ##   #   #      #   ## #    # #   #     #
   ! # #    #   #   ###### #    #  ####  #   #     #
   subroutine RadTrans_v(Inu, s, jnu, anu)
      ! Description:
      !   This function solves the radiative transfer equation.
      !
      implicit none
      real(dp), intent(in) :: s
      real(dp), intent(in), dimension(:) :: jnu, anu
      real(dp), intent(out), dimension(:) :: Inu
      integer :: j, Nf
      real(dp), dimension(size(jnu)) :: tau
      Nf = size(jnu)
      tau = opt_depth(anu, s)
      do j = 1, Nf
         if ( jnu(j) > 1d-100 ) then
            if ( tau(j) > 1d-10 ) then
               Inu(j) = s * jnu(j) * (1d0 - dexp(-tau(j))) / tau(j)
            else
               Inu(j) = s * jnu(j)
            end if
         else
            Inu(j) = 0d0
         end if
      end do
   end subroutine RadTrans_v


   subroutine RadTrans_s(Inu, s, jnu, anu)
      ! Description:
      !   This function solves the radiative transfer equation.
      implicit none
      real(dp), intent(in) :: s, jnu, anu
      real(dp), intent(out) :: Inu
      real(dp) :: tau
      tau = anu * s
      if ( jnu > 1d-100 ) then
         if ( tau > 1d-10 ) then
            Inu = s * jnu * (1d0 - dexp(-tau)) / tau
         else
            Inu = s * jnu
         end if
      else
         Inu = 0d0
      end if
   end subroutine RadTrans_s


   subroutine RadTrans_blob(Inu, s, jnu, anu)
      ! Description:
      !   This function solves the radiative transfer equation.
      !
      implicit none
      real(dp), intent(in) :: s
      real(dp), intent(in), dimension(:) :: jnu, anu
      real(dp), intent(out), dimension(:) :: Inu
      integer :: j, Nf
      Nf = size(jnu, dim=1)
      do j = 1, Nf
         Inu(j) = 2d0 * s * jnu(j) * opt_depth_blob(anu(j), s)
      end do
   end subroutine RadTrans_blob


   ! subroutine RT_line_of_sight
   !    implicit none
   !    integer :: i, ii, i_edge, i_start
   !    real(dp) :: abu, s_max, s_min
   !
   !    ! --> Light path from origin to the observer: 2 s mu_com
   !    ! --> Edge of the blob: 2 R mu_com
   !    ! --> Position at which we will measure the radiation:
   !    i_edge = minloc(R - s, dim = 1, mask = R - s >= 0d0)
   !
   !    !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED)&
   !    !$OMP& PRIVATE(s_min, s_max, abu, i, ii, i_start)
   !    tobs_loop: do i = 1, numdt
   !       if ( i <= i_edge ) then
   !          i_start = 1
   !       else
   !          i_start = i - i_edge
   !       end if
   !       tcom_loop: do ii = i_start, i
   !          if ( ii == 1 ) then
   !             s_min = x_com_f(0d0, t_obs(i), z, Gbulk, mu_obs) * 0.5d0 / mu_com
   !             s_max = x_com_f(t(1), t_obs(i), z, Gbulk, mu_obs) * 0.5d0 / mu_com
   !             abu = s_max - s_min
   !             call RadTrans(Iobs(:, i), abu, jnut(:, 1), anu=anut(:, 1))
   !          else
   !             s_min = x_com_f(t(ii - 1), t_obs(i), z, Gbulk, mu_obs) * 0.5d0 / mu_com
   !             s_max = x_com_f(t(ii), t_obs(i), z, Gbulk, mu_obs) * 0.5d0 / mu_com
   !             abu = s_max - s_min
   !             call RadTrans(Iobs(:, i), abu, jnut(:, ii), anu=anut(:, ii), I0=Iobs(:, i))
   !          end if
   !       end do tcom_loop
   !    end do tobs_loop
   !    !$OMP END PARALLEL DO
   ! end subroutine RT_line_of_sight


   subroutine bolometric_integ(freqs, uu, ubol)
      implicit none
      real(dp), intent(in), dimension(:) :: freqs, uu
      real(dp), intent(out) :: ubol
      integer :: j, Nf
      real(dp) :: uind
      Nf = size(freqs)
      ubol = 0d0
      freqloop: do j = 2, Nf
         if ( uu(j - 1) > 1d-200 .and. uu(j) > 1d-200 ) then
            uind = -dlog(uu(j) / uu(j - 1)) / dlog(freqs(j) / freqs(j - 1))
            if ( uind > 8d0 ) uind = 8d0
            if ( uind < -8d0 ) uind = -8d0
            ubol = ubol + uu(j - 1) * freqs(j - 1) * Pinteg(freqs(j) / freqs(j - 1), uind, 1d-9)
         end if
      end do freqloop
      ubol = dmax1(ubol, 1d-200)
   end subroutine bolometric_integ


   !  #######  #####    #    #####
   !       #  #     #  ##   #     #
   !      #   #       # #         #
   !     #     #####    #    #####
   !    #           #   #   #
   !   #      #     #   #   #
   !  #######  #####  ##### #######
   function ssccZS12(gg, nn, B, R) result(nu0)
      real(dp), intent(in) :: B, R
      real(dp), intent(in), dimension(:) :: nn, gg
      real(dp), dimension(size(gg)) :: nu0
      integer :: k,Ng
      real(dp) :: R15,D0,A0,Ig2n,pp
      Ng = size(gg)
      R15 = R * 1d-15
      D0 = 1.29d-9 * B**2
      A0 = 1.15d-18 * R15 * B**2
      Ig2n = 0d0
      do k = 1, Ng - 1
         if ( nn(k) > 1d-100 .and. nn(k + 1) > 1d-100 ) then
            pp = -dlog(nn(k + 1) * gg(k + 1)**2 / (nn(k) * gg(k)**2)) / dlog(gg(k + 1) / gg(k))
            if ( pp > 8d0 ) pp = 8d0
            if ( pp < -8d0 ) pp = -8d0
            Ig2n = Ig2n + nn(k) * gg(k)**3d0 * Pinteg(gg(k + 1) / gg(k), pp, 1d-9)
         endif
      enddo
      nu0 = D0 + A0 * Ig2n
   end function ssccZS12


   !  ###  #####
   !   #  #     #  ####  #    # #####  #####  ####  #    #
   !   #  #       #    # ##  ## #    #   #   #    # ##   #
   !   #  #       #    # # ## # #    #   #   #    # # #  #
   !   #  #       #    # #    # #####    #   #    # #  # #
   !   #  #     # #    # #    # #        #   #    # #   ##
   !  ###  #####   ####  #    # #        #    ####  #    #

   !
   !   ----------{   Integral over incomming frequencies   }----------
   !
   subroutine IC_emis_full_s(fout, fin, gg, nn, Inu, emiss)
      implicit none
      real(dp), intent(in) :: fout, fin, Inu
      real(dp), intent(in), dimension(:) :: gg, nn
      real(dp), intent(out) :: emiss
      integer :: Ng
      real(dp) :: g1, g2, gmin, gmax, I0
      Ng = size(gg)
      gmin = gg(1)
      gmax = gg(Ng)
      g1 = dmax1(dsqrt(0.25d0 * fout / fin), gmin)
      g2 = dmin1(3d0 * mass_e * cLight**2 / (4d0 * hPlanck * fin), gmax)
      if ( g1 >= g2 ) then
         I0 = 0d0
      else
         I0 = Inu * IC_qromb(fin, fout, g1, g2, gg, nn) / fin**2
      end if
      emiss = 0.25d0 * fout * sigmaT * I0
   end subroutine IC_emis_full_s

   subroutine IC_emis_full_v(fout, fin, gg, nn, Inu, emiss)
      implicit none
      real(dp), intent(in) :: fout
      real(dp), intent(in), dimension(:) :: gg, nn, fin, Inu
      real(dp), intent(out) :: emiss
      integer :: j, Nf, Ng
      real(dp) :: g1, g2, jbol0, gmin, gmax
      real(dp), dimension(size(fin)) :: I0
      Nf = size(fin)
      Ng = size(gg)
      gmin = gg(1)
      gmax = gg(Ng)
      fin_loop: do j = 1, Nf
         g1 = dmax1(dsqrt(0.25d0 * fout / fin(j)), gmin)
         ! g2 = dmin1(3d0 * mass_e * cLight**2 / (4d0 * hPlanck * fin(j)), gmax)
         g2 = dmin1(mass_e * cLight**2 / (hPlanck * fin(j)), gmax)
         if ( g1 >= g2 ) then
            I0(j) = 0d0
         else
            I0(j) = Inu(j) * IC_qromb(fin(j), fout, g1, g2, gg, nn) / fin(j)**2
         end if
      end do fin_loop
      !  :::::  improvised trapezoidal rule  :::::
      jbol0 = fin(1) * I0(1) + fin(Nf) * I0(Nf)
      jbol0 = 0.5d0 * dlog(fin(Nf) / fin(1)) * &
         (jbol0 + 2d0 * sum(fin(2:Nf - 1) * I0(2:Nf - 1))) / dble(Nf - 1)
      emiss = 0.25d0 * fout * sigmaT * jbol0
   end subroutine IC_emis_full_v


   !
   !   ----------{   Romberg integrator   }----------
   !
   function IC_qromb(fin, fout, a, b, gg, nn) result(qromb)
      implicit none
      real(dp), intent(in) :: fin, fout, a, b
      real(dp), intent(in), dimension(:) :: nn, gg
      integer, parameter :: jmax = 50, jmaxp = jmax + 1, kq = 10, km = kq - 1
      real(dp), parameter :: eps = 1d-5
      integer :: jq
      real(dp) :: dqromb, qromb
      real(dp), dimension(jmaxp) :: h, s
      h(1) = 1d0
      do jq = 1, jmax
         call IC_trapzd(fin, fout, dlog(a), dlog(b), s(jq), jq, gg, nn)
         if ( jq >= kq ) then
            call polint(h(jq - km:jq), s(jq - km:jq), 0.0d0, qromb, dqromb)
            if ( dabs(dqromb) <= eps * dabs(qromb) ) return
         end if
         s(jq + 1) = s(jq)
         h(jq + 1) = 0.25d0 * h(jq)
      end do
      print*, fin, fout, a, b, qromb, dqromb
      call an_error('IC_qromb: too many steps')
   end function IC_qromb

   !
   !   ----------{   Trapezoid   }----------
   !
   subroutine IC_trapzd(fin, fout, a, b, s, n, gg, nn)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: a, b, fout, fin
      real(dp), intent(in), dimension(:) :: gg, nn
      real(dp), intent(inout) :: s
      integer :: it, jt
      real(dp) :: del, fsum, x
      if ( n == 1 ) then
         s = 0.5d0 * (b - a) * ( &
         dexp(a) * IC_integrand(fin, fout, dexp(a), gg, nn) + &
         dexp(b) * IC_integrand(fin, fout, dexp(b), gg, nn) )
      else
         it = 2**(n - 2)
         del = (b - a) / dble(it)
         x = a + 0.5d0 * del
         fsum = 0d0
         do jt=1,it
            fsum = fsum + dexp(x) * IC_integrand(fin, fout, dexp(x), gg, nn)
            x = x + del
         end do
         s = 0.5d0 * (s + del * fsum)
      end if
   end subroutine IC_trapzd

   !
   !   ----------{   Dispersion function   }----------
   !
   function IC_integrand(fin, fout, gev, g, n) result(integrand)
      implicit none
      real(dp), intent(in) :: fin, fout, gev
      real(dp), intent(in), dimension(:) :: g, n
      integer :: k, Ng
      real(dp) :: fIC, beta, integrand, q, foutfin, bplus, bminus, x
      Ng = size(g)
      k = minloc(dabs(gev - g), dim=1)

      if ( n(k) < 1d-100 .or. k <= 1 ) then
         integrand = 0d0
         return
      end if

      beta = bofg(gev)
      foutfin = fout / fin
      bplus = 1d0 + beta
      bminus = 1d0 - beta
      x = 0.25 * foutfin / gev**2

      if ( bminus < 1e-5 ) then
         fIC = 3d0 * (2d0 * x * (dlog(x) - x) + x + 1d0)
      else
         if ( bminus / bplus <= foutfin .and. foutfin <= 1d0 ) then
            fIC = bplus * foutfin - bminus
         else if ( 1d0 < foutfin .and. foutfin <= bplus / bminus ) then
            fIC = bplus - bminus * foutfin
         else
            integrand = 0d0
            return
         end if
      end if

      if ( dabs(gev - g(k)) < 1d-15 ) then
         integrand = n(k) * fIC / (beta * gev)**2
      else
         if ( gev < g(k) .or. gev > g(Ng) ) then
            q = dlog(n(k) / n(k - 1)) / dlog(g(k) / g(k - 1))
            if ( q < -8d0 ) q = -8d0
            if ( q >  8d0 ) q =  8d0
            integrand = n(k - 1) * (gev / g(k - 1))**q * fIC / (beta * gev)**2
         else
            q = dlog(n(k + 1) / n(k)) / dlog(g(k + 1) / g(k))
            if ( q < -8d0 ) q = -8d0
            if ( q >  8d0 ) q =  8d0
            integrand = n(k) * (gev / g(k))**q * fIC / (beta * gev)**2
         end if
      end if
   end function IC_integrand


   !
   !  -----  Inverse Compton for isotropic power-law incomming photons  -----
   !
   subroutine IC_iso_powlaw(jnu, nuout, nu, Inu, n, g)
      implicit none
      real(dp), intent(in) :: nuout
      real(dp), intent(in), dimension(:) :: nu, n, g, Inu
      real(dp), intent(out) :: jnu
      integer :: j, k, Ng, Nf
      real(dp) :: w1, w2, gmx_star, gKN, l, q, f1, f2, emis, g1, g2, s1, s2
      Ng = size(g, dim=1)
      Nf = size(nu, dim=1)
      jnu = 0d0
      nu_loop: do j = 1, Nf - 1
         intens_if: if ( Inu(j) > 1d-200 .and. Inu(j + 1) > 1d-200 ) then
            l = -dlog(Inu(j + 1) / Inu(j)) / dlog(nu(j + 1) / nu(j))
            if ( l > 8d0 ) l = 8d0
            if ( l < -8d0 ) l = -8d0
            gKN = mass_e * cLight**2 / (hPlanck * nu(j + 1))
            f1 = nuout / (4d0 * nu(j))
            f2 = nuout / (4d0 * nu(j + 1))
            g2 = dmin1(g(Ng), gKN)
            g1 = dmax1(g(1), dsqrt(f1))
            g1g2_cond: if ( g1 < g2 ) then
               g_loop: do k = 1, Ng - 1
                  if ( g(k) < g1 ) cycle g_loop
                  if ( g(k) > g2 ) exit g_loop
                  e_dist: if ( n(k) > 1d-200 .and. n(k + 1) > 1d-200 ) then
                     q = -dlog(n(k + 1) / n(k)) / dlog(g(k + 1) / g(k))
                     if ( q > 8d0 ) q = 8d0
                     if ( q < -8d0 ) q = -8d0
                     s1 = (q - 1d0) / 2d0
                     s2 = (2d0 * l - q - 1d0) / 2d0
                     if ( s1 > 8d0 ) s1 = 8d0
                     if ( s1 < -8d0 ) s1 = -8d0
                     if ( s2 > 8d0 ) s2 = 8d0
                     if ( s2 < -8d0 ) s2 = -8d0
                     gmx_star = dmin1(g(k + 1), gKN)
                     w1 = dmin1(f1, gmx_star**2)
                     w2 = dmax1(f2, 0.25d0)
                     contrib_if: if ( w1 > w2 ) then
                        if ( 0.25d0 < f1 .and. f1 < g(k)**2 ) then
                           emis = sscG1ISO(gmx_star**(-2), g(k)**(-2), w2, w1, s1, s2)
                        else if ( f1 <= g(k)**2 .and. f2 <= g(k)**2 ) then
                           emis = sscG1ISO(gmx_star**(-2), g(k)**(-2), w2, g(k)**2, s1, s2) &
                                 + sscG2ISO(gmx_star**(-2), 1d0, g(k)**2, w1, s1, s2)
                        else if ( f2 <= gmx_star**2 .and. f2 > g(k)**2 ) then
                           emis = sscG2ISO(gmx_star**(-2), 1d0, w2, w1, s1, s2)
                        else
                           emis = 0d0
                        end if
                        jnu = jnu + emis * n(k) * (g(k)**q) * Inu(j) * sigmaT * (f1**(-l))
                     end if contrib_if
                  end if e_dist
               end do g_loop
            end if g1g2_cond
         end if intens_if
      end do nu_loop
   end subroutine IC_iso_powlaw


   !
   !  -----  Inverse Compton for isotropic power-law incomming photons  -----
   !
   subroutine IC_iso_monochrom(jnu, nuout, uext, nuext, n, g)
      implicit none
      real(dp), intent(in) :: uext, nuext, nuout
      real(dp), intent(in), dimension(:) :: n, g
      real(dp), intent(out) :: jnu
      real(dp), parameter :: eps = 1d-9
      integer :: k, Ng
      real(dp) :: w, gmx_star, gKN, q, q1, q2, emis
      Ng = size(g, dim=1)
      gKN = mass_e * cLight**2 / (hPlanck * nuext)
      w = nuout / (4d0 * nuext)
      jnu = 0d0
      emis = 0d0
      g_loop: do k = 1, Ng - 1
         gmx_star = dmin1(g(k + 1), gKN)
         e_dist: if ( n(k) > 1d-200 .and. n(k + 1) > 1d-200 ) then
            q = -dlog(n(k + 1) / n(k)) / dlog(g(k + 1) / g(k))
            if ( q > 8d0 ) q = 8d0
            if ( q < -8d0 ) q = -8d0
            q1 = 0.5d0 * (q - 1d0)
            q2 = 0.5d0 * (q + 1d0)
            if ( q1 > 8d0 ) q1 = 8d0
            if ( q1 < -8d0 ) q1 = -8d0
            if ( q2 > 8d0 ) q2 = 8d0
            if ( q2 < -8d0 ) q2 = -8d0
            contrib_if: if ( 0.25d0 <= w .and. w <= g(k)**2 .and. g(k) <= gmx_star ) then
               emis = (w / gmx_star**2)**q2 * ( Pinteg((gmx_star / g(k))**2, -q1, eps) &
                     - (w / gmx_star**2) * Pinteg((gmx_star / g(k))**2, -q2, eps) )
            else if ( g(k)**2 < w .and. w <= gmx_star**2 ) then
               emis = (w / gmx_star**2)**q2 * ( Pinteg(gmx_star**2 / w, -q1, eps) &
                     - (w / gmx_star**2) * Pinteg(gmx_star**2 / w, -q2, eps) )
            else
               emis = 0d0
            end if contrib_if
            jnu = jnu + emis * n(k) * g(k)**q * w**(-q1) * cLight * sigmaT * uext / (4d0 * nuext)
         end if e_dist
      end do g_loop
   end subroutine IC_iso_monochrom


   !   ####   ####   ####  #      # #    #  ####
   !  #    # #    # #    # #      # ##   # #    #
   !  #      #    # #    # #      # # #  # #
   !  #      #    # #    # #      # #  # # #  ###
   !  #    # #    # #    # #      # #   ## #    #
   !   ####   ####   ####  ###### # #    #  ####
   !
   !     ----------> Power-law radiation field <----------
   !
   subroutine rad_cool_pwl(dotg, gg, freqs, uu, withKN)
      implicit none
      real(dp), intent(in), dimension(:) :: freqs, gg, uu
      real(dp), intent(out), dimension(:) :: dotg
      logical, intent(in) :: withKN
      integer :: j, k, Ng, Nf
      real(dp) :: uind, urad_const, usum, xi_c, xi_rat, ueval
      real(dp), dimension(size(gg, dim=1), size(freqs, dim=1)) :: xi, uxi
      urad_const = 4d0 * sigmaT * cLight / (3d0 * energy_e)
      xi_c = 4d0 * h_mec2
      Ng = size(gg, dim=1)
      Nf = size(freqs, dim=1)
      do k = 1, Ng
         do j = 1, Nf
            xi(k, j) = xi_c * gg(k) * freqs(j)
            uxi(k, j) = uu(j) * xi_c * gg(k)
         end do
      end do
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(k,j,uind,xi_rat,usum)
      do k = 1, Ng
         usum = 0d0
         ueval = 0d0
         freqloop: do j = 1, Nf - 1
            if ( uxi(k, j + 1) > 1d-100 .and. uxi(k, j) > 1d-100 ) then
               xi_rat = xi(k, j + 1) / xi(k, j)
               uind = dlog(uxi(k, j + 1) / uxi(k, j)) / dlog(xi_rat)
               if ( uind > 8d0 ) uind = 8d0
               if ( uind < -8d0 ) uind = -8d0
               if (withKN) then
                  if (xi(k, j) >= 1d2) then
                     ueval = 4.5d0 * uxi(k, j) &
                           * ( Qinteg(xi_rat, uind + 2d0, 1d-6) &
                           + (dlog(xi(k, j)) - (11d0 / 6d0)) &
                           * Pinteg(xi_rat, uind + 2d0, 1d-6) ) / xi(k, j)
                  else if (xi(k, j) >= 1d0 .and. xi(k, j) < 1d2) then
                     ueval = qromb_w2arg(transKN_fit,xi(k, j), xi(k, j + 1), uind) * uxi(k, j) * xi(k, j)**uind
                  else if (xi(k, j) >= 1d-3 .and. xi(k, j) < 1d0) then
                     ueval = qromb_w2arg(transKN,xi(k, j), xi(k, j + 1), uind) * uxi(k, j) * xi(k, j)**uind
                  else
                     ueval = uxi(k, j) * xi(k, j) * Pinteg(xi_rat, uind, 1d-6)
                  end if
               else
                 ueval = uxi(k, j) * xi(k, j) * Pinteg(xi_rat, uind, 1d-6)
               end if
            end if
            usum = usum + ueval
         end do freqloop
         dotg(k) = urad_const * usum / xi_c**2
      end do
      !$OMP END PARALLEL DO

      contains

      function transKN(x, p) result(integrando)
         implicit none
         real(dp), intent(in) :: p
         real(dp), intent(in), dimension(:) :: x
         real(dp), dimension(size(x, dim=1)) :: integrando
         integrando = x**(-p) * (1d0 + x)**(-1.5d0)
      end function transKN

      function transKN_fit(x, p) result(integrando)
         implicit none
         real(dp), intent(in) :: p
         real(dp), intent(in), dimension(:) :: x
         real(dp), dimension(size(x, dim=1)) :: integrando
         integrando = x**(-p) * dexp(-1.01819432d0 &
               - 0.67980349d0 * dlog(x) &
               - 0.14948459d0 * dlog(x)**2 &
               + 0.00627589d0 * dlog(x)**3)
      end function transKN_fit

   end subroutine rad_cool_pwl

   !
   !     ----------> Monoenergetic radiation field <----------
   !
   subroutine rad_cool_mono(dotg, gg, nu0, u0, withKN)
      implicit none
      real(dp), intent(in) :: u0, nu0
      real(dp), intent(in), dimension(:) :: gg
      logical, intent(in) :: withKN
      real(dp), intent(out), dimension(:) :: dotg
      integer :: Ng,k
      real(dp) :: urad_const
      real(dp), dimension(size(gg, dim=1)) :: xi0
      Ng = size(gg, dim=1)
      urad_const = 4d0 * sigmaT * cLight / (3d0 * energy_e)
      xi0 = 4d0 * gg * nu0 * h_mec2
      do k=1,Ng
         if (withKN) then
            if (xi0(k) >= 1d2) then
               dotg(k) = urad_const * u0 * gg(k)**2 * 4.5d0 * (dlog(xi0(k)) - 11d0 / 6d0) / xi0(k)**2
            else if (xi0(k) >= 1d0 .and. xi0(k) < 1d2) then
               dotg(k) = urad_const * u0 * gg(k)**2 * dexp(-1.01819432d0 &
                     - 0.67980349d0 * LN1(xi0(k), 1d-6) &
                     - 0.14948459d0 * LN2(xi0(k), 1d-6) &
                     + 0.00627589d0 * LN3(xi0(k), 1d-6))
            else if (xi0(k) > 1d-3 .and. xi0(k) < 1d0) then
               dotg(k) = urad_const * u0 * pofg(gg(k))**2 * (1d0 + xi0(k))**(-1.5d0)
            else
               dotg(k) = urad_const * u0 * pofg(gg(k))**2
            end if
         else
            dotg(k) = urad_const * u0 * pofg(gg(k))**2
         end if
      end do
   end subroutine rad_cool_mono


   !  ###### #    #  ####  #      #    # ##### #  ####  #    #
   !  #      #    # #    # #      #    #   #   # #    # ##   #
   !  #####  #    # #    # #      #    #   #   # #    # # #  #
   !  #      #    # #    # #      #    #   #   # #    # #  # #
   !  #       #  #  #    # #      #    #   #   # #    # #   ##
   !  ######   ##    ####  ######  ####    #   #  ####  #    #
   subroutine photons_evol(dt, nin, nout, nu, QQ, tesc, Loss)
      implicit none
      real(dp), intent(in) :: dt, tesc
      real(dp), intent(in), dimension(:) :: nu, QQ, Loss, nin
      real(dp), intent(out), dimension(:) :: nout
      integer :: Nf
      real(dp), dimension(size(nu)) :: b, r, zero
      Nf = size(nu)
      zero = zeros1D(Nf,.true.)
      r = nin + dt * QQ
      b = 1d0 + dt * (Loss + 1d0 / tesc)
      call tridag_ser(zero(2:), b, zero(2:), r, nout)
   end subroutine photons_evol


   subroutine intens_evol(dt, Iin, Iout, nu, jnu, anu)
      implicit none
      real(dp), intent(in) :: dt
      real(dp), intent(in), dimension(:) :: nu, jnu, anu, Iin
      real(dp), intent(out), dimension(:) :: Iout
      integer :: Nf
      real(dp), dimension(size(nu)) :: b, r, zero
      Nf = size(nu)
      zero = zeros1D(Nf,.true.)
      b = 1d0 + dt * anu * cLight
      r = Iin + dt * jnu * cLight
      call tridag_ser(zero(2:), b, zero(2:), r, Iout)
   end subroutine intens_evol


   subroutine Kompaneets_FinDif(dt, nin, nout, nu, th_e, ne, tesc)
      implicit none
      real(dp), intent(in) :: dt, th_e, ne, tesc
      real(dp), intent(in), dimension(:) :: nu, nin!, QQ
      real(dp), intent(out), dimension(:) :: nout
      real(dp), parameter :: eps = 1e-3
      integer :: i, Ng, Ng1
      real(dp) :: dBB, dtc, fth_e
      real(dp), dimension(size(nu)) :: dx, dxp2, dxm2, CCp2, CCm2, BBp2, BBm2, YYp2, YYm2, WWp2, WWm2, ZZp2, ZZm2, a, b, c, r, x, DD

      Ng = size(nu)
      Ng1 = Ng - 1

      dtc = ne * sigmaT * cLight * dt
      fth_e = th_e * (1d0 + 3.683 * th_e + 4d0 * th_e**2) / (1d0 + th_e)
      x = hPlanck * nu / (th_e * mass_e * cLight**2)
      DD = x**4
      dxp2(:Ng1) = x(2:) - x(:Ng1)
      dxp2(Ng) = dxp2(Ng1)
      dxm2(2:) = dxp2(:Ng1)
      dxm2(1) = dxm2(2)
      dx = dsqrt(dxp2 * dxm2)

      CCp2(:Ng1) = 0.25d0 * (DD(2:) + DD(:Ng1))
      CCm2(2:) = CCp2(:Ng1)
      CCp2(Ng) = 0.25d0 * DD(Ng)
      CCm2(1) = 0.25d0 * DD(1)
      call polint(x(2:), CCm2(2:), x(1), CCm2(1), dBB)

      BBp2 = CCp2
      BBm2 = CCm2

      WWp2 = dxp2 * BBp2 / CCp2
      WWm2 = dxm2 * BBm2 / CCm2

      do i = 1, Ng

         if ( 0.5d0 * WWp2(i) > 200d0 ) then
            ZZp2(i) = 200d0
         else if ( 0.5d0 * WWp2(i) < -200d0 ) then
            ZZp2(i) = -200d0
         else
            ZZp2(i) = 0.5d0 * WWp2(i)
         end if

         if ( 0.5d0 * WWm2(i) > 200d0 ) then
            ZZm2(i) = 200d0
         else if ( 0.5d0 * WWm2(i) < -200d0 ) then
            ZZm2(i) = -200d0
         else
            ZZm2(i) = 0.5d0 * WWm2(i)
         end if

         ! if ( 127d0 * WWp2(i)**8 < eps * 154828800d0 ) then
         if ( dabs(WWp2(i)) < 0.1d0 ) then
            YYp2(i) = 1d0 - WWp2(i)**2 / 24d0 + 7d0 * WWp2(i)**4 / 5760d0 - 31d0 * WWp2(i)**6 / 967680d0
         else
            YYp2(i) = dabs(WWp2(i)) * dexp(-dabs(ZZp2(i))) / ( 1d0 - dexp(-2d0 * dabs(ZZp2(i))) )
            ! YYp2(i) = dabs(WWp2(i)) / ( (1d0 - 1d0 / ZZp2(i)**2 ) * ZZp2(i) )
         end if

         ! if ( 127d0 * WWm2(i)**8 < eps * 154828800d0 ) then
         if ( dabs(WWm2(i)) < 0.1d0 ) then
            YYm2(i) = 1d0 - WWm2(i)**2 / 24d0 + 7d0 * WWm2(i)**4 / 5760d0 - 31d0 * WWm2(i)**6 / 967680d0
         else
            YYm2(i) = dabs(WWm2(i)) / ( dexp(dabs(ZZm2(i))) - dexp(-dabs(ZZm2(i))) )
         end if

      end do

      r = nin !+ dtc * QQ
      c = -dtc * fth_e * CCp2 * YYp2 * dexp(ZZp2)/ (dx * dxp2 * x**2)
      a = -dtc * fth_e * CCm2 * YYm2 * dexp(-ZZm2) / (dx * dxm2 * x**2)
      b = 1d0 + dtc * fth_e * ( CCp2 * YYp2 * dexp(-ZZp2) / dxp2 + CCm2 * YYm2 * dexp(ZZm2) / dxm2 ) / (dx * x**2) + dt / tesc !* (Loss + 1d0 / tesc)

      call tridag_ser(a(2:), b, c(:Ng1), r, nout)

   end subroutine Kompaneets_FinDif


   ! ######                #######
   ! #     #   ##   #####  #       #      #    # #    #
   ! #     #  #  #  #    # #       #      #    #  #  #
   ! ######  #    # #    # #####   #      #    #   ##
   ! #   #   ###### #    # #       #      #    #   ##
   ! #    #  #    # #    # #       #      #    #  #  #
   ! #     # #    # #####  #       ######  ####  #    #
   subroutine observed_Fnu(r, doppler, jnut, z, dlum, flux)
      implicit none
      real(dp), intent(in) :: z, dlum
      real(dp), intent(in), dimension(:) :: r, doppler
      real(dp), intent(in), dimension(:,:) :: jnut
      real(dp), intent(out), dimension(:,:) :: flux
      integer :: i, j, numf, numr
      numf = size(jnut, dim=1)
      numr = size(jnut, dim=2)
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(i, j)
      do j=1, numf
         do i=2, numr
            flux(j, i) = dmax1(1d-200, (1d0 + z) * &
                  qromb_arr(r(:i), r(1), r(i), jnut(j, :i) * (r(:i) * doppler(:i))**2) &
                  / (fourpi * dlum**2))
         end do
      end do
      !$OMP END PARALLEL DO
   end subroutine observed_Fnu

   subroutine observed_nuFnu(r, doppler, jnut, nu, dlum, flux)
      implicit none
      real(dp), intent(in) :: dlum
      real(dp), intent(in), dimension(:) :: r, doppler
      real(dp), intent(in), dimension(:,:) :: jnut, nu
      integer :: i, j, numf, numr
      real(dp), dimension(size(jnut, dim=1), size(jnut, dim=2)) :: flux
      numf = size(jnut, dim=1)
      numr = size(jnut, dim=2)
      flux = 1d-200
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) PRIVATE(i, j)
      do j=1, numf
         do i=2, numr
            flux(j, i) = dmax1(1d-200, &
                  qromb_arr(r(:i), r(1), r(i), nu(j, :i) * jnut(j, :i) * r(:i)**2 * doppler(:i)**3) / (fourpi * dlum**2))
         end do
      end do
      !$OMP END PARALLEL DO
   end subroutine observed_nuFnu

   ! function observed_flux_struc(r, doppler, th, phi, jnut, dlum, z) result(flux)
   !    implicit none
   !    real(dp) :: r, doppler, th, phi, jnut, dlum, z
   !    real(dp) :: flux
   ! end function observed_flux_struc

end module radiation
