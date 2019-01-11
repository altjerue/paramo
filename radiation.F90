module radiation
   use data_types
   use constants
   use misc
   use pwl_integ
   use SRtoolkit
   use anaFormulae
   !$ use omp_lib
   implicit none

   interface RadTrans
      module procedure RadTrans_v
      module procedure RAdTrans_m
   end interface RadTrans

   interface opt_depth
      module procedure opt_depth_s
      module procedure opt_depth_v
   end interface opt_depth

contains
   ! #    # #####   ####     ###### #    # #  ####   ####
   ! ##  ## #    # #         #      ##  ## # #      #
   ! # ## # #####   ####     #####  # ## # #  ####   ####
   ! #    # #    #      #    #      #    # #      #      #
   ! #    # #    # #    #    #      #    # # #    # #    #
   ! #    # #####   ####     ###### #    # #  ####   ####
   subroutine mbs_emissivity(jnu, freqs, gg, nn, B)
      implicit none
      real(dp), intent(in) :: B
      real(dp), intent(in), dimension(:) :: freqs, gg, nn
      real(dp), intent(out), dimension(:) :: jnu
      integer :: j, k, Ng, Nf
      real(dp) :: qq
      Ng = size(gg, dim=1)
      Nf = size(freqs, dim=1)
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(qq, j, k)
      freqs_loop: do j = 1, Nf
         jnu(j) = 0d0
         calc_jnu: do k = 1, Ng - 1
            if ( nn(k) > 1d-100 .and. nn(k + 1) > 1d-100) then
               qq = -dlog(nn(k + 1) / nn(k)) / dlog(gg(k + 1) / gg(k))
               if ( qq > 8d0 ) qq = 8d0
               if ( qq < -8d0 ) qq = -8d0
               jnu(j) = jnu(j) + j_mb(freqs(j), B, nn(k), gg(k), gg(k + 1), qq, RMA_new)
            end if
         end do calc_jnu
         if ( jnu(j) < 1d-200 ) jnu(j) = 0d0
      end do freqs_loop
      !$OMP END PARALLEL DO
   end subroutine mbs_emissivity


   ! #    # #####   ####       ##   #####   ####   ####  #####
   ! ##  ## #    # #          #  #  #    # #      #    # #    #
   ! # ## # #####   ####     #    # #####   ####  #    # #    #
   ! #    # #    #      #    ###### #    #      # #    # #####
   ! #    # #    # #    #    #    # #    # #    # #    # #   #
   ! #    # #####   ####     #    # #####   ####   ####  #    #
   subroutine mbs_absorption(anu, freqs, gg, nn, B)
      implicit none
      real(dp), intent(in) :: B
      real(dp), intent(in), dimension(:) :: freqs, gg, nn
      real(dp), intent(out), dimension(:) :: anu
      integer :: j, k, Ng, Nf
      real(dp) :: qq
      Ng = size(gg, dim=1)
      Nf = size(freqs, dim=1)
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(qq, j, k)
      freqs_loop: do j = 1, Nf
         anu(j) = 0d0
         calc_anu: do k = 1, Ng - 1
            if ( nn(k) > 1d-100 .and. nn(k + 1) > 1d-100) then
               qq = -dlog(nn(k + 1) / nn(k)) / dlog(gg(k + 1) / gg(k))
               if ( qq > 8d0 ) qq = 8d0
               if ( qq < -8d0 ) qq = -8d0
               anu(j) = anu(j) + a_mb(freqs(j), B, nn(k), gg(k), gg(k + 1), qq, RMA_new)
            end if
         end do calc_anu
         if ( anu(j) < 1d-200 ) anu(j) = 0d0
      end do freqs_loop
      !$OMP END PARALLEL DO
   end subroutine mbs_absorption


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
      real(dp), dimension(size(absor, dim=1)) :: tau
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

   function opt_depth_blob(tau) result(u)
      implicit none
      real(dp), intent(in) :: tau
      real(dp) :: u
      if ( tau > 100d0 ) then
         u = 0.5d0 - 1d0 / tau**2
      else if ( tau >= 0.01 .and. tau <= 100.0 ) then
         u = 0.5d0 * (1d0 - 2d0 * (1d0 - (1d0 + tau) * dexp(-tau)) / tau**2)
      else
         u = (tau / 3d0) - 0.125d0 * tau**2
      end if
   end function opt_depth_blob



   ! # #    # ##### ###### #    #  ####  # ##### #   #
   ! # ##   #   #   #      ##   # #      #   #    # #
   ! # # #  #   #   #####  # #  #  ####  #   #     #
   ! # #  # #   #   #      #  # #      # #   #     #
   ! # #   ##   #   #      #   ## #    # #   #     #
   ! # #    #   #   ###### #    #  ####  #   #     #
   subroutine RadTrans_v(Inu, s, jnu, anu, I0)
      ! Description:
      !   This function solves the radiative transfer equation.
      !
      implicit none
      real(dp), intent(in) :: s
      real(dp), intent(in), dimension(:) :: I0, jnu, anu
      real(dp), intent(out), dimension(:) :: Inu
      optional :: anu, I0
      integer :: j, Nf
      real(dp) :: Snu
      real(dp), dimension(size(jnu, dim=1)) :: tau
      Nf = size(jnu, dim=1)
      if ( present(I0) ) then
         Inu = I0
      else
         Inu = 0d0
      end if
      if ( present(anu) ) then
         tau = opt_depth(anu, s)
      else
         tau = 0d0
      end if
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j, Snu)
      do j = 1, Nf
         if ( jnu(j) > 1d-100 ) then
            if ( tau(j) >= 1d0 ) then
               Snu = jnu(j) / anu(j)
               Inu(j) = Snu + dexp(-tau(j)) * (Inu(j) - Snu)
            else
               Inu(j) = Inu(j) + s * jnu(j)
            end if
         else
            if ( tau(j) >= 1d0 ) then
               Inu(j) = Inu(j) * dexp(-tau(j))
            else
               Inu(j) = Inu(j)
            end if
         end if
      end do
      !$OMP END PARALLEL DO
   end subroutine RadTrans_v


   subroutine RadTrans_m(Inu, s, jnu, anu, I0)
      ! Description:
      !   This function solves the radiative transfer equation.
      !
      implicit none
      real(dp), intent(in), dimension(:) :: s, I0
      real(dp), intent(in), dimension(:, :) :: jnu, anu
      real(dp), intent(out), dimension(:) :: Inu
      optional :: anu, I0
      integer :: j, Nf, Ns, i
      real(dp) :: Snu
      real(dp), dimension(size(jnu, dim=1)) :: tau
      Nf = size(jnu, dim=1)
      Ns = size(s, dim=1)
      if ( present(I0) ) then
         Inu = I0
      else
         Inu = 0d0
      end if
      if ( present(anu) ) then
         tau = opt_depth(anu, s)
      else
         tau = 0d0
      end if
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(i, j, Snu)
      do j = 1, Nf
         do i = 1, Ns
            if ( jnu(j, i) > 1d-100 ) then
               if ( tau(j) >= 1d0 ) then
                  Snu = jnu(j, i) / anu(j, i)
                  Inu(j) = Snu + dexp(-tau(j)) * (Inu(j) - Snu)
               else
                  Inu(j) = Inu(j) + s(i) * jnu(j, i)
               end if
            else
               if ( tau(j) >= 1d0 ) then
                  Inu(j) = Inu(j) * dexp(-tau(j))
               else
                  Inu(j) = Inu(j)
               end if
            end if
         end do
      end do
      !$OMP END PARALLEL DO
   end subroutine RadTrans_m


   subroutine RadTrans_blob(Inu, s, jnu, anu)
      implicit none
      real(dp), intent(in) :: s
      real(dp), intent(in), dimension(:) :: jnu, anu
      real(dp), intent(out), dimension(:) :: Inu
      integer :: j, Nf
      real(dp) :: Snu
      real(dp), dimension(size(jnu, dim=1)) :: tau
      Nf = size(jnu, dim=1)
      tau = opt_depth(anu, s)
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j, Snu)
      do j = 1, Nf
         if ( jnu(j) > 1d-100 ) then
            if ( tau(j) > 1d-100 ) then
               Snu = jnu(j) / anu(j)
               Inu(j) = 0.25d0 * opt_depth_blob(tau(j)) * jnu(j) * s / (pi * tau(j))
            else
               Inu(j) = 0.25d0 * s * jnu(j) / pi
            end if
         end if
      end do
      !$OMP END PARALLEL DO
   end subroutine RadTrans_blob


   !  ###  #####
   !   #  #     #     ####   ####   ####  #      # #    #  ####
   !   #  #          #    # #    # #    # #      # ##   # #    #
   !   #  #          #      #    # #    # #      # # #  # #
   !   #  #          #      #    # #    # #      # #  # # #  ###
   !   #  #     #    #    # #    # #    # #      # #   ## #    #
   !  ###  #####      ####   ####   ####  ###### # #    #  ####
   function IC_cool(gg, freqs, Inu) result(urad)
      implicit none
      real(dp), intent(in), dimension(:) :: freqs, gg, Inu
      integer :: j, k
      real(dp) :: nuKN, urad, Iind, Ibol
      do k = 1, size(gg)
         nuKN = 3d0 * mass_e * cLight**2 / (4d0 * hPlanck * gg(k))
         Ibol = 0d0
         freqloop: do j = 2, size(freqs)
            if ( freqs(j) >= nuKN ) exit freqloop
            if ( Inu(j - 1) > 1d-200 .and. Inu(j) > 1d-200) then
               Iind = -dlog(Inu(j) / Inu(j - 1)) / dlog(freqs(j) / freqs(j - 1))
               if ( Iind > 8d0 ) Iind = 8d0
               if ( Iind < -8d0 ) Iind = -8d0
               Ibol = Ibol + Inu(j - 1) * freqs(j - 1) * Pinteg(freqs(j) / freqs(j - 1), Iind, 1d-9)
            end if
         end do freqloop
      end do
      urad = 4d0 * pi * Ibol / cLight
   end function IC_cool


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
   subroutine ssc_emissivity(fout, gg, nn, Imbs, emiss)
      implicit none
      real(dp), intent(in), dimension(:) :: gg, nn, Imbs, fout
      real(dp), intent(out), dimension(:) :: emiss
      integer :: j, k, Nf
      real(dp) :: g1, g2, jbol0, gmin, gmax
      real(dp), dimension(size(fout)) :: fin, I0
      Nf = size(fout)
      gmin = gg(1)
      gmax = gg(size(gg, dim=1))
      fin = fout
      !$OMP  PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j, g1, g2, I0, jbol0)
      fout_loop: do k=1, Nf
         fin_loop: do j=1, Nf
            g1 = dmax1(dsqrt(0.25d0 * fout(k) / fin(j)), gmin)
            g2 = dmin1(0.75d0 * mass_e * cLight**2 / (hPlanck * fin(j)), gmax)
            if ( g1 >= g2 ) then
               I0(j) = 0d0
            else
               I0(j) = Imbs(j) * ssc_qromb(fin(j), fout(k), g1, g2, gg, nn) / fin(j)**2
            end if
         end do fin_loop
         !  :::::  improvised trapezoidal rule  :::::
         jbol0 = fin(1) * I0(1) + fin(Nf) * I0(Nf)
         jbol0 = 0.5d0 * dlog(fin(Nf) / fin(1)) * &
            (jbol0 + 2d0 * sum(fin(2:Nf - 1) * I0(2:Nf - 1))) / dble(Nf - 1)
         emiss(k) = 0.25d0 * fout(k) * sigmaT * jbol0
      end do fout_loop
      !$OMP END PARALLEL DO

   end subroutine ssc_emissivity

   !
   !   ----------{   Romberg integrator   }----------
   !
   function ssc_qromb(fin, fout, a, b, gg, nn) result(qromb)
      implicit none
      real(dp), intent(in) :: fin, fout, a, b
      real(dp), intent(in), dimension(:) :: nn, gg
      integer, parameter :: jmax = 25, jmaxp = jmax + 1, kq = 6, km = kq - 1
      real(dp), parameter :: eps = 1d-4
      integer :: jq
      real(dp) :: dqromb, qromb
      real(dp), dimension(jmaxp) :: h, s
      h(1) = 1d0
      do jq = 1, jmax
         call ssc_trapzd(fin, fout, dlog(a), dlog(b), s(jq), jq, gg, nn)
         if ( jq >= kq ) then
            call polint(h(jq - km:jq), s(jq - km:jq), 0.0d0, qromb, dqromb)
            if ( dabs(dqromb) <= eps * dabs(qromb) ) return
         end if
         s(jq + 1) = s(jq)
         h(jq + 1) = 0.25d0 * h(jq)
      end do
      print*, fin, fout, a, b, qromb, dqromb
      call an_error('ssc_qromb: too many steps')
   end function ssc_qromb

   !
   !   ----------{   Trapezoid   }----------
   !
   subroutine ssc_trapzd(fin, fout, a, b, s, n, gg, nn)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: a, b, fout, fin
      real(dp), intent(in), dimension(:) :: gg, nn
      real(dp), intent(inout) :: s
      integer :: it, jt
      real(dp) :: del, fsum, x
      if ( n == 1 ) then
         s = 0.5d0 * (b - a) * ( &
         dexp(a) * ssc_integrand(fin, fout, dexp(a), gg, nn) + &
         dexp(b) * ssc_integrand(fin, fout, dexp(b), gg, nn) )
      else
         it = 2**(n - 2)
         del = (b - a) / dble(it)
         x = a + 0.5d0 * del
         fsum = 0d0
         do jt=1,it
            fsum = fsum + dexp(x) * ssc_integrand(fin, fout, dexp(x), gg, nn)
            x = x + del
         end do
         s = 0.5d0 * (s + del * fsum)
      end if
   end subroutine ssc_trapzd

   !
   !   ----------{   Dispersion function   }----------
   !
   function ssc_integrand(fin, fout, gev, gg, nn) result(integrand)
      implicit none
      real(dp), intent(in) :: fin, fout, gev
      real(dp), intent(in), dimension(:) :: gg, nn
      integer :: kk, Ng
      real(dp) :: fIC, nnev, beta, integrand, qq, foutfin, bplus, bminus
      Ng = size(gg, dim=1)
      kk = minloc(dabs(gev - gg), dim=1)

      if ( nn(kk) < 1d-100 ) then
         integrand = 0d0
         return
      end if

      beta = bofg(gev)
      foutfin = fout / fin
      bplus = 1d0 + beta
      bminus = 1d0 - beta

      if ( bminus / bplus <= foutfin .and. foutfin <= 1d0 ) then
         fIC = bplus * foutfin - bminus
      else if ( 1d0 <= foutfin .and. foutfin <= bplus / bminus ) then
         fIC = bplus - bminus * foutfin
      else
         integrand = 0d0
         return
      end if

      if ( any(dabs(gev - gg) == 0d0) ) then
         integrand = one_over_gb2(gev) * nn(kk) * fIC
         return
      else
         if ( gev < gg(kk) .or. gev > gg(Ng) ) then
            qq = dlog(nn(kk) / nn(kk - 1)) / dlog(gg(kk) / gg(kk - 1))
         else
            qq = dlog(nn(kk + 1) / nn(kk)) / dlog(gg(kk + 1) / gg(kk))
         end if
         if ( qq < -8d0 ) qq = -8d0
         if ( qq >  8d0 ) qq =  8d0

         nnev = nn(kk) * (gev / gg(kk))**qq
         integrand = one_over_gb2(gev) * nnev * fIC
         return
      end if
   end function ssc_integrand


   subroutine SSC_pwlEED(jSSC, nu, Imbs, n, g)
      implicit none
      real(dp), intent(in), dimension(:) :: nu, n, g, Imbs
      real(dp), intent(out), dimension(:) :: jSSC
      integer :: i, j, k, Ng, Nf
      real(dp) :: w1, w2, gmx_star, e0, l, q, f1, f2, emis, g1, g2, s1, s2
      Ng = size(g, dim=1)
      Nf = size(nu, dim=1)
      !$OMP  PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(i, j, k, emis, w1, w2, q, l, g1, g2, e0, f1, f2, s1, s2)
      nu_loop: do j = 1, Nf
         jSSC(j) = 0d0
         nu1_loop: do i = 1, Nf - 1
            intens_if: if ( Imbs(i) > 1d-200 .and. Imbs(i + 1) > 2d-100 ) then
               l = -dlog(Imbs(i + 1) / Imbs(i)) / dlog(nu(i + 1) / nu(i))
               if ( l > 8d0 ) l = 8d0
               if ( l < -8d0 ) l = -8d0
               e0 = mass_e * cLight**2 / (hPlanck * nu(i + 1))
               f1 = 0.25d0 * nu(j) / nu(i)
               f2 = 0.25d0 * nu(j) / nu(i + 1)
               g2 = dmin1(g(Ng), e0)
               g1 = dmax1(g(1), dsqrt(f1))
               g1g2_cond: if ( g1 < g2 ) then
                  g_loop: do k = 1, Ng - 1
                     if ( g(k) < g1 ) cycle g_loop
                     if ( g(k) > g2 ) exit g_loop
                     e_dist: if ( n(k) > 1d-200 .and. n(k + 1) > 1d-200 ) then
                        q = -dlog(n(k + 1) / n(k)) / dlog(g(k + 1) / g(k))
                        if ( q > 8d0 ) q = 8d0
                        if ( q < -8d0 ) q = -8d0
                        s1 = 0.5d0 * (q - 1d0)
                        s2 = 0.5d0 * (2d0 * l - q - 1d0)
                        if ( s1 > 8d0 ) s1 = 8d0
                        if ( s1 < -8d0 ) s1 = -8d0
                        if ( s2 > 8d0 ) s2 = 8d0
                        if ( s2 < -8d0 ) s2 = -8d0
                        gmx_star = dmin1(g(k + 1), e0)
                        w1 = dmin1(f1, gmx_star**2)
                        w2 = dmax1(f2, 0.25d0)
                        contrib_if: if ( w1 > w2 ) then
                           if ( 0.25d0 < f1 .and. f1 < g(k)**2 ) then
                              emis = sscG1ISO(1d0 / gmx_star**2, 1d0 / g(k)**2, w2, w1, s1, s2)
                           else if ( f2 <= g(k)**2 .and. g(k)**2 <= f1 ) then
                              emis = sscG1ISO(1d0 / gmx_star**2, 1d0 / g(k)**2, w2, g(k)**2, s1, s2) + &
                                 sscG2ISO(1d0 / gmx_star**2, 1d0, g(k)**2, w1, s1, s2)
                           else if ( g(k)**2 < f2 .and. f2 <= gmx_star**2 ) then
                              emis = sscG2ISO(1d0 / gmx_star**2, 1d0, w2, w1, s1, s2)
                           else
                              emis = 0d0
                           end if
                           jSSC(j) = jSSC(j) + emis * n(k) * g(k)**q * Imbs(i) * f1**(-l) * sigmaT
                        end if contrib_if
                     end if e_dist
                  end do g_loop
               end if g1g2_cond
            end if intens_if
         end do nu1_loop   
      end do nu_loop
      !$OMP END PARALLEL DO
   end subroutine SSC_pwlEED


   subroutine EIC_pwlEED(jEIC, nu, uext, nuext, n, g)
      implicit none
      real(dp), intent(in) :: uext, nuext
      real(dp), intent(in), dimension(:) :: nu, n, g
      real(dp), intent(out), dimension(:) :: jEIC
      real(dp), parameter :: eps = 1d-9
      integer :: j, k, Ng, Nf
      real(dp) :: w, gmx_star, e0, q, q1, q2, emis
      Ng = size(g, dim=1)
      Nf = size(nu, dim=1)
      e0 = mass_e * cLight**2 / (hPlanck * nuext)
      !$OMP  PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j, k, w, emis, q, q1, q2, gmx_star)
      nu_loop: do j = 1, Nf
         w = 0.25d0 * nu(j) / nuext
         jEIC(j) = 0d0
         g_loop: do k = 1, Ng - 1
            gmx_star = dmin1(g(k + 1), e0)
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
                  emis = (w / gmx_star**2)**q2 * ( Pinteg((gmx_star / g(k))**2, -q1, eps) - (w / gmx_star**2) * Pinteg((gmx_star / g(k))**2, -q2, eps) )
               else if ( g(k)**2 < w .and. w <= gmx_star**2 ) then
                  emis = (w / gmx_star**2)**q2 * ( Pinteg(gmx_star**2 / w, -q1, eps) - (w / gmx_star**2) * Pinteg(gmx_star**2 / w, -q2, eps) )
               else
                  cycle g_loop
               end if contrib_if
               jEIC(j) = jEIC(j) + emis * n(k) * g(k)**q / w**q1
            end if e_dist
         end do g_loop
         jEIC(j) = jEIC(j) * cLight * sigmaT * uext * 0.25d0 / (pi * nuext)
      end do nu_loop
      !$OMP END PARALLEL DO
   end subroutine EIC_pwlEED

end module radiation
