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
   ! #    # #####   ####     ###### #    # #  ####   ####
   ! ##  ## #    # #         #      ##  ## # #      #
   ! # ## # #####   ####     #####  # ## # #  ####   ####
   ! #    # #    #      #    #      #    # #      #      #
   ! #    # #    # #    #    #      #    # # #    # #    #
   ! #    # #####   ####     ###### #    # #  ####   ####
   subroutine mbs_emissivity(jnu, freq, gg, nn, B)
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
            jnu = jnu + j_mb(freq, B, nn(k), gg(k), gg(k + 1), qq, RMA_new)
         end if
      end do calc_jnu
      if ( jnu < 1d-200 ) jnu = 0d0
   end subroutine mbs_emissivity


   !  #    # #####   ####       ##   #####   ####   ####  #####
   !  ##  ## #    # #          #  #  #    # #      #    # #    #
   !  # ## # #####   ####     #    # #####   ####  #    # #    #
   !  #    # #    #      #    ###### #    #      # #    # #####
   !  #    # #    # #    #    #    # #    # #    # #    # #   #
   !  #    # #####   ####     #    # #####   ####   ####  #    #
   subroutine mbs_absorption(anu, freq, gg, nn, B)
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
            anu = anu + a_mb(freq, B, nn(k), gg(k), gg(k + 1), qq, RMA_new)
         end if
      end do calc_anu
      if ( anu < 1d-200 ) anu = 0d0
   end subroutine mbs_absorption


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
            u = 0.5d0 - 1d0 / tau**2
         else if ( tau >= 0.01d0 .and. tau <= 100d0 ) then
            u = 0.5d0 * (1d0 - 2d0 * (1d0 - (1d0 + tau) * dexp(-tau)) / tau**2)
         else
            u = (tau / 3d0) - 0.125d0 * tau**2
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
         g2 = dmin1(3d0 * mass_e * cLight**2 / (4d0 * hPlanck * fin(j)), gmax)
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
            gKN = mass_e * cLight**2 / (4d0 * hPlanck * nu(j + 1))
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
                           emis = sscG1ISO(1d0 / gmx_star**2, 1d0 / g(k)**2, w2, w1, s1, s2)
                        else if ( f2 <= g(k)**2 .and. g(k)**2 <= f1 ) then
                           emis = sscG1ISO(1d0 / gmx_star**2, 1d0 / g(k)**2, w2, g(k)**2, s1, s2) + &
                              sscG2ISO(1d0 / gmx_star**2, 1d0, g(k)**2, w1, s1, s2)
                        else if ( g(k)**2 < f2 .and. f2 <= gmx_star**2 ) then
                           emis = sscG2ISO(1d0 / gmx_star**2, 1d0, w2, w1, s1, s2)
                        else
                           cycle g_loop
                        end if
                        jnu = jnu + emis * n(k) * g(k)**q * Inu(j) * sigmaT * f1**(-l)
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
      gKN = mass_e * cLight**2 / (4d0 * hPlanck * nuext)
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
               emis = (w / gmx_star**2)**q2 * ( Pinteg((gmx_star / g(k))**2, -q1, eps) - (w / gmx_star**2) * Pinteg((gmx_star / g(k))**2, -q2, eps) )
            else if ( g(k)**2 < w .and. w <= gmx_star**2 ) then
               emis = (w / gmx_star**2)**q2 * ( Pinteg(gmx_star**2 / w, -q1, eps) - (w / gmx_star**2) * Pinteg(gmx_star**2 / w, -q2, eps) )
            else
               cycle g_loop
            end if contrib_if
            jnu = jnu + emis * n(k) * g(k)**q * w**(-q1)
         end if e_dist
      end do g_loop
      jnu = jnu * cLight * sigmaT * uext / (4d0 * pi * nuext)
   end subroutine IC_iso_monochrom


   !   ####   ####   ####  #      # #    #  ####
   !  #    # #    # #    # #      # ##   # #    #
   !  #      #    # #    # #      # # #  # #
   !  #      #    # #    # #      # #  # # #  ###
   !  #    # #    # #    # #      # #   ## #    #
   !   ####   ####   ####  ###### # #    #  ####
   subroutine rad_cool(dotg, gg, freqs, uu, withKN)
      implicit none
      real(dp), intent(in), dimension(:) :: freqs, gg, uu
      real(dp), intent(out), dimension(:) :: dotg
      logical, intent(in) :: withKN
      integer :: j, k, Ng, Nf
      real(dp) :: uind, urad_const, usum, xi_c, xi_rat
      real(dp), dimension(size(gg, dim=1), size(freqs, dim=1)) :: xi, uxi
      urad_const = 4d0 * sigmaT * cLight / (3d0 * energy_e)
      xi_c = 4d0 * hPlanck / energy_e
      Ng = size(gg, dim=1)
      Nf = size(freqs, dim=1)
      do k = 1, Ng
         do j = 1, Nf
            xi(k, j) = xi_c * gg(k) * freqs(j)
            uxi(k, j) = uu(j) * xi_c * gg(k)
         end do
      end do
      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(k, j, uind, xi_rat, usum)
      do k = 1, Ng
         usum = 0d0
         freqloop: do j = 1, Nf - 1
            if ( uxi(k, j + 1) > 1d-100 .and. uxi(k, j) > 1d-100 ) then
               xi_rat = xi(k, j + 1) / xi(k, j)
               uind = dlog(uxi(k, j + 1) / uxi(k, j)) / dlog(xi_rat)
               if ( uind > 8d0 ) uind = 8d0
               if ( uind < -8d0 ) uind = -8d0
               if ( xi_c * gg(k) * freqs(j) > 1d0 .and. withKN ) then
                  usum = usum + 4.5d0 * uxi(k, j) &
                        * (Qinteg(xi_rat, uind + 2d0, 1d-6) &
                        + (dlog(xi(k, j)) - (11d0 / 6d0)) &
                        * Pinteg(xi_rat, uind + 2d0, 1d-6)) &
                        / xi(k, j)
               else if ( xi_c * gg(k) * freqs(j) > 1d0 .and. .not. withKN ) then
                  cycle freqloop
               else
                  usum = usum + uxi(k, j) * xi(k, j) * Pinteg(xi_rat, uind, 1d-6)
               end if
            end if
         end do freqloop
         dotg(k) = urad_const * usum / xi_c**2
      end do
      !$OMP END PARALLEL DO
   end subroutine rad_cool


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
      zero = zeros1D(Nf)
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
      zero = zeros1D(Nf)
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

end module radiation
