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
   function mbs_emissivity(freqs, gg, nn, B) result(jnu)
      implicit none
      real(dp), intent(in) :: B
      real(dp), intent(in), dimension(:) :: freqs, gg, nn
      integer :: j, k, Ng, Nf
      real(dp) :: qq
      real(dp), dimension(size(freqs)) :: jnu
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
      end do freqs_loop
      !$OMP END PARALLEL DO
   end function mbs_emissivity


   ! #    # #####   ####       ##   #####   ####   ####  #####
   ! ##  ## #    # #          #  #  #    # #      #    # #    #
   ! # ## # #####   ####     #    # #####   ####  #    # #    #
   ! #    # #    #      #    ###### #    #      # #    # #####
   ! #    # #    # #    #    #    # #    # #    # #    # #   #
   ! #    # #####   ####     #    # #####   ####   ####  #    #
   function mbs_absorption(freqs, gg, nn, B) result(anu)
      implicit none
      real(dp), intent(in) :: B
      real(dp), intent(in), dimension(:) :: freqs, gg, nn
      integer :: j, k, Ng, Nf
      real(dp) :: qq
      real(dp), dimension(size(freqs)) :: anu
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
      end do freqs_loop
      !$OMP END PARALLEL DO
   end function mbs_absorption


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


   !  ####   ####   ####   ####      ####   ####  ###### ######
   ! #      #      #    # #    #    #    # #    # #      #
   !  ####   ####  #      #         #      #    # #####  #####
   !      #      # #      #         #      #    # #      #
   ! #    # #    # #    # #    #    #    # #    # #      #
   !  ####   ####   ####   ####      ####   ####  ###### #
   function ssc_cool_coef(nu0B,gg,freqs,Inu) result(nu0)
      implicit none
      real(dp), intent(in) :: nu0B
      real(dp), intent(in), dimension(:) :: freqs,gg,Inu
      integer :: j,k
      real(dp) :: nuKN,urad,Iind,Ibol
      real(dp), dimension(size(gg)) :: nu0
      do k = 1, size(gg)
         nuKN = 3d0 * mass_e * cLight**2 / (4d0 * hPlanck * gg(k))
         Ibol = 0d0
         freqloop: do j = 2, size(freqs)
            if ( freqs(j) >= nuKN ) exit freqloop
            if ( Inu(j - 1) > 1d-50 .and. Inu(j) > 1d-50) then
               Iind = -dlog(Inu(j) / Inu(j - 1)) / dlog(freqs(j) / freqs(j - 1))
               if ( Iind > 8d0 ) cycle freqloop !Iind = 8d0
               if ( Iind < -8d0 ) cycle freqloop !Iind = -8d0
               Ibol = Ibol + Inu(j - 1) * freqs(j - 1) * Pinteg(freqs(j) / freqs(j - 1),Iind,1d-9)
            end if
         end do freqloop
         urad = 4d0 * pi * Ibol / cLight
         nu0(k) = nu0B + 4d0 * sigmaT * urad / (3d0 * mass_e * cLight)
      end do
   end function ssc_cool_coef


   !  ####   ####   ####     ###### #    # #  ####   ####
   ! #      #      #    #    #      ##  ## # #      #
   !  ####   ####  #         #####  # ## # #  ####   ####
   !      #      # #         #      #    # #      #      #
   ! #    # #    # #    #    #      #    # # #    # #    #
   !  ####   ####   ####     ###### #    # #  ####   ####
   subroutine ssc_emissivity(fout, gmin, gmax, gg, nn, Imbs, emiss)
      implicit none
      real(dp), intent(in) :: gmin, gmax
      real(dp), intent(in), dimension(:) :: gg, nn, Imbs, fout
      real(dp), intent(out), dimension(:) :: emiss
      integer :: j, k, Nf
      real(dp) :: g1, g2, jbol0
      real(dp), dimension(size(fout)) :: fin, I0

      Nf = size(fout)
      fin = fout

      !$OMP  PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(j, g1, g2, I0, jbol0)
      fout_loop: do k=1,Nf
         fin_loop: do j=1,Nf
            g1 = dmax1(dsqrt(0.25d0 * fout(k) / fin(j)), gmin)
            g2 = dmin1(0.75d0 * mass_e * cLight2 / (hPlanck * fin(j)), gmax)
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


   !  #######  #####    #    #####
   !       #  #     #  ##   #     #
   !      #   #       # #         #
   !     #     #####    #    #####
   !    #           #   #   #
   !   #      #     #   #   #
   !  #######  #####  ##### #######
   function ssccZS12(gg,nn,B,R) result(nu0)
      implicit none
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


   ! ###### #      #    # #    #
   ! #      #      #    #  #  #
   ! #####  #      #    #   ##
   ! #      #      #    #   ##
   ! #      #      #    #  #  #
   ! #      ######  ####  #    #
   function flux_dens(Inu,dL,z,dopp,R) result(Fnu)
      ! ************************************************************************
      !  Description:
      !     Using Eq. (13) in Zheng & Zhang, 2011, ApJ, 728, 105.
      !  Arguments:
      !     - Inu : intensity
      !     - dL  : luminosity distance
      !     - z   : redshift
      !     - dopp: Doppler factor
      !     - R   : radius of the emision region
      !     - Fnu : flux density measured by a distant observer
      ! ************************************************************************
      implicit none
      real(dp), intent(in) :: dL,z,dopp,R
      real(dp), intent(in), dimension(:) :: Inu
      real(dp), dimension(size(Inu)) :: Fnu
      Fnu = pi * R**2 * dopp**3 * (1d0 + z) * Inu / dL**2
   end function flux_dens

end module radiation
