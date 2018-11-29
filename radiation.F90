module radiation
   use data_types
   use constants
   use misc
   use pwl_integ
   use SRtoolkit
   use anaFormulae
   !$ use omp_lib
   use magnetobrem
   implicit none

contains
   !
   ! #    # #####   ####     ###### #    # #  ####   ####
   ! ##  ## #    # #         #      ##  ## # #      #
   ! # ## # #####   ####     #####  # ## # #  ####   ####
   ! #    # #    #      #    #      #    # #      #      #
   ! #    # #    # #    #    #      #    # # #    # #    #
   ! #    # #####   ####     ###### #    # #  ####   ####
   !
   function mbs_emissivity(freqs, gg, nn, B, mbs_or_syn) result(jnu)
      implicit none
      real(dp), intent(in) :: B
      real(dp), intent(in), dimension(:) :: freqs, gg, nn
      logical, intent(in) :: mbs_or_syn
      integer :: j, k, kglob, Ng, Nf
      real(dp) :: qq, n_globgmx
      real(dp), dimension(size(freqs)) :: jnu
      jnu = 0d0
      Ng = size(gg, dim=1)
      Nf = size(freqs, dim=1)
      !$OMP PARALLEL DO ORDERED COLLAPSE(2) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(qq)
      freqs_loop: do j = 1, Nf
         calc_jnu: do k = 1, Ng - 1
            if ( nn(k) > 1d-100 .and. nn(k + 1) > 1d-100) then
               qq = -dlog(nn(k + 1) / nn(k)) / dlog(gg(k + 1) / gg(k))
               if ( mbs_or_syn ) then
                  if ( qq < SS(1) ) qq = SS(1)
                  if ( qq > SS(numSS) ) qq = SS(numSS)
                  if ( freqs(j) > dexp(XX(numXX)) * nuconst * B ) then
                     jnu(j) = jnu(j) + &
                        j_mb_qromb(freqs(j), B, nn(k), gg(k), gg(k + 1), qq, chunche_c100g20, RMA_new)
                  else
                     if ( gg(k) < globgmax .and. gg(k + 1) <= globgmax ) then
                        jnu(j) = jnu(j) + &
                           j_mb(freqs(j), B, nn(k), gg(k), gg(k + 1), qq, RMA_new, chunche_c100g20, jtable)
                     else if (gg(k) < globgmax .and. gg(k + 1) > globgmax) then
                        kglob = locate(gg, globgmax, in_bounds = .true.)
                        n_globgmx = nn(kglob) * (gg(kglob + 1) / gg(kglob))**qq
                        jnu(j) = jnu(j) + &
                           j_mb(freqs(j), B, nn(k), gg(k), globgmax, qq, RMA_new, chunche_c100g20, jtable) + &
                           j_mb_qromb(freqs(j), B, n_globgmx, globgmax, gg(k + 1), qq, chunche_c100g20, RMA_new)
                     else
                        jnu(j) = jnu(j) + &
                           j_mb_qromb(freqs(j), B, nn(k), gg(k), gg(k + 1), qq, chunche_c100g20, RMA_new)
                     end if
                  end if
               else
                  if ( qq > 8d0 ) qq = 8d0
                  if ( qq < -8d0 ) qq = -8d0
                  jnu(j) = jnu(j) + j_mb_qromb(freqs(j), B, nn(k), gg(k), gg(k + 1), qq, 1d0, RMA_new)
               end if
            end if
         end do calc_jnu
      end do freqs_loop
      !$OMP END PARALLEL DO
   end function mbs_emissivity


   !
   ! #    # #####   ####       ##   #####   ####   ####  #####
   ! ##  ## #    # #          #  #  #    # #      #    # #    #
   ! # ## # #####   ####     #    # #####   ####  #    # #    #
   ! #    # #    #      #    ###### #    #      # #    # #####
   ! #    # #    # #    #    #    # #    # #    # #    # #   #
   ! #    # #####   ####     #    # #####   ####   ####  #    #
   !
   function mbs_absorption(freqs, gg, nn, B, mbs_or_syn) result(anu)
      implicit none
      real(dp), intent(in) :: B
      real(dp), intent(in), dimension(:) :: freqs,gg,nn
      logical, intent(in) :: mbs_or_syn
      integer :: j,k,kglob
      real(dp) :: qq,n_globgmx
      real(dp), dimension(size(freqs)) :: anu,aa

      !$OMP PARALLEL DO ORDERED COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(qq,k)
      freqs_loop: do j=1, size(freqs)
         aa(j) = 0d0
         !$OMP ORDERED
         calc_jnu: do k=1, size(gg) - 1
            if ( nn(k) > 1d-100 ) then
               qq = -dlog(nn(k + 1) / nn(k)) / dlog(gg(k + 1) / gg(k))
               if ( mbs_or_syn ) then
                  if ( qq < SS(1) ) qq = SS(1)
                  if ( qq > SS(numSS) ) qq = SS(numSS)
                  if ( freqs(j) > dexp(XX(numXX)) * nuconst * B ) then
                     aa(j) = aa(j) + a_mb_qromb(freqs(j), B, nn(k), gg(k), gg(k + 1), qq, chunche_c100g20, RMA_new)
                  else
                     if ( gg(k) < globgmax .and. gg(k + 1) <= globgmax ) then
                        aa(j) = aa(j) + &
                           a_mb(freqs(j), B, nn(k), gg(k), gg(k + 1), qq, RMA_new, chunche_c100g20, jtable)
                     else if ( gg(k) < globgmax .and. gg(k + 1) > globgmax ) then
                        kglob = locate(gg, globgmax, in_bounds=.true.)
                        n_globgmx = nn(kglob) * ( gg(kglob + 1) / gg(kglob) )**qq
                        aa(j) = aa(j) + &
                           a_mb(freqs(j), B, nn(k), gg(k), globgmax, qq, RMA_new, chunche_c100g20, jtable) + &
                           a_mb_qromb(freqs(j), B, n_globgmx, globgmax, gg(k + 1), qq, chunche_c100g20, RMA_new)
                     else
                        aa(j) = aa(j) + &
                           a_mb_qromb(freqs(j), B, nn(k), gg(k), gg(k + 1), qq, chunche_c100g20, RMA_new)
                     end if
                  end if
               else
                  if ( qq > 8d0 ) qq = 8d0
                  if ( qq < 8d0 ) qq = -8d0
                  aa(j) = aa(j) + a_mb_qromb(freqs(j), B, nn(k), gg(k), gg(k + 1), qq, 1d0, RMA_new)
               end if
            end if
         end do calc_jnu
         !$OMP END ORDERED
         if ( aa(j) < 1d-100 ) aa(j) = 0d0
         anu(j) = aa(j)
      end do freqs_loop
      !$OMP END PARALLEL DO
   end function mbs_absorption


   !
   ! # #    # ##### ###### #    #  ####  # ##### #   #
   ! # ##   #   #   #      ##   # #      #   #    # #
   ! # # #  #   #   #####  # #  #  ####  #   #     #
   ! # #  # #   #   #      #  # #      # #   #     #
   ! # #   ##   #   #      #   ## #    # #   #     #
   ! # #    #   #   ###### #    #  ####  #   #     #
   !
   subroutine solRadTrans_pathIndep(Inu, s, jnu, anu, Is0)
      ! Description:
      !
      !   This function solves the radiative transfer equation assuming that
      !   the emissivity and absorption are independent of the path. We take
      !   Eq. (1.30) of Rybicki & Lightman (1979) and make
      !
      !        tau_nu = s * alpha_nu
      !
      !   for the optically thick case.
      implicit none
      real(dp), intent(in) :: s
      real(dp), intent(in), dimension(:) :: jnu, anu, Is0
      real(dp), intent(out), dimension(:) :: Inu
      optional :: jnu, anu, Is0
      integer :: j
      real(dp) :: Snu
      logical :: with_emiss, with_absor, with_inten
      with_emiss = present(jnu)
      with_absor = present(anu)
      with_inten = present(Is0)
      do j=1,size(jnu)
         if ( with_emiss .and. jnu(j) > 1d-100 ) then
            if ( with_absor .and. anu(j) > 1d-100) then
               Snu = jnu(j) / anu(j)
               if ( with_inten ) then
                  Inu(j) = Snu + dexp(-anu(j) * s) * (Is0(j) - Snu)
               else
                  Inu(j) = Snu * (1d0 - dexp(-anu(j) * s))
               end if
            else
               if ( with_inten ) then
                  Inu(j) = Is0(j) + s * jnu(j)
               else
                  Inu(j) = s * jnu(j)
               end if
            end if
         else
            if ( with_absor .and. anu(j) > 1d-100) then
               if ( with_inten ) then
                  Inu(j) = Is0(j) * dexp(-anu(j) * s)
               else
                  Inu(j) = 0d0
               end if
            else
               if ( with_inten ) then
                  Inu(j) = Is0(j)
               else
                  Inu(j) = 0d0
               end if
            end if
         end if
      end do
   end subroutine solRadTrans_pathIndep

#if 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine solRadTrans_lumPath(Inu, s, jnu, anu)
      ! Description:
      !
      !   This function solves the radiative transfer equation assuming that
      !   the emissivity and absorption are independent of the path. We take
      !   Eq. (1.30) of Rybicki & Lightman (1979) and make
      !
      !        tau_nu = s * alpha_nu
      !
      !   for the optically thick case.
      implicit none
      real(dp), intent(in), dimension(:) :: s
      real(dp), intent(in), dimension(:, :) :: jnu, anu, Is0
      real(dp), intent(out), dimension(:, :) :: Inu
      optional :: jnu, anu, Is0
      integer :: j, i
      real(dp) :: Snu
      logical :: with_emiss, with_absor, with_inten
      with_emiss = present(jnu)
      with_absor = present(anu)
      with_inten = present(Is0)
      do j=1,size(jnu)
         do i=1,size(s)
         if ( with_emiss .and. jnu(j) > 1d-100 ) then
            if ( with_absor .and. anu(j) > 1d-100) then
               Snu = jnu(j) / anu(j)
               if ( with_inten ) then
                  Inu(j) = Snu + dexp(-anu(j) * s) * (Is0(j) - Snu)
               else
                  Inu(j) = Snu * (1d0 - dexp(-anu(j) * s))
               end if
            else
               if ( with_inten ) then
                  Inu(j) = Is0(j) + s * jnu(j)
               else
                  Inu(j) = s * jnu(j)
               end if
            end if
         else
            if ( with_absor .and. anu(j) > 1d-100) then
               if ( with_inten ) then
                  Inu(j) = Is0(j) * dexp(-anu(j) * s)
               else
                  Inu(j) = 0d0
               end if
            else
               if ( with_inten ) then
                  Inu(j) = Is0(j)
               else
                  Inu(j) = 0d0
               end if
            end if
         end if
      end do
   end do

   contains

      !
      !   ----------{   Romberg integrator   }----------
      !
      function blob_qromb(tini, tfin, emiss, tem) result(qromb)
         implicit none
         real(dp), intent(in) :: tini, tfin
         real(dp), intent(in), dimension(:) :: tem, emiss
         integer, parameter :: jmax = 25, jmaxp = jmax + 1, kq = 6, km = kq - 1
         integer :: jq
         real(dp), parameter :: eps = 1d-3
         real(dp), dimension(jmaxp) :: h, s
         real(dp) :: dqromb, qromb
         if ( tini == tfin ) then
            qromb = 0d0
            return
         end if
         h(1) = 1d0
         do jq=1,jmax
            call blob_trapzd(dlog(tini), dlog(tfin), s(jq), jq, emiss, tem)
            ! call blob_trapzd(dlog(tini), dlog(tfin), s(jq), jq, dlog(emiss), dlog(tem))
            if ( jq >= kq ) then
               call polint(h(jq - km:jq), s(jq - km:jq), 0.0d0, qromb, dqromb)
               if ( dabs(dqromb) <= eps * dabs(qromb) ) return
            end if
            s(jq + 1) = s(jq)
            h(jq + 1) = 0.25d0 * h(jq)
         end do
         print*,qromb,dqromb
         call an_error('blob_qromb: too many steps')
      end function blob_qromb

      !
      !   ----------{   Trapezoid   }----------
      !
      subroutine blob_trapzd(tini, tfin, s, n, emiss, tem)
         implicit none
         real(dp), intent(in) :: tini, tfin
         real(dp), intent(in), dimension(:) :: tem, emiss
         real(dp), intent(inout) :: s
         integer, intent(in) :: n
         real(dp) :: del, fsum, x, jini, jfin, jx, dj
         integer :: it, jt
         if ( n == 1 ) then
            ! call polint(tem, emiss, dexp(tini), jini, dj)
            ! call polint(tem, emiss, dexp(tfin), jfin, dj)
            ! call polint(tem, emiss, tini, jini, dj)
            ! call polint(tem, emiss, tfin, jfin, dj)
            s = 0.5d0 * (tfin - tini) * ( &
               dexp(tini) * blob_integrand(dexp(tini), emiss, tem) + &
               dexp(tfin) * blob_integrand(dexp(tfin), emiss, tem) )
               ! dexp(tini) * jini + dexp(tfin) * jfin )
               ! dexp(tini) * dexp(jini) + dexp(tfin) * dexp(jfin) )
         else
            it = 2**(n - 2)
            del = (tfin - tini) / dble(it)
            x = tini + 0.5d0 * del
            fsum = 0d0
            do jt = 1, it
               ! call polint(tem, emiss, dexp(x), jx, dj)
               ! call polint(tem, emiss, x, jx, dj)
               ! fsum = fsum + dexp(x) * jx
               ! fsum = fsum + dexp(x) * dexp(jx)
               fsum = fsum + dexp(x) * blob_integrand(dexp(x), emiss, tem)
               x = x + del
            end do
            s = 0.5d0 * (s + del * fsum)
         end if
      end subroutine blob_trapzd

      function blob_integrand(t, em, te) result(res)
         implicit none
         double precision, intent(in) :: t
         double precision, intent(in), dimension(:) :: te, em
         integer :: t_pos, Nt
         double precision :: res, sind
         Nt = size(te)
         t_pos = locate(te, t, .true.)
         if ( t_pos >= Nt ) then
            t_pos = Nt - 1
            ! res = 0d0
         end if
         if ( em(t_pos) > 1d-100 .and. em(t_pos + 1) > 1d-100 ) then
            sind = -dlog(em(t_pos + 1) / em(t_pos)) / dlog(te(t_pos + 1) / te(t_pos))
            res = em(t_pos) * (t / te(t_pos))**sind
         else
            res = 0d0
         end if
      end function blob_integrand
   end subroutine solRadTrans_lumPath
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif


   !
   !  ####   ####   ####   ####      ####   ####  ###### ######
   ! #      #      #    # #    #    #    # #    # #      #
   !  ####   ####  #      #         #      #    # #####  #####
   !      #      # #      #         #      #    # #      #
   ! #    # #    # #    # #    #    #    # #    # #      #
   !  ####   ####   ####   ####      ####   ####  ###### #
   !
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


   !
   !  ####   ####   ####     ###### #    # #  ####   ####
   ! #      #      #    #    #      ##  ## # #      #
   !  ####   ####  #         #####  # ## # #  ####   ####
   !      #      # #         #      #    # #      #      #
   ! #    # #    # #    #    #      #    # # #    # #    #
   !  ####   ####   ####     ###### #    # #  ####   ####
   !
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


   !
   !  #######  #####    #    #####
   !       #  #     #  ##   #     #
   !      #   #       # #         #
   !     #     #####    #    #####
   !    #           #   #   #
   !   #      #     #   #   #
   !  #######  #####  ##### #######
   !
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


   !
   ! ###### #      #    # #    #
   ! #      #      #    #  #  #
   ! #####  #      #    #   ##
   ! #      #      #    #   ##
   ! #      #      #    #  #  #
   ! #      ######  ####  #    #
   !
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
