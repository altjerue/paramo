module dist_evol
   use data_types
   use constants
   use misc
   use pwl_integ
   use SRtoolkit
   implicit none

   integer :: d_Ng, d_Nt
   real(dp) :: d_dtinj, d_g1, d_g2, d_tcurr, d_tesc
   real(dp), allocatable, dimension(:) :: d_gg, d_qind, d_Q0, d_t
   real(dp), allocatable, dimension(:, :) :: d_nu0, d_dist
   public :: d_ng, d_t, d_dtinj, d_g1, d_g2, d_tcurr, d_qind, &
      d_Q0, d_Nt, d_dist, d_tesc

contains

   !
   ! # #    #      # ######  ####  ##### #  ####  #    #
   ! # ##   #      # #      #    #   #   # #    # ##   #
   ! # # #  #      # #####  #        #   # #    # # #  #
   ! # #  # #      # #      #        #   # #    # #  # #
   ! # #   ## #    # #      #    #   #   # #    # #   ##
   ! # #    #  ####  ######  ####    #   #  ####  #    #
   !
   function injection(t, dtinj, g, g1, g2, qind, th, Qth, Qnth) result(Qinj)
      implicit none
      real(dp), intent(in) :: g1, g2, th, Qth, Qnth, dtinj, t, qind
      real(dp), intent(in), dimension(:) :: g
      integer :: k
      real(dp), dimension(size(g)) :: Qinj
      
      if ( t <= dtinj ) then
         do k = 1, size(g)
            if ( Qth < 1d-100 ) then
               Qinj(k) = Qnth * powlaw_dis(g(k), g1, g2, qind)
            else if ( Qnth < 1d-100 ) then
               Qinj(k) = Qth * RMaxwell(g(k), th)
            else
               if ( g(k) >= g1 .and. g(k) <= g2 ) then
                  Qinj(k) = Qnth * powlaw_dis(g(k), g1, g2, qind) + &
                     Qth * RMaxwell(g(k), th)
               else
                  Qinj(k) = Qth * RMaxwell(g(k), th)
               end if
            end if
         end do
      else
         Qinj = 0d0
      end if
   end function injection


   !  #######          ######
   !  #       # #    # #     # # ######
   !  #       # ##   # #     # # #
   !  #####   # # #  # #     # # #####
   !  #       # #  # # #     # # #
   !  #       # #   ## #     # # #
   !  #       # #    # ######  # #
   subroutine FP_FinDif_difu(dt, g, nin, nout, nu0, DD, QQ, tesc)
      implicit none
      real(dp), intent(in) :: dt, tesc
      real(dp), intent(in), dimension(:) :: g, nin, nu0, DD, QQ
      real(dp), intent(out), dimension(:) :: nout
      real(dp), parameter :: eps = 1e-2
      integer :: i, Ng, Ng1
      real(dp), dimension(size(g)) :: gdot, dx, dxp2, dxm2, CCp2, CCm2, BBp2, BBm2, YYp2, YYm2, WWp2, WWm2, ZZp2, ZZm2, a, b, c, r

      Ng = size(g)
      Ng1 = Ng - 1

      gdot = nu0 * g**2
      dxp2(:Ng1) = g(2:) - g(:Ng1)
      dxp2(Ng) = dxp2(Ng1)
      dxm2(2:) = dxp2(:Ng1)
      dxm2(1) = dmin1(g(1) - 1d0, dxm2(2))
      dx = dsqrt(dxp2 * dxm2)

      CCp2(:Ng1) = 0.25d0 * (DD(2:) + DD(:Ng1))
      CCp2(Ng) = 0.25d0 * DD(Ng)
      CCm2(2:) = CCp2(:Ng1)
      CCm2(1) = 0.25d0 * DD(1)

      BBp2(:Ng1) = 0.5d0 * ( (DD(2:) - DD(:Ng1)) / dxp2(:Ng1) - (gdot(2:) + gdot(:Ng1)) )
      BBm2(2:) = 0.5d0 * ( (DD(2:) - DD(:Ng1)) / dxm2(2:) - (gdot(2:) + gdot(:Ng1)) )
      BBp2(Ng) = 0.5d0 * ( -DD(Ng) / dxp2(Ng) - gdot(Ng) )
      BBm2(1) = 0.5d0 * ( DD(1) / dxm2(1) - gdot(1) )

      WWp2 = dxp2 * BBp2 / CCp2
      WWm2 = dxm2 * BBm2 / CCm2

      do i = 1, Ng

         if ( 0.5d0 * WWp2(i) > 200d0 ) then
            ZZp2(i) = dexp(200d0)
         else if ( 0.5d0 * WWp2(i) < -200d0 ) then
            ZZp2(i) = dexp(-200d0)
         else
            ZZp2(i) = dexp(0.5d0 * WWp2(i))
         end if

         if ( 0.5d0 * WWm2(i) > 200d0 ) then
            ZZm2(i) = dexp(200d0)
         else if ( 0.5d0 * WWm2(i) < -200d0 ) then
            ZZm2(i) = dexp(-200d0)
         else
            ZZm2(i) = dexp(0.5d0 * WWm2(i))
         end if

         if ( 31d0 * WWp2(i)**6 < eps * 967680d0 ) then
            YYp2(i) = 1d0 - WWp2(i)**2 / 24d0 + 7d0 * WWp2(i)**4 / 5760d0
         else
            YYp2(i) = WWp2(i) / ( (1d0 - 1d0 / ZZp2(i)**2 ) * ZZp2(i) )
         end if

         if ( 31d0 * WWm2(i)**6 < eps * 967680d0 ) then
            YYm2(i) = 1d0 - WWm2(i)**2 / 24d0 + 7d0 * WWm2(i)**4 / 5760d0
         else
            YYm2(i) = WWm2(i) / ( (1d0 - 1d0 / ZZm2(i)**2 ) * ZZm2(i) )
         end if

      end do

      r = nin + dt * QQ
      c = dt * CCm2 * YYm2 / (dx * dxm2 * ZZm2)
      a = dt * CCp2 * YYp2 * ZZp2 / (dx * dxp2)
      b = 1d0 + dt * ( CCm2 * YYm2 * ZZm2 / dxm2 + CCp2 * YYp2 / (dxp2 * ZZp2) ) / dx + dt / tesc

      call tridag_ser(-a(:Ng1), b, -c(2:), r, nout)

   end subroutine FP_FinDif_difu


   subroutine FP_FinDif_cool(dt, g, nin, nout, nu0, QQ, tesc)
      implicit none
      real(dp), intent(in) :: dt, tesc
      real(dp), intent(in), dimension(:) :: g, nin, nu0, QQ
      real(dp), intent(out), dimension(:) :: nout
      integer :: Ng, Ng1
      real(dp), dimension(size(g)) :: gdot, dx, dxp2, dxm2, BBp2, BBm2, a, b, c, r

      Ng = size(g)
      Ng1 = Ng - 1

      gdot = nu0 * g**2
      dxp2(:Ng1) = g(2:) - g(:Ng1)
      dxp2(Ng) = dxp2(Ng1)
      dxm2(2:) = dxp2(:Ng1)
      dxm2(1) = dmin1(g(1) - 1d0, dxm2(2))
      dx = dsqrt(dxp2 * dxm2)

      BBp2(:Ng1) = 0.5d0 * (gdot(2:) + gdot(:Ng1))
      BBp2(Ng) = 0.5d0 * gdot(Ng)
      BBm2(2:) = BBp2(:Ng1)
      BBm2(1) = 0.5d0 * gdot(1)

      r = nin + dt * QQ
      a = zeros1D(Ng)
      c = -dt * BBp2 / dx
      b = 1d0 + dt / tesc + BBm2 * dt / dx
      ! b(1) = 1d0 + dt / tesc

      call tridag_ser(a(:Ng1), b, c(2:), r, nout)

   end subroutine FP_FinDif_cool

   ! ===========================================================================
   !
   !   ####   ####   ####  #            #      # #    # ######  ####
   !  #    # #    # #    # #            #      # ##   # #      #
   !  #      #    # #    # #      ##### #      # # #  # #####   ####
   !  #      #    # #    # #            #      # #  # # #           #
   !  #    # #    # #    # #            #      # #   ## #      #    #
   !   ####   ####   ####  ######       ###### # #    # ######  ####
   !
   subroutine distrib_setup(nn, gg, nu0, t, dtinj, tesc, g1, g2, qind, theta_e, Qth, Qnth)
      implicit none
      real(dp), intent(in) :: dtinj, tesc, g1, g2, theta_e, Qth, Qnth, qind
      real(dp), intent(in), dimension(:) :: t, gg
      real(dp), intent(in), dimension(:, :) :: nu0, nn
      integer :: k, Ng, Nt
      real(dp), dimension(size(gg, dim=1)) :: Qinj

      Ng = size(gg, dim=1)
      Nt = size(t, dim=1)
      call realloc(d_t, Nt)
      call realloc(d_gg, Ng)
      call realloc(d_Q0, Ng)
      call realloc(d_qind, Ng)
      call realloc(d_nu0, Nt, Ng)
      call realloc(d_dist, Nt, Ng)

      d_Ng = Ng
      d_Nt = Nt
      d_t = t
      d_tcurr = t(Nt)
      d_dtinj = dtinj
      d_tesc = tesc
      d_g1 = g1
      d_g2 = g2
      d_qind = qind
      d_nu0 = nu0
      d_dist = nn
      d_gg = gg

      Qinj = injection(d_tcurr, d_dtinj, gg, g1, g2, qind, theta_e, Qth, Qnth)

      do k = 1, Ng - 1
         if ( Qinj(k + 1) > 1d-100 .and. Qinj(k) > 1d-100 ) then
            d_qind(k) = -dlog(Qinj(k + 1) / Qinj(k)) / dlog(gg(k + 1) / gg(k))
            if ( d_qind(k) > 8d0 ) d_qind(k) = 8d0
            if ( d_qind(k) < -8d0 ) d_qind(k) = -8d0
            d_Q0(k) = Qinj(k) * gg(k)**d_qind(k)
         else
            d_qind(k) = 0d0
            d_Q0(k) = 0d0
         endif
      enddo
      d_qind(Ng) = d_qind(Ng - 1)
      d_Q0(Ng) = Qinj(Ng) * gg(Ng)**d_qind(Ng)

   end subroutine distrib_setup


   function distrib(gamma) result(nn)
      implicit none
      real(dp), intent(in) :: gamma
      integer :: k, i
      real(dp) :: nn, g0, ghigh, glow, t, g, tt, t0
      logical :: injecting

      injecting = .true.
      t = d_tcurr
      g = gamma
      k = locate(d_gg, g, .true.)
      i = d_Nt - 1
      do while ( i > 0 .and. k < d_Ng )
         g0 = g / (1d0 - d_nu0(i, k) * (t - d_t(i)) * g)
         if ( g0 >= d_gg(k + 1) ) then
            t = t - (d_gg(k + 1) - g) / (d_nu0(i, k) * g * d_gg(k + 1))
            g = d_gg(k + 1)
            k = k + 1
         else
            t = d_t(i)
            g = g0
            i = i - 1
         end if
      end do

      ! if ( k == d_Ng .and. t > 0d0 ) return
      ! k = locate(d_gg, g0, .true.)
      ! i = locate(d_t, t, .true.)

      nn = 0d0
      t = d_t(i)
      tt = dmin1( d_t(i + 1), t + d_tesc )
      contrib_loop: do while ( t < d_tcurr .and. g0 > d_gg(1) )

         g = g0 / (1d0 + d_nu0(i, k) * (tt - t) * g0)

         if ( g0 > d_g1 .and. g0 <= d_g2 .and. t < d_dtinj ) then

            injec_cond: if ( t < d_dtinj .and. tt < d_dtinj ) then

               ghigh = g0
               if ( g < d_gg(k) ) then
                  glow = d_gg(k)
                  t = t + (g0 - d_gg(k)) / (d_nu0(i, k) * g0 * d_gg(k))
                  nn = nn + d_Q0(k) * ghigh**(1d0 - d_qind(k)) * Pinteg(ghigh / glow, 2d0 - d_qind(k), 1d-9) / (d_nu0(i, k) * g0**2)
                  g0 = d_gg(k)
                  k = k - 1
               else
                  glow = g
                  nn = nn + d_Q0(k) * ghigh**(1d0 - d_qind(k)) * Pinteg(ghigh / glow, 2d0 - d_qind(k), 1d-9) / (d_nu0(i, k) * g0**2)
                  g0 = g
                  t = tt
                  if ( t >= d_t(i + 1) ) i = i + 1
               end if

            else if ( t < d_dtinj .and. tt >= d_dtinj ) then

               g = g0 / (1d0 + d_nu0(i, k) * (d_dtinj - t) * g0)
               ghigh = g0

               if ( g < d_gg(k) ) then
                  glow = d_gg(k)
                  t = t + (g0 - d_gg(k)) / (d_nu0(i, k) * g0 * d_gg(k))
                  nn = nn + d_Q0(k) * ghigh**(1d0 - d_qind(k)) * Pinteg(ghigh / glow, 2d0 - d_qind(k), 1d-9) / (d_nu0(i, k) * g0**2)
                  g0 = d_gg(k)
                  k = k - 1
               else
                  glow = dmax1(g, g / ( 1d0 + d_nu0(i, k) * (d_dtinj - t) * g ))
                  nn = nn + d_Q0(k) * ghigh**(1d0 - d_qind(k)) * Pinteg(ghigh / glow, 2d0 - d_qind(k), 1d-9) / (d_nu0(i, k) * g0**2)
                  t = dmin1( d_dtinj, t + d_tesc)
                  g0 = g
                  if ( t >= d_dtinj ) i = i + 1
               end if

            end if injec_cond

            tt = dmin1( d_t(i + 1), t + d_tesc )
            injecting = .true.

         else

            if ( injecting ) then
               t0 = t
               injecting = .false.
            end if

            if ( g < d_gg(k) .and. k > 1 ) then
               t = t + (g0 - d_gg(k)) / (d_nu0(i, k) * g0 * d_gg(k))
               nn = nn / (1d0 - d_nu0(i, k) * (tt - t) * d_gg(k))**2
               g0 = d_gg(k)
               k = k - 1
            else
               nn = nn / (1d0 - d_nu0(i, k) * (tt - t) * g)**2
               t = tt
               g0 = g
               if ( t >= t0 + d_tesc ) then
                  exit contrib_loop
               else
                  i = i + 1
               end if
            end if

            tt = dmin1( d_t(i + 1), t0 + d_tesc )

         end if

      end do contrib_loop

   end function distrib


end module dist_evol
