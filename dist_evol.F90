module dist_evol
   use data_types
   use constants
   use misc
   use pwl_integ
   use SRtoolkit
   use K2
   implicit none

   ! integer :: d_Ng, d_Nt
   ! real(dp) :: d_dtinj, d_g1, d_g2, d_tcurr, d_tesc
   ! real(dp), allocatable, dimension(:) :: d_gg, d_qind, d_Q0, d_t
   ! real(dp), allocatable, dimension(:, :) :: d_nu0, d_dist
   ! public :: d_ng, d_t, d_dtinj, d_g1, d_g2, d_tcurr, d_qind, &
   !    d_Q0, d_Nt, d_dist, d_tesc

   interface RMaxwell
      module procedure RMaxwell_s
      module procedure RMaxwell_v
   end interface RMaxwell

   interface powlaw_dis
      module procedure powlaw_dis_s
      module procedure powlaw_dis_v
   end interface powlaw_dis

contains

      ! #####  #  ####  ##### #####  # #####  #    # ##### #  ####  #    #  ####
      ! #    # # #        #   #    # # #    # #    #   #   # #    # ##   # #
      ! #    # #  ####    #   #    # # #####  #    #   #   # #    # # #  #  ####
      ! #    # #      #   #   #####  # #    # #    #   #   # #    # #  # #      #
      ! #    # # #    #   #   #   #  # #    # #    #   #   # #    # #   ## #    #
      ! #####  #  ####    #   #    # # #####   ####    #   #  ####  #    #  ####
      !
      !  :::::   Relativistic Maxwell distribution   :::::
      !
      function RMaxwell_s(g, th) result(rm)
         implicit none
         real(dp), intent(in) :: g, th
         real(dp) :: rm
         call K2_init
         rm = dmax1(1d-200, g**2 * bofg(g) * dexp(-g / th) / (th * dexp(K2_func(-dlog(th)))))
      end function RMaxwell_s

      function RMaxwell_v(g, th) result(rm)
         implicit none
         doubleprecision, intent(in) :: th
         real(dp), intent(in), dimension(:) :: g
         integer :: k
         real(dp), dimension(size(g)) :: rm
         call K2_init
         do k=1,size(g)
            rm(k) = dmax1( 1d-200, g(k)**2 * bofg(g(k)) * dexp(-g(k) / th) &
            / (th * dexp(K2_func(-dlog(th)))) )
         enddo
      end function RMaxwell_v

      !
      !   :::::   Power-law distribution   :::::
      !
      function powlaw_dis_s(g, g1, g2, s) result(pwl)
         implicit none
         doubleprecision, intent(in) :: g,g1,g2,s
         doubleprecision :: pwl
         if ( g >= g1 .and. g <=g2 ) then
            pwl = g**(-s)
         else
            pwl = 0d0
         end if
         ! pwl = g**(-s) * dexp(-(g1 / g)**2) * dexp(-(g / g2)**2)
      end function powlaw_dis_s

      function powlaw_dis_v(g, g1, g2, s) result(pwl)
         implicit none
         doubleprecision, intent(in) :: g1,g2,s
         doubleprecision, intent(in), dimension(:) :: g
         integer :: k
         doubleprecision, dimension(size(g)) :: pwl
         do k = 1, size(g)
            if ( g(k) >= g1 .and. g(k) <= g2 ) then
               pwl(k) = g(k)**(-s)
            else
               pwl(k) = 0d0
            end if
         end do
         ! pwl = g**(-s) * dexp(-(g1 / g)**2) * dexp(-(g / g2)**2)
      end function powlaw_dis_v


   ! # #    #      # ######  ####  ##### #  ####  #    #
   ! # ##   #      # #      #    #   #   # #    # ##   #
   ! # # #  #      # #####  #        #   # #    # # #  #
   ! # #  # #      # #      #        #   # #    # #  # #
   ! # #   ## #    # #      #    #   #   # #    # #   ##
   ! # #    #  ####  ######  ####    #   #  ####  #    #
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
   !
   !   ----------{    Power-law distribution Normalization    }----------
   !
   function pwl_norm(norm, index, e1, e2) result(k)
      implicit none
      real(dp), intent(in) :: norm, index, e1, e2
      real(dp), parameter :: eps = 1d-6
      real(dp) :: k, integ
      integ = e1**(1d0 - index) * Pinteg(e2 / e1, index, eps)
      k = 1d0 / (norm * integ)
   end function pwl_norm


   !  #######          ######
   !  #       # #    # #     # # ######
   !  #       # ##   # #     # # #
   !  #####   # # #  # #     # # #####
   !  #       # #  # # #     # # #
   !  #       # #   ## #     # # #
   !  #       # #    # ######  # #
   subroutine FP_FinDif_difu(dt, g, nin, nout, nu0, DD, QQ, AA, tesc)
      implicit none
      real(dp), intent(in) :: dt, tesc, AA
      real(dp), intent(in), dimension(:) :: g, nin, nu0, DD, QQ
      real(dp), intent(out), dimension(:) :: nout
      real(dp), parameter :: eps = 1e-3
      integer :: i, Ng, Ng1
      real(dp) :: dBB
      real(dp), dimension(size(g)) :: gdot, dx, dxp2, dxm2, CCp2, CCm2, BBp2, BBm2, YYp2, YYm2, WWp2, WWm2, ZZp2, ZZm2, a, b, c, r, x

      Ng = size(g)
      Ng1 = Ng - 1

      gdot = - nu0 * g**2 - AA * g
      x = g - 1d0
      dxp2(:Ng1) = x(2:) - x(:Ng1)
      dxp2(Ng) = dxp2(Ng1)
      dxm2(2:) = dxp2(:Ng1)
      dxm2(1) = dxm2(2)
      dx = dsqrt(dxp2 * dxm2)

      CCp2(:Ng1) = 0.25d0 * (DD(2:) + DD(:Ng1))
      CCm2(2:) = CCp2(:Ng1)
      CCp2(Ng) = 0.25d0 * DD(Ng)
      CCm2(1) = 0.25d0 * DD(1)

      BBp2(:Ng1) = -0.5d0 * ((DD(2:) - DD(:Ng1)) / dxp2(:Ng1) + (gdot(2:) + gdot(:Ng1)))
      BBm2(2:) = -0.5d0 * ((DD(2:) - DD(:Ng1)) / dxm2(2:) + (gdot(2:) + gdot(:Ng1)))
      BBp2(Ng) = -0.5d0 * ((DD(Ng) - DD(Ng1)) / dxp2(Ng) + (gdot(Ng) + gdot(Ng1)))
      call polint(x(2:), BBm2(2:), x(1), BBm2(1), dBB)

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

      r = nin + dt * QQ
      c = -dt * CCp2 * YYp2 * dexp(ZZp2)/ (dx * dxp2)
      a = -dt * CCm2 * YYm2 * dexp(-ZZm2) / (dx * dxm2)
      b = 1d0 + dt * ( CCp2 * YYp2 * dexp(-ZZp2) / dxp2 + CCm2 * YYm2 * dexp(ZZm2) / dxm2 ) / dx + dt / tesc

      call tridag_ser(a(2:), b, c(:Ng1), r, nout)

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
      dxm2(1) = dxm2(2)
      dx = dsqrt(dxp2 * dxm2)

      BBp2(:Ng1) = 0.5d0 * (gdot(2:) + gdot(:Ng1))
      BBm2(2:) = BBp2(:Ng1)
      BBp2(Ng) = 0.5d0 * (gdot(Ng) + gdot(Ng1))
      BBm2(1) = 0.5d0 * (gdot(2) + gdot(1))

      r = nin + dt * QQ
      a = zeros1D(Ng)
      c = -dt * BBp2 / dx
      b = 1d0 + BBm2 * dt / dx + dt / tesc

      call tridag_ser(a(2:), b, c(:Ng1), r, nout)

   end subroutine FP_FinDif_cool


   subroutine time_step(dt, g, D, gdot, tstep, lct)
      implicit none
      real(dp), intent(in) :: tstep, lct
      real(dp), intent(in), dimension(:) :: g, D, gdot
      real(dp), intent(out) :: dt
      integer :: Ng
      real(dp), dimension(size(g)) :: H, dg
      Ng = size(g)
      dg(2:) = g(2:) - g(:Ng - 1)
      dg(1) = dg(2)
      H(2:) = dmax1(1d-200, dabs(D(2:) - D(:Ng - 1)) / dg(2:)) + dabs(gdot(2:))
      H(1) = dmax1(1d-200, dabs(D(2) - D(1)) / dg(1)) + dabs(gdot(1))
      dt = dmax1(tstep * lct, 0.5d0 * minval(dg / H, mask=H>1d-50))
   end subroutine time_step






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
   !
   !     Slope of the Relativistic Maxwell distribution
   !
   function dRMaxwell(g, th) result(drm)
      implicit none
      real(dp), intent(in) :: g, th
      real(dp) :: drm
      call K2_init
      !!$drm = -((2d0 * g**(2) - 1d0) / (g**(2) - 1d0) - g / th)
      drm = (g - g**3 - th + 2 * th * g**2) * dexp(-g/th) / ( dsqrt(g**2 - 1d0) * th**2 * dexp(K2_func(-dlog(th))) )
   end function dRMaxwell


   subroutine FP_FinDif_difu2(dt, g, nin, nout, nu0, DD, QQ, tesc)
      implicit none
      real(dp), intent(in) :: dt, tesc
      real(dp), intent(in), dimension(:) :: g, nin, nu0, DD, QQ
      real(dp), intent(out), dimension(:) :: nout
      real(dp), parameter :: eps = 1e-6
      integer :: i, Ng, Ng1
      real(dp), dimension(size(g)) :: gdot, dx, dxp2, dxm2, a, b, c, r, &
         CCp2, CCm2, BBp2, BBm2, YYp2, YYm2, XXp2, XXm2, WWp2, WWm2, ZZp2, ZZm2

      Ng = size(g)
      Ng1 = Ng - 1

      gdot = nu0 * g**2
      dxp2(:Ng1) = g(2:) - g(:Ng1)
      dxp2(Ng) = dxp2(Ng1)
      dxm2(2:) = dxp2(:Ng1)
      dxm2(1) = dmin1(g(1) - 1d0, dxm2(2))
      dx = dsqrt(dxp2 * dxm2)

      CCp2(:Ng1) = 0.25d0 * (DD(2:) + DD(:Ng1))
      CCm2(2:) = CCp2(:Ng1)
      CCp2(Ng) = 0.25d0 * DD(Ng)
      CCm2(1) = 0.25d0 * DD(1)

      BBp2(:Ng1) = 0.5d0 * ( (DD(2:) - DD(:Ng1)) / dxp2(:Ng1) - (gdot(2:) + gdot(:Ng1)) )
      BBm2(2:) = 0.5d0 * ( (DD(2:) - DD(:Ng1)) / dxm2(2:) - (gdot(2:) + gdot(:Ng1)) )
      BBp2(Ng) = 0d0 !0.5d0 * ( -DD(Ng) / dxp2(Ng) - gdot(Ng) )
      BBm2(1) = 0d0 !0.5d0 * ( DD(1) / dxm2(1) - gdot(1) )

      WWp2(:Ng1) = dxp2(:Ng1) * BBp2(:Ng1) / CCp2(:Ng1)
      WWm2(2:) = dxm2(2:) * BBm2(2:) / CCm2(2:)
      WWp2(Ng) = 0d0
      WWm2(1) = 0d0

      do i = 1, Ng

         if ( 0.5d0 * WWp2(i) > -lzero ) then
            XXp2(i) = dexp(-lzero)
            YYp2(i) = dexp(lzero)
         else if ( 0.5d0 * WWp2(i) < lzero ) then
            XXp2(i) = dexp(lzero)
            YYp2(i) = dexp(-lzero)
         else
            XXp2(i) = dexp(0.5d0 * WWp2(i))
            YYp2(i) = dexp(-0.5d0 * WWp2(i))
         end if

         if ( 0.5d0 * WWm2(i) > -lzero ) then
            XXm2(i) = dexp(-lzero)
            YYm2(i) = dexp(lzero)
         else if ( 0.5d0 * WWm2(i) < lzero ) then
            XXm2(i) = dexp(lzero)
            YYm2(i) = dexp(-lzero)
         else
            XXm2(i) = dexp(0.5d0 * WWm2(i))
            YYm2(i) = dexp(-0.5d0 * WWm2(i))
         end if

         if ( 127d0 * WWp2(i)**8 < eps * 154828800d0 ) then
            ZZp2(i) = 1d0 - WWp2(i)**2 / 24d0 + 7d0 * WWp2(i)**4 / 5760d0 - 31d0 * WWp2(i)**6 / 967680d0
         else
            ZZp2(i) = 0.5d0 * WWp2(i) / dsinh(0.5d0 * WWp2(i))
         end if

         if ( 127d0 * WWm2(i)**8 < eps * 154828800d0 ) then
            ZZm2(i) = 1d0 - WWm2(i)**2 / 24d0 + 7d0 * WWm2(i)**4 / 5760d0 - 31d0 * WWm2(i)**6 / 967680d0
         else
            ZZm2(i) = 0.5d0 * WWm2(i) / dsinh(0.5d0 * WWm2(i))
         end if
      
      end do

      r = nin + dt * QQ
      a = -dt * CCm2 * ZZm2 * YYm2 / (dx * dxm2)
      c = -dt * CCp2 * ZZp2 * XXp2 / (dx * dxp2)
      b = 1d0 + dt * ( CCp2 * ZZp2 * YYp2 / dxp2 + CCm2 * ZZm2 * XXm2 / dxm2 ) / dx !+ dt / tesc

      call tridag_ser(a(2:), b, c(:Ng1), r, nout)

   end subroutine FP_FinDif_difu2

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
#endif

end module dist_evol
