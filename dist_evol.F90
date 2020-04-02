module dist_evol
   use data_types
   use constants
   use misc
   use pwl_integ
   use SRtoolkit
   use K2
   implicit none

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
         doubleprecision, intent(in) :: g, g1, g2, s
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
   function injection_pwl(t, dtinj, g, g1, g2, qind, Q0) result(Qinj)
      implicit none
      real(dp), intent(in) :: g1, g2, Q0, dtinj, t, qind
      real(dp), intent(in), dimension(:) :: g
      real(dp), dimension(size(g)) :: Qinj
      if ( t <= dtinj ) then
         Qinj = Q0 * powlaw_dis(g, g1, g2, qind)
      else
         Qinj = 0d0
      end if
   end function injection_pwl


   function injection_hyb(t, dtinj, g, g1, g2, qind, th, Qth, Qnth) result(Qinj)
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
   end function injection_hyb
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
   subroutine FP_FinDif_difu(dt_in, g, nin, nout, gdot_in, Din, Qin, tesc_in, tlc)
      implicit none
      real(dp), intent(in) :: dt_in, tesc_in, tlc
      real(dp), intent(in), dimension(:) :: g, nin, gdot_in, Din, Qin
      real(dp), intent(out), dimension(:) :: nout
      real(dp), parameter :: eps = 1e-3
      integer :: i, Ng, Ng1
      real(dp) :: dBB, dt, tesc
      real(dp), dimension(size(g)) :: dx, dxp2, dxm2, CCp2, CCm2, BBp2, BBm2, &
         YYp2, YYm2, WWp2, WWm2, ZZp2, ZZm2, a, b, c, r, DD, QQ, gdot

      Ng = size(g)
      Ng1 = Ng - 1
      !-----> Dimensionless
      dt = dt_in / tlc
      if ( tesc_in < 1d100 .and. tesc_in > 1d-100 ) then
         tesc = tesc_in / tlc
      else
         tesc = tesc_in
      end if
      gdot = gdot_in * tlc
      DD = Din * tlc
      QQ = Qin * tlc

      dxp2(:Ng1) = g(2:) - g(:Ng1)
      dxp2(Ng) = dxp2(Ng1)
      dxm2(2:) = dxp2(:Ng1)
      dxm2(1) = dxm2(2)
      dx = dsqrt(dxp2 * dxm2)

      CCp2(:Ng1) = 0.25d0 * (DD(2:) + DD(:Ng1))
      CCm2(2:) = CCp2(:Ng1)
      CCp2(Ng) = 0.25d0 * DD(Ng)
      CCm2(1) = 0.25d0 * DD(1)

      BBp2(:Ng1) = 0.5d0 * ((DD(2:) - DD(:Ng1)) / dxp2(:Ng1) + (gdot(2:) + gdot(:Ng1)))
      BBm2(2:) = 0.5d0 * ((DD(2:) - DD(:Ng1)) / dxm2(2:) + (gdot(2:) + gdot(:Ng1)))
      BBp2(Ng) = 0.5d0 * ((DD(Ng) - DD(Ng1)) / dxp2(Ng) + (gdot(Ng) + gdot(Ng1)))
      call polint(g(2:), BBm2(2:), g(1), BBm2(1), dBB)

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
      c = -dt * CCp2 * YYp2 * dexp(ZZp2) / (dx * dxp2)
      a = -dt * CCm2 * YYm2 * dexp(-ZZm2) / (dx * dxm2)
      b = 1d0 + dt * ( CCp2 * YYp2 * dexp(-ZZp2) / dxp2 + CCm2 * YYm2 * dexp(ZZm2) / dxm2 ) / dx + dt / tesc
      call tridag_ser(a(2:), b, c(:Ng1), r, nout)

      nout = nout

   end subroutine FP_FinDif_difu


   subroutine FP_FinDif_cool(dt, g, nin, nout, gdot, QQ, tesc)
      implicit none
      real(dp), intent(in) :: dt, tesc
      real(dp), intent(in), dimension(:) :: g, nin, gdot, QQ
      real(dp), intent(out), dimension(:) :: nout
      integer :: Ng, Ng1
      real(dp), dimension(size(g)) :: dx, dxp2, dxm2, BBp2, BBm2, a, b, c, r

      Ng = size(g)
      Ng1 = Ng - 1

      dxp2(:Ng1) = g(2:) - g(:Ng1)
      dxp2(Ng) = dxp2(Ng1)
      dxm2(2:) = dxp2(:Ng1)
      dxm2(1) = dxm2(2)
      dx = dsqrt(dxp2 * dxm2)

      BBp2(:Ng1) = -0.5d0 * (gdot(2:) + gdot(:Ng1))
      BBm2(2:) = BBp2(:Ng1)
      BBp2(Ng) = -0.5d0 * (gdot(Ng) + gdot(Ng1))
      BBm2(1) = -0.5d0 * (gdot(2) + gdot(1))

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
      real(dp) :: Hmin
      real(dp), dimension(size(g)) :: H, dg
      Ng = size(g)
      dg(2:) = g(2:) - g(:Ng - 1)
      dg(1) = dg(2)
      H(2:) = dmax1(1d-200, dabs(D(2:) - D(:Ng - 1)) / dg(2:)) + dabs(gdot(2:))
      H(1) = dmax1(1d-200, dabs(D(2) - D(1)) / dg(1)) + dabs(gdot(1))
      Hmin = minval(dg / H, mask=H>1d-50)
      if ( Hmin > lct ) then
         dt = tstep * lct
      else
         dt = dmax1(tstep * lct, 0.5d0 * Hmin)
      end if
   end subroutine time_step

end module dist_evol





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
#endif
