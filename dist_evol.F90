module dist_evol
   use data_types
   use constants
   use misc
   use pwl_integ
   use SRtoolkit
   implicit none


contains

   !
   ! # #    #      # ######  ####  ##### #  ####  #    #
   ! # ##   #      # #      #    #   #   # #    # ##   #
   ! # # #  #      # #####  #        #   # #    # # #  #
   ! # #  # #      # #      #        #   # #    # #  # #
   ! # #   ## #    # #      #    #   #   # #    # #   ##
   ! # #    #  ####  ######  ####    #   #  ####  #    #
   !
   function injection(t, dtinj, Ng, gmin, gmax, g1, g2, qind, th, Qth, Qnth) result(Qinj)
      implicit none
      integer, intent(in) :: Ng
      real(dp), intent(in) :: gmin, gmax, g1, g2, th, Qth, Qnth, dtinj, t, qind
      real(dp) :: g
      integer :: k
      real(dp), dimension(Ng) :: Qinj
      if ( t <= dtinj ) then
         do k = 1, Ng
            g = gmin * (gmax / gmin)**(dble(k - 1) / dble(Ng - 1))
            if ( Qth < 1d-100 ) then
               Qinj(k) = Qnth * powlaw_dis(g, g1, g2, qind)
            else if ( Qnth < 1d-100 ) then
               Qinj(k) = Qth * RMaxwell(g, th)
            else
               if ( g >= g1 .and. g <= g2 ) then
                  Qinj(k) = Qnth * powlaw_dis(g, g1, g2, qind) + Qth * RMaxwell(g, th)
               else
                  Qinj(k) = Qth * RMaxwell(g, th)
               end if
            end if
         end do
      else
         Qinj = 0d0
      end if
   end function injection


   !
   !   ####   ####   ####  #            #      # #    # ######  ####
   !  #    # #    # #    # #            #      # ##   # #      #
   !  #      #    # #    # #      ##### #      # # #  # #####   ####
   !  #      #    # #    # #            #      # #  # # #           #
   !  #    # #    # #    # #            #      # #   ## #      #    #
   !   ####   ####   ####  ######       ###### # #    # ######  ####
   !
   subroutine distrib_setup(nn, gg, nu0, t, dtinj, tesc, Ng, gmin, gmax, g1, g2, qind, theta_e, Qth, Qnth)
      implicit none
      integer, intent(in) :: Ng
      real(dp), intent(in) :: dtinj, tesc, g1, g2, theta_e, Qth, Qnth, qind, &
         gmin, gmax
      real(dp), intent(in), dimension(:) :: t
      real(dp), intent(in), dimension(:) :: nu0
      real(dp), intent(inout), dimension(:, :) :: nn, gg
      integer :: k, kk
      double precision :: g0, ghigh, glow, g, tup, tdw, t_curr
      real(dp), dimension(Ng) :: q, Qinj, Q0

      if ( t(1) <= 0d0 ) then
         gg(2, :) = gg(1, :)
         nn(2, :) = t(2) * injection(t(2), dtinj, Ng, gmin, gmax, g1, g2, qind, theta_e, Qth, Qnth)
         return
      end if

      Qinj = injection(t(2), dtinj, Ng, gmin, gmax, g1, g2, qind, theta_e, Qth, Qnth)
      do k = 1, Ng - 1
         if ( Qinj(k + 1) > 1d-100 .and. Qinj(k) > 1d-100 ) then
            q(k) = -dlog(Qinj(k + 1) / Qinj(k)) / dlog(gg(2, k + 1) / gg(2, k))
            if ( q(k) > 8d0 ) q(k) = 8d0
            if ( q(k) < -8d0 ) q(k) = -8d0
            Q0(k) = Qinj(k) * gg(2, k)**q(k)
         else
            q(k) = 0d0
            Q0(k) = 0d0
         endif
      enddo
      q(Ng) = q(Ng - 1)
      Q0(Ng) = Qinj(Ng) * gg(2, Ng)**q(Ng)

   end subroutine distrib_setup


   function distrib(g)
      implicit none
      
      ! t_curr = t(2)
      t_curr = dmin1( t(2), t(1) + tesc )
      
      gloop: do k = 1, Ng
         nn(2, k) = nn(1, k)
         kk = max0(k - 1, 1)
         tup = t(1)
         contrib_loop: do while ( tup < t_curr )
            
            tdw = t_curr
            
            if ( g0 > g1 .and. g0 <= g2 .and. tup < dtinj ) then
               
               injection: if ( tup < dtinj .and. t_curr < dtinj ) then
                  g = g0 / (1d0 + nu0(kk) * (tdw - tup) * g0)
                  ghigh = g0
                  if ( g < gg(1, kk) ) then
                     glow = gg(1, kk)
                     tdw = tup + (g0 - gg(1, kk)) / (nu0(kk) * g0 * gg(1, kk))
                     nn(2, k) = nn(2, k) + Q0(kk) * ghigh**(1d0 - q(kk)) * Pinteg(ghigh / glow, 2d0 - q(kk), 1d-9) / (nu0(kk) * g0**2)
                     g0 = gg(1, kk)
                     g = gg(1, kk)
                     kk = kk - 1
                     tup = tdw
                  else
                     glow = g
                     nn(2, k) = nn(2, k) + Q0(kk) * ghigh**(1d0 - q(kk)) * Pinteg(ghigh / glow, 2d0 - q(kk), 1d-9) / (nu0(kk) * g0**2)
                     tup = t_curr
                     g0 = g
                  end if
               else if ( tup < dtinj .and. tdw >= dtinj ) then
                  g = g0 / (1d0 + nu0(kk) * (dtinj - tup) * g0)
                  ghigh = g0
                  if ( g < gg(1, kk) ) then
                     glow = gg(1, kk)
                     tdw = tup + (g0 - gg(1, kk)) / (nu0(kk) * g0 * gg(1, kk))
                     nn(2, k) = nn(2, k) + Q0(kk) * ghigh**(1d0 - q(kk)) * Pinteg(ghigh / glow, 2d0 - q(kk), 1d-9) / (nu0(kk) * g0**2)
                     g0 = gg(1, kk)
                     g = gg(1, kk)
                     kk = kk - 1
                     tup = tdw
                  else
                     glow = g / ( 1d0 + nu0(kk) * (dtinj - tup) * g )
                     nn(2, k) = nn(2, k) + Q0(kk) * ghigh**(1d0 - q(kk)) * Pinteg(ghigh / glow, 2d0 - q(kk), 1d-9) / (nu0(kk) * g0**2)
                     tup = dtinj
                     g0 = g
                  end if
               end if injection
               
            else
               
               g0 = gg(1, k)
               g = g0 / (1d0 + nu0(kk) * (tdw - tup) * g0)
               
               if ( g < gg(1, kk) .and. kk > 1 ) then
                  tdw = tup + (g0 - gg(1, kk)) / (nu0(kk) * g0 * gg(1, kk))
                  g = gg(1, kk)
                  nn(2, k) = nn(2, k) / (1d0 - nu0(kk) * (tdw - tup) * g)**2
                  g0 = gg(1, kk)
                  kk = kk - 1
                  tup = tdw
               else
                  nn(2, k) = nn(2, k) / (1d0 - nu0(kk) * (tdw - tup) * g)**2
                  tup = t_curr
                  g0 = g
               end if
               
               gg(2, k) = g
               
            end if
            
         end do contrib_loop
      end do gloop
      
   end function distrib


end module dist_evol
