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
   function injection(t, dtinj, gg, g1, g2, qind, th, Qth, Qnth) result(Qinj)
      implicit none
      real(dp), intent(in) :: g1,g2,th,Qth,Qnth,dtinj,t,qind
      real(dp), intent(in), dimension(:) :: gg
      integer :: k
      real(dp), dimension(size(gg)) :: Qinj
      if ( t <= dtinj ) then
         if ( Qth < 1d-100 ) then
            Qinj = Qnth * powlaw_dis(gg,g1,g2,qind)
         else if ( Qnth < 1d-100 ) then
            Qinj = Qth * RMaxwell(gg, th)
         else
            do k = 1, size(gg)
               if ( gg(k) >= g1 .and. gg(k) <= g2 ) then
                  Qinj(k) = Qnth * powlaw_dis(gg(k),g1,g2,qind) + Qth * RMaxwell(gg(k), th)
               else
                  Qinj(k) = Qth * RMaxwell(gg(k), th)
               end if
            end do
         end if
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
   subroutine cooling_lines(nn, gg, nu0, t, dtinj, tesc, g1, g2, qind, theta_e, Qth, Qnth)
      real(dp), intent(in) :: dtinj, tesc, g1, g2, theta_e, Qth, Qnth, qind
      real(dp), intent(in), dimension(:) :: t
      real(dp), intent(in), dimension(:) :: nu0
      real(dp), intent(inout), dimension(:, :) :: gg, nn
      integer :: k, kk, Ng
      double precision :: g0, ghigh, glow, g, tup, tdw, t_curr
      real(dp), dimension(size(gg, dim=2)) :: q, Qinj, Q0
      Ng = size(nu0, dim=1)
      if ( t(1) <= 0d0 ) then
         gg(2, :) = gg(1, :) / (1d0 + nu0 * (t(2) - t(1)) * gg(1, :))
         nn(2, :) = t(2) * injection(t(2), dtinj, gg(2, :), g1, g2, qind, theta_e, Qth, Qnth)
         return
      end if

      Qinj = injection(t(1), dtinj, gg(1, :), g1, g2, qind, theta_e, Qth, Qnth)
      do kk = 1, Ng - 1
         if ( Qinj(kk + 1) > 1d-100 .and. Qinj(kk) > 1d-100 ) then
            q(kk) = -dlog(Qinj(kk + 1) / Qinj(kk)) / dlog(gg(1, kk + 1) / gg(1, kk))
            if ( q(kk) > 8d0 ) q(kk) = 8d0
            if ( q(kk) < -8d0 ) q(kk) = -8d0
            Q0(kk) = Qinj(kk) * gg(1, kk)**q(kk)
         else
            q(kk) = 0d0
            Q0(kk) = 0d0
         endif
      enddo
      q(Ng) = q(Ng - 1)
      Q0(Ng) = Qinj(Ng) * gg(1, Ng)**q(Ng)

      t_curr = t(2)

      gloop: do k = 1, Ng
         nn(2, k) = nn(1, k)
         kk = max0(k - 1, 1)
         g0 = gg(1, k)
         g = g0
         tup = t(1)
         contrib_loop: do while ( tup < t_curr .and. kk >= 1 )

            tdw = t_curr

            injection: if ( tup < dtinj .and. t_curr < dtinj ) then

               g = g0 / (1d0 + nu0(kk) * (tdw - tup) * g0)
               ghigh = g0
               if ( g < gg(1, kk) ) then
                  glow = gg(1, kk)
                  tdw = tup + (g0 - gg(1, kk)) / (nu0(kk) * g0 * gg(1, kk))
                  ! nn(2, k) = nn(2, k) * (1d0 + nu0(kk) * (tdw - tup) * g0)**2
                  nn(2, k) = nn(2, k) + Q0(kk) * ghigh**(1d0 - q(kk)) * Pinteg(ghigh / glow, 2d0 - q(kk), 1d-9) / (nu0(kk) * g0**2)
                  g0 = gg(1, kk)
                  g = gg(1, kk)
                  kk = kk - 1
                  tup = tdw
               else
                  glow = g
                  ! nn(2, k) = nn(2, k) * (1d0 + nu0(kk) * (tdw - tup) * g0)**2
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
                  ! nn(2, k) = nn(2, k) * (1d0 + nu0(kk) * (tdw - tup) * g0)**2
                  nn(2, k) = nn(2, k) + Q0(kk) * ghigh**(1d0 - q(kk)) * Pinteg(ghigh / glow, 2d0 - q(kk), 1d-9) / (nu0(kk) * g0**2)
                  g0 = gg(1, kk)
                  g = gg(1, kk)
                  kk = kk - 1
                  tup = tdw
               else
                  glow = g / ( 1d0 + nu0(kk) * (dtinj - tup) * g )
                  ! nn(2, k) = nn(2, k) * (1d0 + nu0(kk) * (dtinj - tup) * g0)**2
                  nn(2, k) = nn(2, k) + Q0(kk) * ghigh**(1d0 - q(kk)) * Pinteg(ghigh / glow, 2d0 - q(kk), 1d-9) / (nu0(kk) * g0**2)
                  tup = dtinj
                  g0 = g
               end if

            else

               g = g0 / (1d0 + nu0(kk) * (tdw - tup) * g0)
               if ( g < gg(1, kk) ) then
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

            end if injection

            if ( g * tesc * nu0(kk) < 1d0 ) exit contrib_loop

         enddo contrib_loop
         gg(2, k) = g
      end do gloop

   end subroutine cooling_lines

end module dist_evol
