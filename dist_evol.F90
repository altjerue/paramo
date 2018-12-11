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
   subroutine cooling_lines(nn, gg, nu0, times, dtinj, tesc, g1, g2, theta_e, Qth, Qnth)
      real(dp), intent(in) :: dtinj, tesc, g1, g2, theta_e, Qth, Qnth
      real(dp), intent(in), dimension(:) :: times
      real(dp), intent(in), dimension(:,:) :: nu0, gg
      real(dp), intent(out), dimension(:) :: nn
      integer :: k, i, kk, Ng, Nt, iup, kup, istart
      double precision :: g0, ghigh, glow, g, tup, tdw, t_curr
      real(dp), dimension(size(gg, dim=2)) :: q, Q0

      Ng = size(gg, dim = 2)
      Nt = size(times, dim = 1)
      t_curr = dmin1(times(Nt), tesc)
      istart = minloc(t_curr - times, mask = t_curr >= times)
      Q0 = injection(t(istart), dtinj, gg(istart, :), g1, g2, qind, theta_e, Qth, Qnth)
      do kk = 1, Ng - 1
         if ( Qinj(kk + 1) > 1d-100 .and. Qinj(kk) > 1d-100 ) then
            q(kk) = -dlog(Qinj(kk + 1) / Qinj(kk)) / dlog(gg(kk + 1) / gg(kk))
            if ( q(kk) > 8d0 ) q(kk) = 8d0
            if ( q(kk) < -8d0 ) q(kk) = -8d0
            Q0(kk) = Qinj(kk) * gg(kk)**q(kk)
         else
            q(kk) = 0d0
            Q0(kk) = 0d0
         endif
      enddo

         if ( Qinj(Ng - 1) > 1d-100 .and. Qinj(Ng) > 1d-100 ) then
            q(Ng) = -dlog(Qinj(Ng) / Qinj(Ng - 1)) / dlog(gg(Ng) / gg(Ng - 1))
            if ( q(i, Ng) > 8d0 ) q(Ng) = 8d0
            if ( q(i, Ng) < -8d0 ) q(Ng) = -8d0
            Q0(Ng) = Qinj(Ng) * gg(Ng)**q(Ng)
         else
            q(Ng) = 0d0
            Q0(Ng) = 0d0
         endif
         
      enddo

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED) &
      !$OMP& PRIVATE(i, k, iup, kup, t_cool, tup, tdw, glow, ghigh, g, g0)
      gloop: do kk = 1, Ng
         !
         !  #    # #####  #    #   ##   #####  #####   ####
         !  #    # #    # #    #  #  #  #    # #    # #
         !  #    # #    # #    # #    # #    # #    #  ####
         !  #    # #####  # ## # ###### #####  #    #      #
         !  #    # #      ##  ## #    # #   #  #    # #    #
         !   ####  #      #    # #    # #    # #####   ####
         !
         tdw = t_curr
         i = i_end
         g = gg(kk)
         k = max(1, kk - 1)
         find_g0: do while( i >= 1 .and. k < Ng )

            tup = times(i)
            g0 = g / ( 1d0 - nu0(i, k) * (tdw - tup) * g )

            if ( g0 > gg(k + 1) ) then
               tdw = tup - (gg(k + 1) - g0) / (nu0(i, k) * g0 * gg(k + 1))
               g = gg(k + 1)
               k = k + 1
            else
               tdw = times(i)
               g = g0
               i = i - 1
            end if
            kup = k
         end do find_g0
         !
         !  #####   ####  #    # #    # #    #   ##   #####  #####   ####
         !  #    # #    # #    # ##   # #    #  #  #  #    # #    # #
         !  #    # #    # #    # # #  # #    # #    # #    # #    #  ####
         !  #    # #    # # ## # #  # # # ## # ###### #####  #    #      #
         !  #    # #    # ##  ## #   ## ##  ## #    # #   #  #    # #    #
         !  #####   ####  #    # #    # #    # #    # #    # #####   ####
         !
         nn(kk) = 0d0
         i = iup != locate(times, tup, .true.)
         k = kup ! k = locate(gg, g0, .true.)
         contrib_loop: do while ( i < Nt .and. g >= gg(1) )

            tdw = times(i + 1)

            injection: if ( tup < dtinj .and. tdw < dtinj ) then

               g = g0 / (1d0 + nu0(i, k) * (tdw - tup) * g0)

               ghigh = dmax1(1d0, g0)

               if ( g < gg(k) ) then

                  glow = gg(k)
                  tup = tup + (g0 - gg(k)) / (nu0(i, k) * g0 * gg(k))
                  g0 = gg(k)
                  nn(kk) = nn(kk) + Q0(i, k) * ghigh**(1d0 - q(i, k)) * & 
                  Pinteg(ghigh / glow, 2d0 - q(i, k), 1d-9) / (nu0(i, k) * g0**2)
                  k = k - 1

               else

                  glow = g
                  tup = times(i + 1)
                  g0 = g
                  nn(kk) = nn(kk) + Q0(i, k) * ghigh**(1d0 - q(i, k)) * Pinteg(ghigh / glow, 2d0 - q(i, k), 1d-9) / (nu0(i, k) * g0**2)
                  i = i + 1

               end if

            else if ( tup < dtinj .and. tdw >= dtinj ) then

               g = g0 / (1d0 + nu0(i, k) * (dtinj - tup) * g0)

               ghigh = dmax1(1d0, g0)

               if ( g < gg(k) ) then

                  glow = gg(k)
                  tup = tup + (g0 - gg(k)) / (nu0(i, k) * g0 * gg(k))
                  g0 = gg(k)
                  nn(kk) = nn(kk) + Q0(i, k) * ghigh**(1d0 - q(i, k)) * Pinteg(ghigh / glow, 2d0 - q(i, k), 1d-9) / (nu0(i, k) * g0**2)
                  k = k - 1

               else

                  glow = g / ( 1d0 + nu0(i, k) * (dtinj - tup) * g )
                  tup = dtinj
                  g0 = g
                  nn(kk) = nn(kk) + Q0(i, k) * ghigh**(1d0 - q(i, k)) * Pinteg(ghigh / glow, 2d0 - q(i, k), 1d-9) / (nu0(i, k) * g0**2)
                  i = locate(times, dtinj, .true.)

               end if

            else

               g = g0 / (1d0 + nu0(i, k) * (tdw - tup) * g0)

               if ( g < gg(k) ) then

                  tdw = tup + (g0 - gg(k)) / (g0 * nu0(i, k) * gg(k))
                  nn(kk) = nn(kk) / (1d0 - nu0(i, k) * (tdw - tup) * gg(k))**2
                  tup = tdw
                  g0 = gg(k)
                  k = k - 1

               else

                  nn(kk) = nn(kk) / (1d0 - nu0(i, k) * (tdw - tup) * g)**2
                  tup = times(i + 1)
                  g0 = g
                  i = i + 1

               end if

            end if injection

         enddo contrib_loop
         !nn(kk) = dmax1(1d-200, nn(kk))
      end do gloop
      !$OMP END PARALLEL DO

   end subroutine cooling_lines

end module dist_evol
