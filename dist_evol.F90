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
   function injection(t,dtinj,gg,g1,g2,qind,th,Qth,Qnth) result(Qinj)
      implicit none
      real(dp), intent(in) :: g1,g2,th,Qth,Qnth,dtinj,t,qind
      real(dp), intent(in), dimension(:) :: gg
      integer :: k
      real(dp), dimension(size(gg)) :: Qinj
      if ( t <= dtinj ) then
         if ( Qth < 1d-100 ) then
            Qinj = dmax1(1d-200, Qnth * powlaw_dis(gg,g1,g2,qind))
         else if ( Qnth < 1d-100 ) then
            Qinj = dmax1(1d-200, Qth * RMaxwell(gg, th))
         else
            do k = 1, size(gg)
               if ( gg(k) >= g1 .and. gg(k) <= g2 ) then
                  Qinj(k) = dmax1(1d-200, Qnth * powlaw_dis(gg(k),g1,g2,qind) + Qth * RMaxwell(gg(k), th))
               else
                  Qinj(k) = dmax1(1d-200, Qth * RMaxwell(gg(k), th))
               end if
            end do
         end if
      else
         Qinj = 1d-200
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
   subroutine cooling_lines(nn,gg,nu0,times,dtinj,Qinj)
      real(dp), intent(in) :: dtinj
      real(dp), intent(in), dimension(:) :: gg,times
      real(dp), intent(in), dimension(:,:) :: nu0,Qinj
      real(dp), intent(out), dimension(:) :: nn
      integer :: k,i,kk,Ng,Nt
      double precision :: g0,ghigh,glow,g,tup,tdw,t_curr
      real(dp), dimension(size(times), size(gg)) :: q, Q0

      Ng = ubound(gg, dim=1)
      Nt = ubound(times, dim=1)
      t_curr = times(Nt)
      do i = 1, Nt

         do kk = 1, Ng - 1
            if ( Qinj(i, kk + 1) > 1d-100 .and. Qinj(i, kk) > 1d-100 ) then
               q(i, kk) = -dlog(Qinj(i, kk + 1) / Qinj(i, kk)) / dlog(gg(kk + 1) / gg(kk))
               if ( q(i, kk) > 8d0 ) q(i, kk) = 8d0
               if ( q(i, kk) < -8d0 ) q(i, kk) = -8d0
               Q0(i, kk) = Qinj(i, kk) * gg(kk)**q(i, kk)
            else
               q(i, kk) = 0d0
               Q0(i, kk) = 0d0
            endif
         enddo

         if ( Qinj(i, Ng - 1) > 1d-100 .and. Qinj(i, Ng) > 1d-100 ) then
            q(i, Ng) = -dlog(Qinj(i, Ng) / Qinj(i, Ng - 1)) / dlog(gg(Ng) / gg(Ng - 1))
            if ( q(i, Ng) > 8d0 ) q(i, Ng) = 8d0
            if ( q(i, Ng) < -8d0 ) q(i, Ng) = -8d0
            Q0(i, Ng) = Qinj(i, Ng) * gg(Ng)**q(i, Ng)
         else
            q(i, Ng) = 0d0
            Q0(i, Ng) = 0d0
         endif

      enddo

      !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(STATIC) DEFAULT(SHARED) &
      !$OMP& PRIVATE(i,k,tup,tdw,glow,ghigh,g,g0)
      gloop: do kk = 1, Ng
         !
         !  #    # #####  #    #   ##   #####  #####   ####
         !  #    # #    # #    #  #  #  #    # #    # #
         !  #    # #    # #    # #    # #    # #    #  ####
         !  #    # #####  # ## # ###### #####  #    #      #
         !  #    # #      ##  ## #    # #   #  #    # #    #
         !   ####  #      #    # #    # #    # #####   ####
         !
         tdw= t_curr
         g = gg(kk)
         i = Nt - 1
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
         i = locate(times, tup, .true.)
         ! k = locate(gg, g0, .true.)
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
