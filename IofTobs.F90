program IofTobs

   
   implicit none


   write(*, *) ''
   write(*, *) '--> Moving to the observer frame'

   Iobs = 0d0
   i_edge = minloc(abs(2d0 * R * mu_com - sen_lum), dim=1, mask=2d0 * R * mu_com < sen_lum)
   if ( i_edge == 0 ) i_edge = numdt + 1

   !$OMP PARALLEL DO ORDERED COLLAPSE(2) SCHEDULE(AUTO) DEFAULT(SHARED) &
   !$OMP& PRIVATE(tob_min, tob_max, ii, sind, i_start, at_the_edge)
   freq_loop: do j=1, numdf
      obs_loop: do i = 1, numdt
         at_the_edge = .false.
         if ( 2d0 * R * mu_com - sen_lum(i) < 0d0 ) at_the_edge = .true.
         if ( .not. at_the_edge ) then
            i_start = 1
         else
            i_start = i - i_edge
         end if
         int_loop: do ii = i_start, i
            tob_min = t_com_f(t_obs(i), z, gamma_bulk, t(ii - 1) * cspeed * mu_com, mu_obs)
            tob_max = t_com_f(t_obs(i), z, gamma_bulk, t(ii    ) * cspeed * mu_com, mu_obs)
            if ( ii == 1 ) then
               Iobs(j, i) = dabs(tob_max - tob_min) * jnut(j, 1)
            else
               if ( jnut(j, ii) > 1d-100 .and. jnut(j, ii - 1) > 1d-100 ) then
                  sind = -dlog(jnut(j, ii) / jnut(j, ii - 1)) / dlog(tob_max / tob_min)
                  if ( sind < -8d0 ) sind = -8d0
                  if ( sind > 8d0 ) sind = 8d0
                  Iobs(j, i) = Iobs(j, i) + jnut(j, ii - 1) * tob_min * ibpfunc(tob_max / tob_min, sind, 1d-6) / (gamma_bulk * mu_com * (mu_obs - bofg(gamma_bulk)) * D)
               end if
            end if
         end do int_loop
      end do obs_loop
   end do freq_loop
   !$OMP END PARALLEL DO


end program IofTobs
