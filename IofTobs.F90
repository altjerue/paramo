program IofTobs
   use data_types
   use constants
   use misc
   use hdf5
   use h5_inout
   !$ use omp_lib
   use SRtoolkit
   use pwl_integ
   implicit none

   integer :: herror, numArgs
   integer :: i, j, ii, numdf, numdt, i_edge, i_start
   integer(HID_T) :: file_id, group_id

   real(dp) :: B, R, d_L, z, Gbulk, theta, mu_obs, mu_com, D, tob_min, &
      tob_max, sind
   real(dp), allocatable, dimension(:) :: t, t_obs, nu, s
   real(dp), allocatable, dimension(:, :) :: jnut, Iobs

   character(len=256) :: paramo_fname
   character(len=*), parameter :: args_error = "Usage:"//new_line('A')//&
   "   xITobs part-evol-file"//new_line('A')//&
   "Options:"//new_line('A')//&
   "   part-evol-file: File in HDF5 format from Paramo"

   logical :: at_the_edge, Iobs_exists

   numArgs = command_argument_count()
   ! call get_command_argument(0, program_name)

   if (numArgs /= 1) call an_error(args_error)
   call get_command_argument(1, paramo_fname)

   ! Opening Paramo file
   call h5open_f(herror)
   call h5io_openf(trim(paramo_fname), file_id, herror)
   call h5io_openg(file_id, 'Params', group_id, herror)
   call h5io_rint0(group_id, 'numdf', numdf, herror)
   call h5io_rint0(group_id, 'numdt', numdt, herror)
   call h5io_rdble0(group_id, 'Bfield', B, herror)
   call h5io_rdble0(group_id, 'R', R, herror)
   call h5io_rdble0(group_id, 'd_lum', d_L, herror)
   call h5io_rdble0(group_id, 'redshift', z, herror)
   call h5io_rdble0(group_id, 'gamma_bulk', Gbulk, herror)
   call h5io_rdble0(group_id, 'view-angle', theta, herror)
   call h5io_closeg(group_id, herror)

   allocate(t(numdt), nu(numdf), s(numdt), t_obs(numdt))
   allocate(jnut(numdf, numdt), Iobs(numdf, numdt))

   ! Reading data
   call h5io_rdble1(file_id, 'time', t, herror)
   call h5io_rdble1(file_id, 't_obs', t_obs, herror)
   call h5io_rdble1(file_id, 'sen_lum', s, herror)
   call h5io_rdble1(file_id, 'nu_obs', nu, herror)
   call h5io_rdble2(file_id, 'jnut', jnut, herror)


   ! ----->>   Viewing angle
   mu_obs = dcos(theta * pi / 180d0)
   mu_com = mu_com_f(Gbulk, mu_obs)
   D = Doppler(Gbulk, mu_obs)

   write(*, *) ''
   write(*, *) '--> Moving to the observer frame'

   Iobs = 0d0
   i_edge = minloc(abs(2d0 * R * mu_com - s), dim=1, mask=2d0 * R * mu_com < s)
   if ( i_edge == 0 ) i_edge = numdt + 1

   !$OMP PARALLEL DO ORDERED COLLAPSE(2) SCHEDULE(AUTO) DEFAULT(SHARED) &
   !$OMP& PRIVATE(tob_min, tob_max, ii, sind, i_start, at_the_edge)
   freq_loop: do j=1, numdf
      obs_loop: do i = 1, numdt
         at_the_edge = .false.
         if ( 2d0 * R * mu_com - s(i) < 0d0 ) at_the_edge = .true.
         if ( .not. at_the_edge ) then
            i_start = 1
         else
            i_start = i - i_edge
         end if
         int_loop: do ii = i_start, i
            if ( ii == 1 ) then
               tob_min = t_com_f(t_obs(i), z, Gbulk, 0d0, mu_obs)
               tob_max = t_com_f(t_obs(i), z, Gbulk, t(ii) * cLight * mu_com, mu_obs)
               Iobs(j, i) = dabs(tob_max - tob_min) * jnut(j, 1)
            else
               tob_min = t_com_f(t_obs(i), z, Gbulk, t(ii - 1) * cLight * mu_com, mu_obs)
               tob_max = t_com_f(t_obs(i), z, Gbulk, t(ii    ) * cLight * mu_com, mu_obs)
               if ( jnut(j, ii) > 1d-100 .and. jnut(j, ii - 1) > 1d-100 ) then
                  sind = -dlog(jnut(j, ii) / jnut(j, ii - 1)) / dlog(tob_max / tob_min)
                  if ( sind < -8d0 ) sind = -8d0
                  if ( sind > 8d0 ) sind = 8d0
                  Iobs(j, i) = Iobs(j, i) + jnut(j, ii - 1) * tob_min * Pinteg(tob_max / tob_min, sind, 1d-6) / (Gbulk * mu_com * (mu_obs - bofg(Gbulk)) * D)
               end if
            end if
         end do int_loop
      end do obs_loop
   end do freq_loop
   !$OMP END PARALLEL DO

   call h5lexists_f(file_id, 'Iobs', Iobs_exists, herror)
   if ( Iobs_exists ) call h5ldelete_f(file_id, 'Iobs', herror)
   call h5io_wdble2(file_id, 'Iobs', Iobs, herror)

   ! Closing everything
   call h5io_closef(file_id, herror)
   call h5close_f(herror)

   if ( Iobs_exists ) then
      call system('h5repack '//trim(paramo_fname)//' tmp.h5')
      call system('mv tmp.h5 '//trim(paramo_fname))
   end if

end program IofTobs
