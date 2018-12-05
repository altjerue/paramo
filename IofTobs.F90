program IofTobs
   use data_types
   use constants
#ifdef FBAR
   use, intrinsic :: iso_fortran_env, only : I4P => int32, R8P => real64
   use forbear, only : bar_object
#endif
   use misc
   use hdf5
   use h5_inout
   !$ use omp_lib
   use SRtoolkit
   use pwl_integ
   use radiation
   implicit none

#ifdef FBAR
   type(bar_object) :: bar
   integer(I4P), volatile :: cur
#endif
   integer(HID_T) :: file_id, group_id
   integer :: herror, numArgs, i, ii, numdf, numdt, i_edge, i_start
   real(dp) :: B, R, d_L, z, Gbulk, theta, mu_obs, mu_com, D, factor, &
      abu, s_min, s_max
   real(dp), allocatable, dimension(:) :: t, t_obs, nu, s, pos
   real(dp), allocatable, dimension(:, :) :: jnut, anut, Iobs

   character(len=256) :: paramo_fname
   character(len=*), parameter :: args_error = "Usage:"//new_line('A')//&
   "   xITobs part-evol-file"//new_line('A')//&
   "Options:"//new_line('A')//&
   "   part-evol-file: File in HDF5 format from Paramo"

   logical :: Iobs_exists!, at_the_edge

   numArgs = command_argument_count()
   ! call get_command_argument(0, program_name)

   if (numArgs /= 1) call an_error(args_error)
   call get_command_argument(1, paramo_fname)

   write(*, *) ''
   write(*,*) '=======  I_nu(t_obs)  ======='
   write(*, *) '-> Working with file: '//trim(paramo_fname)

   !   ----->   Opening Paramo file
   call h5open_f(herror)
   call h5io_openf(paramo_fname, file_id, herror)
   call h5io_openg(file_id, 'Parameters', group_id, herror)
   call h5io_rint0(group_id, 'numdf', numdf, herror)
   call h5io_rdble0(group_id, 'Bfield', B, herror)
   call h5io_rdble0(group_id, 'R', R, herror)
   call h5io_rdble0(group_id, 'd_lum', d_L, herror)
   call h5io_rdble0(group_id, 'redshift', z, herror)
   call h5io_rdble0(group_id, 'gamma_bulk', Gbulk, herror)
   call h5io_rdble0(group_id, 'view-angle', theta, herror)
   call h5io_closeg(group_id, herror)
   call h5io_rint0(file_id, 't_stop', numdt, herror)

   allocate(t(numdt), nu(numdf), s(numdt), t_obs(numdt), pos(numdt))
   allocate(jnut(numdf, numdt), anut(numdf, numdt), Iobs(numdf, numdt))

   !   ----->   Reading data
   call h5io_rdble1(file_id, 'time', t, herror)
   call h5io_rdble1(file_id, 't_obs', t_obs, herror)
   call h5io_rdble1(file_id, 'sen_lum', s, herror)
   call h5io_rdble1(file_id, 'nu_obs', nu, herror)
   call h5io_rdble2(file_id, 'jnut', jnut, herror)
   call h5io_rdble2(file_id, 'anut', anut, herror)


   mu_obs = dcos(theta * pi / 180d0)
   mu_com = mu_com_f(Gbulk, mu_obs)
   D = Doppler(Gbulk, mu_obs)

   ! --> Light path from origin to the observer: 2 s mu_com
   ! --> Edge of the blob: 2 R mu_com
   ! --> Position at which we will measure the radiation:
   pos = R - s
   factor = 1d0 / (2d0 * Gbulk * mu_com * (mu_obs - bofg(Gbulk)) * D) ! <--- ds' / (c dt')
   i_edge = minloc(2d0 * R - s, dim = 1, mask = 2d0 * R - s >= 0d0)
   pos = 2d0 * mu_com * pos
   write(*, *) '-> Solving Radiative Transfer Equation'

#ifdef FBAR
   cur = 0
   call bar%initialize(width=48, max_value=real(numdt, R8P),                          &
      bracket_left_string='|', bracket_left_color_fg='blue',                          &
      empty_char_string=' ', empty_char_color_fg='blue', empty_char_color_bg='white', &
      filled_char_string=' ', filled_char_color_bg='blue',                            &
      bracket_right_string='|', bracket_right_color_fg='blue',                        &
      prefix_string='progress ', prefix_color_fg='red',                               &
      add_progress_percent=.true., progress_percent_color_fg='yellow',                &
      add_progress_speed=.false., progress_speed_color_fg='green',                    &
      add_date_time=.true., date_time_color_fg='magenta',                             &
      add_scale_bar=.false., scale_bar_color_fg='blue', scale_bar_style='underline_on')
   call bar%start
#endif

   !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(AUTO) DEFAULT(SHARED)&
   !$OMP& PRIVATE(s_min, s_max, abu, i, j, ii, sind, i_start)
   tobs_loop: do i = 1, numdt
      if ( i <= i_edge ) then
         i_start = 1
      else
         i_start = i - i_edge
      end if
      tcom_loop: do ii = i_start, i
         if ( ii == 1 ) then
            s_min = x_com_f(0d0, t_obs(i), z, Gbulk, mu_obs) / (2d0 * R * mu_com)
            s_max = x_com_f(t(1), t_obs(i), z, Gbulk, mu_obs) / (2d0 * R * mu_com)
            abu = s_max - s_min
            call RadTrans(Iobs(:, i), abu, jnu=jnut(:, 1), anu=anut(:, 1))
         else
            s_min = x_com_f(t(ii - 1), t_obs(i), z, Gbulk, mu_obs) / (2d0 * R * mu_com)
            s_max = x_com_f(t(ii), t_obs(i), z, Gbulk, mu_obs) / (2d0 * R * mu_com)
            abu = s_max - s_min
            call RadTrans(Iobs(:, i), abu, jnu=jnut(:, ii), anu=anut(:, ii), I0=Iobs(:, i))
         end if
      end do tcom_loop
#ifdef FBAR
   !$OMP CRITICAL
   cur = cur + 1
   if ( cur < numdt ) call bar%update(current=real(cur, R8P))
   !$OMP END CRITICAL
#endif
   end do tobs_loop
   !$OMP END PARALLEL DO

   write(*, *) '-> Saving'
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
