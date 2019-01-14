program IofTobs
   use data_types
   use constants
   use misc
   use hdf5
   use h5_inout
   !$ use omp_lib
   use SRtoolkit
   use pwl_integ
   use radiation
   implicit none

   integer(HID_T) :: file_id, group_id
   integer :: herror, numArgs, j, i, ii, numdf, numdt, i_edge, i_start
   real(dp) :: B, R, d_L, z, Gbulk, theta, mu_obs, mu_com, D, abu, volume, tau,&
      s_min, s_max
   real(dp), allocatable, dimension(:) :: t, t_obs, nu, s, pos
   real(dp), allocatable, dimension(:, :) :: jmbs, jssc, jeic, ambs, flux

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
   write(*,*) '=======  nu F_nu(t_obs)  ======='
   write(*, *) '-> Working with file: '//trim(paramo_fname)

   !   ----->   Opening Paramo file
   call h5open_f(herror)
   call h5io_openf(paramo_fname, file_id, herror)
   call h5io_openg(file_id, 'Parameters', group_id, herror)
   call h5io_rint0(group_id, 'numdf', numdf, herror)
   call h5io_rdble0(group_id, 'R', R, herror)
   call h5io_rdble0(group_id, 'd_lum', d_L, herror)
   call h5io_rdble0(group_id, 'redshift', z, herror)
   call h5io_rdble0(group_id, 'gamma_bulk', Gbulk, herror)
   call h5io_rdble0(group_id, 'view-angle', theta, herror)
   call h5io_closeg(group_id, herror)
   call h5io_rint0(file_id, 't_stop', numdt, herror)
   call h5io_rdble0(file_id, 'Bfield', B, herror)

   allocate(t(numdt), nu(numdf), s(numdt), t_obs(numdt), pos(numdt))
   allocate(jmbs(numdf, numdt), ambs(numdf, numdt), flux(numdf, numdt), &
      jssc(numdf, numdt), jeic(numdf, numdt))

   !   ----->   Reading data
   call h5io_rdble1(file_id, 'time', t, herror)
   call h5io_rdble1(file_id, 't_obs', t_obs, herror)
   call h5io_rdble1(file_id, 'sen_lum', s, herror)
   call h5io_rdble1(file_id, 'frequency', nu, herror)
   call h5io_rdble2(file_id, 'jmbs', jmbs, herror)
   call h5io_rdble2(file_id, 'jssc', jssc, herror)
   call h5io_rdble2(file_id, 'jeic', jeic, herror)
   call h5io_rdble2(file_id, 'anut', ambs, herror)


   mu_obs = dcos(theta * pi / 180d0)
   mu_com = mu_com_f(Gbulk, mu_obs)
   D = Doppler(Gbulk, mu_obs)
   volume = 4d0 * pi * R**3 / 3d0

   write(*, *) '-> Solving Radiative Transfer Equation'
   do i = 1, numdt
      do j = 1, numdf
         tau = dmax1(1d-200, 2d0 * R * ambs(j, i))
         flux(j, i) = D**4 * volume * ( 3d0 * nu(j) * opt_depth_blob(tau) * jmbs(j, i) / tau + nu(j) * (jssc(j, i) + jeic(j, i)) ) / (4d0 * pi * d_L**2)
      end do
   end do

   write(*, *) '-> Saving'
   call h5lexists_f(file_id, 'flux', Iobs_exists, herror)
   if ( Iobs_exists ) call h5ldelete_f(file_id, 'flux', herror)
   call h5io_wdble2(file_id, 'flux', flux, herror)

   ! Closing everything
   call h5io_closef(file_id, herror)
   call h5close_f(herror)

   if ( Iobs_exists ) then
      call system('h5repack '//trim(paramo_fname)//' tmp.h5')
      call system('mv tmp.h5 '//trim(paramo_fname))
   end if

end program IofTobs
