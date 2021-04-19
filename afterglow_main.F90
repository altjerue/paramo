program afterglow_main
   use data_types
   use misc
   implicit none
   character(len=*), parameter :: args_error = "Usage:"//new_line('A')//&
      '  xAglow with-cooling params_file output_file'//new_line('A')//&
      'Options:'//new_line('A')//&
      '  params_file     Parameters file'//new_line('A')//&
      '  KN-cool         Klein-Nishina cooling: T/F (True/False)'//new_line('A')//&
      '  SSA heating     with self-absorption heating: T/F (True/False)'//new_line('A')//&
      '  with-wind       constant or wind-like external medium:'//new_line('A')//&
      '  flow-geom       isotropic or beamed outflow: T/F (True/False)'//new_line('A')//&
      '                    -1   isotropic (spherical) blast wave'//new_line('A')//&
      '                     0   beamed jet with half-opening angle 1 / Gamma'//new_line('A')//&
      '                     1   beamed jet with initial half-opening angle 0.2 and'//new_line('A')//&
      '                         no sideways expansion'//new_line('A')//&
      '                     2   beamed jet with initial half-opening angle 0.2 and'//new_line('A')//&
      '                         non-relativistic sideways expansion'//new_line('A')//&
      '                     3   beamed jet with initial half-opening angle 0.2 and'//new_line('A')//&
      '                         relativistic sideways expansion'//new_line('A')//&
      '  blob            is the emission region a blob?: T/F (True/False)'//new_line('A')//&
      '  output-name     Name of the output file'//new_line('A')
   integer :: numArgs, fgeom
   character(len=256) :: program_name, params_file, output_file, wCoolKN, &
         wAbs, wWind, flow_geom, wBlob
   logical :: with_coolKN, with_wind, with_abs, with_blob

   numArgs = command_argument_count()
   call get_command_argument(0, program_name)

   if ( numArgs /= 7 ) call an_error(args_error)

   call get_command_argument(1, params_file)
   call get_command_argument(2, wCoolKN)
   call get_command_argument(3, wAbs)
   call get_command_argument(4, wWind)
   call get_command_argument(5, flow_geom)
   call get_command_argument(6, wBlob)
   call get_command_argument(7, output_file)

   !----->   With or without KN cooling
   if ( wCoolKN == 'True' ) then
      with_coolKN = .true.
   elseif ( wCoolKN == 'False' ) then
      with_coolKN = .false.
   else
      call an_error(args_error)
   end if

   !----->   With or without SSA boiler
   if ( wAbs == 'True' ) then
      with_abs = .true.
   elseif ( wAbs == 'False' ) then
      with_abs = .false.
   else
      call an_error(args_error)
   end if

   !----->   With or without SSC emissivity
   if ( wWind == 'True' ) then
      with_wind = .true.
   elseif ( wWind == 'False' ) then
      with_wind = .false.
   else
      call an_error(args_error)
   end if

   !----->   With or without SSC emissivity
   read(flow_geom, *) fgeom
   if ( fgeom > 3 .or. fgeom < -1 ) call an_error(args_error)

   !----->   With or without SSC emissivity
   if ( wBlob == 'True' ) then
      with_blob = .true.
   elseif ( wBlob == 'False' ) then
      with_blob = .false.
   else
      call an_error(args_error)
   end if

#ifdef MEZCAL
   call mezcal(trim(params_file), output_file, with_coolKN)
#else
   call bw1D_afterglow(trim(params_file), output_file, with_coolKN, with_abs, &
         with_wind, fgeom, with_blob)
#endif
end program afterglow_main
