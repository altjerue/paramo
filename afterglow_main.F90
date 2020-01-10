program afterglow_main
   use data_types
   use misc
   implicit none

   character(len=*), parameter :: args_error = "Usage:"//new_line('A')//&
      '  xAglow with-cooling params_file output_file'//new_line('A')//&
      'Options:'//new_line('A')//&
      '  params_file       Parameters file'//new_line('A')//&
      '  output_file       Name of output file'//new_line('A')//&
      '  with cooling      T/F (True/False)'//new_line('A')
   integer :: numArgs
   character(len=256) :: program_name, params_file, output_file, wCool, wIC, wAbs
   logical :: with_cool, with_ic, with_abs

   numArgs = command_argument_count()
   call get_command_argument(0, program_name)

   if ( numArgs /= 5 ) call an_error(args_error)

   !    -----> With or without absorption
   call get_command_argument(1, wAbs)
   if ( wAbs == 'T' ) then
      with_abs = .true.
   elseif ( wAbs == 'F' ) then
      with_abs = .false.
   else
      call an_error(args_error)
   end if

   !    -----> With or without cooling
   call get_command_argument(2, wCool)
   if ( wCool == 'T' ) then
      with_cool = .true.
   elseif ( wCool == 'F' ) then
      with_cool = .false.
   else
      call an_error(args_error)
   end if

   !    -----> With or without SSC emissivity
   call get_command_argument(3, wIC)
   if ( wIC == 'T' ) then
      with_ic = .true.
   elseif ( wIC == 'F' ) then
      with_ic = .false.
   else
      call an_error(args_error)
   end if

   call get_command_argument(4, params_file)
   call get_command_argument(5, output_file)

   call afterglow(trim(params_file), trim(output_file), with_abs, with_cool, with_ic)

end program afterglow_main
