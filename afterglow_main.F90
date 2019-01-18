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
   character(len=256) :: program_name, params_file, output_file, wCool
   logical :: with_cool

   numArgs = command_argument_count()
   call get_command_argument(0, program_name)

   if ( numArgs /= 3 ) call an_error(args_error)

   !    -----> With or without cooling
   call get_command_argument(1, wCool)
   if ( wCool == 'T' ) then
      with_cool = .true.
   elseif ( wCool == 'F' ) then
      with_cool = .false.
   else
      call an_error(args_error)
   end if

   call get_command_argument(2, params_file)
   call get_command_argument(3, output_file)

   call afterglow(trim(params_file), trim(output_file), with_cool)

end program afterglow_main
