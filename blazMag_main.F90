program blazMag_main
   use data_types
   use misc
   implicit none

   character(len=*), parameter :: args_error = "Usage:"//new_line('A')//&
      '  xParamo params_file with-SSA with-IC with-cooling'//new_line('A')//&
      'Options:'//new_line('A')//&
      '  params_file: Parameters file'//new_line('A')//&
      '  with MBS self-absorption: T/F (True/False)'//new_line('A')//&
      '  with IC: T/F (True/False)'//new_line('A')//&
      '  with cooling: T/F (True/False)'

   integer :: numArgs
   character(len=256) :: program_name, params_file, output_file, wCool, &
      wAbs, wIC
   logical :: with_cool, with_abs, with_ic

   numArgs = command_argument_count()
   call get_command_argument(0, program_name)

   if ( numArgs /= 5 ) call an_error(args_error)

   !     ::::::::::     Opening parameters file     ::::::::::
   call get_command_argument(1, params_file)

   !    -----> With or without absorption
   call get_command_argument(2, wAbs)
   if ( wAbs == 'T' ) then
      with_abs = .true.
   elseif ( wAbs == 'F' ) then
      with_abs = .false.
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

   !    -----> With or without cooling
   call get_command_argument(4, wCool)
   if ( wCool == 'T' ) then
      with_cool = .true.
   elseif ( wCool == 'F' ) then
      with_cool = .false.
   else
      call an_error(args_error)
   end if

   !    -----> Reading output file
   call get_command_argument(5, output_file)

   call blazMag(trim(params_file), trim(output_file), with_cool, with_abs, with_ic)

end program blazMag_main
