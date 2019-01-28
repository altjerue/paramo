program blazMag_main
   use data_types
   use misc
   implicit none

   character(len=*), parameter :: args_error = "Usage:"//new_line('A')//&
      '  xParamo params_file with-cooling'//new_line('A')//&
      'Options:'//new_line('A')//&
      '  params_file: Parameters file'//new_line('A')//&
      '  with cooling: T/F (True/False)'//new_line('A')//&
      '  with SSC: T/F (True/False)'//new_line('A')//&
      '  with MBS self-absorption: T/F (True/False)'

   integer :: numArgs
   character(len=256) :: program_name, params_file, output_file, wCool, &
      wMBSabs, wSSCem, HybDis
   logical :: with_cool, with_abs, with_ssc, hyb_dis

   numArgs = command_argument_count()
   call get_command_argument(0, program_name)

   if ( numArgs /= 6 ) call an_error(args_error)

   !   ::::::::::   Opening parameters file   ::::::::::
   call get_command_argument(1, params_file)

   !    -----> With or without cooling
   call get_command_argument(2, HybDis)
   if ( HybDis == 'T' ) then
      hyb_dis = .true.
   elseif ( HybDis == 'F' ) then
      hyb_dis = .false.
   else
      call an_error(args_error)
   end if

   !    -----> With or without absorption
   call get_command_argument(3, wMBSabs)
   if ( wMBSabs == 'T' ) then
      with_abs = .true.
   elseif ( wMBSabs == 'F' ) then
      with_abs = .false.
   else
      call an_error(args_error)
   end if

   !    -----> With or without SSC emissivity
   call get_command_argument(4, wSSCem)
   if ( wSSCem == 'T' ) then
      with_ssc = .true.
   elseif ( wSSCem == 'F' ) then
      with_ssc = .false.
   else
      call an_error(args_error)
   end if

   !    -----> With or without cooling
   call get_command_argument(5, wCool)
   if ( wCool == 'T' ) then
      with_cool = .true.
   elseif ( wCool == 'F' ) then
      with_cool = .false.
   else
      call an_error(args_error)
   end if

   !    -----> Reading output file
   call get_command_argument(6, output_file)

   call blazMag(trim(params_file), trim(output_file), hyb_dis, with_cool, with_abs, with_ssc)

end program blazMag_main
