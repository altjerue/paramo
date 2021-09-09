program turBlaz_main
    use data_types
    use misc
    implicit none
    character(len=*), parameter :: args_error = "Usage:"//new_line('A')//&
       '  xTurBlaz params_file KN-cool SSAbs with-wind flow-geom blob output-name'//new_line('A')//&
       'Options:'//new_line('A')//&
       '  params_file     Parameters file'//new_line('A')//&
       '  KN-cool         Klein-Nishina cooling: T/F (True/False)'//new_line('A')//&
       '  SSAbs           Synchrotron Self-Absorption: T/F (True/False)'//new_line('A')//&
       '  output-name     Name of the output file'//new_line('A')
    integer :: numArgs
    character(len=256) :: program_name, params_file, output_file, wCool, wAbs
    logical :: with_cool, with_abs
 
    numArgs = command_argument_count()
    call get_command_argument(0, program_name)
 
    if ( numArgs /= 4 ) call an_error(args_error)
 
    call get_command_argument(1, params_file)
    call get_command_argument(2, wCool)
    call get_command_argument(3, wAbs)
    call get_command_argument(4, output_file)
 
    !    -----> With or without absorption
    if ( wAbs == 'True' ) then
       with_abs = .true.
    elseif ( wAbs == 'False' ) then
       with_abs = .false.
    else
       call an_error(args_error)
    end if
 
    !    -----> With or without SSC emissivity
    if ( wCool == 'True' ) then
       with_cool = .true.
    elseif ( wCool == 'False' ) then
       with_cool = .false.
    else
       call an_error(args_error)
    end if
 
    call turBlaz(trim(params_file), output_file, with_cool, with_abs)
 
 end program turBlaz_main
