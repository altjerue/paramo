program blazMag_main
   use data_types
   use misc
   implicit none
   character(len=*), parameter :: args_error = "Usage:"//new_line('A')//&
      '  xBlazMasg params_file KN-cool SSA-boiler with-wind flow-geom blob output-name'//new_line('A')//&
      'Options:'//new_line('A')//&
      '  params_file     Parameters file'//new_line('A')//&
      '  output-name     Name of the output file'//new_line('A')//&
      '  KN-cool         Klein-Nishina cooling: T/F (True/False)'//new_line('A')//&
      '  with_SSA        with self-absorption: T/F (True/False)'//new_line('A')
   integer :: numArgs, i
   character(len=256) :: program_name
   character(len=256), allocatable, dimension(:) :: args
   logical, allocatable, dimension(:) :: with_arg

   numArgs = command_argument_count()
   if ( numArgs /= 4 ) call an_error(args_error)
   call get_command_argument(0, program_name)
   allocate(args(numArgs), with_arg(numArgs - 2))
   do i=1, numArgs
      call get_command_argument(i, args(i))
      if ( i > 2 ) then
         if ( args(i) == 'T' ) then
            with_arg(i) = .true.
         else if ( args(i) == 'F' ) then
            with_arg(i) = .false.
         else
            call an_error(args_error)
         end if
      end if
   end do

#ifdef INTERJETS
   call interJets(trim(args(1)), args(2))
#else
   call blazMag(trim(args(1)), args(2), with_arg(1), with_arg(2))
#endif
end program blazMag_main


! arguments:
!
!   - full_ssc_cooling
!   - ssa_boiler
