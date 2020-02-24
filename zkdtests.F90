!some sudo code
!gram looks like a name to call the program by so make this.
  !--looked it up still don't know time to ask jesus
gram testszkd
!add all dependencies
use hdf5
!not sure what implicit is
  !---implicit none is to keep i,j,k,l,m,n to be assumed to be integeres
  !--- implicit none should always be used
implicit none

!add subroutine to do are calculation:
call steadysh_state
!write finised ie learn how write function works
write(*,*) '===== Finished ====== '
!create the functrion you called earlier to do the calculation
contains
!contains is some way of organizing subroutines ask jesus
   subroutine steadyish_state
     implicit none
     write(*,*) "compile tests"
!again not sure what implicit is
!reserve memory for all needed variable and outputs? so learn how variable are made in fortran
  !--outputs are managed by using the intent(out)
  !--intent doesn't appear to be needed if it is not a parameter or variable that leaves
!then initiallize all inputs taht are scalar
!next initiallize all matrix and vectors
!learn how for loops and if statements work in fortran
  !--looks like you can name for loops adn if you do you end with end do NAME
!end function and in theory feel proud
  end subroutine steadysh_state

end program testszkd
