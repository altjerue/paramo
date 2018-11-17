module data_types
   use, intrinsic :: iso_fortran_env
   implicit none

   integer, parameter :: sp=kind(1.0)  ! double precision
   integer, parameter :: dp=kind(1.0d0)  ! double precision
   ! * Note: kind(X) gives the kind of number X. Such value is compiler and processor dependent.
   ! integer, parameter :: sp = REAL32
   ! integer, parameter :: dp = REAL64
   ! integer, parameter :: qp = REAL128
end module data_types
