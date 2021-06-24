#define SINGLE_PREC      (0)
#define DOUBLE_PREC      (1)
#define ISO_FORT_REAL32  (2)
#define ISO_FORT_REAL64  (3)
#define ISO_FORT_REAL128 (4)
#define PRECISION_TYPE (DOUBLE_PREC)
!> Set the double precision
module data_types
   use, intrinsic :: iso_fortran_env
   implicit none
#if (PRECISION_TYPE==DOUBLE_PREC)
   integer, parameter :: dp=kind(1.0d0)
#elif (PRECISION_TYPE==ISO_FORT_REAL32)
   integer, parameter :: dp = REAL32
#elif (PRECISION_TYPE==ISO_FORT_REAL64)
   integer, parameter :: dp = REAL64
#elif (PRECISION_TYPE==ISO_FORT_REAL128)
   integer, parameter :: dp = REAL128
#else
   integer, parameter :: dp=kind(1.0)
#endif
end module data_types
