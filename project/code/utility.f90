! File: utility.f90
! Author: Lam tran
! Purpose: Define useful constants

module utility
  
  implicit none
  
  integer, parameter :: fp = selected_real_kind(15)
  integer, parameter :: maxFileLen=50, maxStrLen=200
  real (fp), parameter :: pi = acos(-1.0_fp)
  
end module utility
