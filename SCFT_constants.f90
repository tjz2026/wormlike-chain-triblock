!-------------------------------------------------------------------------

! type    : module

! purpose : define the numerical and physical constants 
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

  module constants
     implicit none

!=========================================================================
!>>> integer constants: numerical precision                            <<<
!=========================================================================
!>>> integer constants: file unit handler                              <<<
!=========================================================================

! standard console output
     integer, public, parameter :: mystd = 6

! log file out
     integer, public, parameter :: myout = 99

! common file output
     integer, public, parameter :: mytmp = 100
end module constants
